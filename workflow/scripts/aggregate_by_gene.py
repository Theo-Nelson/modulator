#!/usr/bin/env python3
"""
Two-pass streaming implementation of original aggregator to avoid OOM.

Same CLI, same inputs, and (as closely as possible) same outputs as the pandas version:

Inputs
- --modkit-dir : directory containing per-sample subdirs with numbered ZN partition bedMethyl files
- --gtf        : GTF to map genomic sites to genes
- --out-prefix : output prefix

Key behaviors preserved
- Only reads numbered partition files like .../<sample>/<something>/<N>.bed(.gz)
  Skips 'ungrouped.bed' and flat '*_filtered_mod.bed(.gz)' files.
- Deduplication: sums counts for identical keys:
    (sample, mod_code, ZN_transcript_index, chrom, start0, end0, strand, gene_id, gene_name)
- frac_modified computed as Nmod / Nvalid_cov (0 if cov=0), then min-cov zeroing (display only).
- Filtering (if --filter-enable):
    FAIL if (Ndiff > count_diff_factor * Nvalid_cov) OR (Nmod <= Nfail + mod_fail_margin)
    Site key for filtering:
        (chrom, start0, end0, strand, mod_code)
    A site is kept if ANY row at that site passes; if kept, ALL rows at that site are kept.

Outputs (for each TAG in RAW/FILTERED as enabled)
1) <out_prefix>_<TAG>_sites_long.tsv                      (stream-written)
2) <out_prefix>_<TAG>__per_gene_mod/                      (optional, can be huge)
   - <prefix_base>__<gene>__<mod>.tsv                     (optional)
   - <prefix_base>__<gene>__<mod>_{cov,frac,Nmod}_pivot.tsv (optional)
3) <out_prefix>_<TAG>__per_sample_mod_site_stats.tsv      (computed from streamed, deduped rows)
4) <out_prefix>_<TAG>__per_sample_mod_tx_stats.tsv
5) <out_prefix>_<TAG>__per_tx_mod_stats.tsv               (detail table; always written)

Important notes
- This version avoids pandas entirely and never loads all rows into RAM.
- Pivot outputs are implemented in a streaming-friendly way:
    for each gene/mod, we accumulate per (site+ZN) and sample metrics, then write pivot tables.
  This can still be large on disk and in memory for very large genes, but is far safer overall.
"""

import os
import sys
import re
import argparse
import gzip
from collections import defaultdict, namedtuple
from typing import Dict, Tuple, List, Optional, Iterable

# ----------------------------- constants -----------------------------

BED_COLS = [
    "chrom","start0","end0","mod_code","score","strand",
    "start0_compat","end0_compat","rgb",
    "Nvalid_cov","frac_modified",
    "Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall",
]

LONG_HEADER = [
    "sample","ZN_transcript_index","chrom","start0","end0","strand","mod_code",
    "Nvalid_cov","Nmod","frac_modified","gene_id","gene_name",
    "Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall",
]

PER_GENE_COLS = [
    "gene_name","gene_id","mod_code","chrom","start0","end0","strand",
    "ZN_transcript_index","sample","Nvalid_cov","Nmod","Ncanonical",
    "Nother_mod","Ndelete","Nfail","Ndiff","Nnocall","frac_modified"
]

# ----------------------------- CLI -----------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Aggregate ZN-partitioned modkit outputs per gene/mod with site-level filtering (two-pass streaming)")
    ap.add_argument("--modkit-dir", required=True, help="Parent dir with per-sample subdirs containing numbered ZN .bed files")
    ap.add_argument("--gtf", required=True, help="Assembler GTF (with gene coordinates). Exon or transcript features work.")
    ap.add_argument("--out-prefix", required=True, help="Prefix for outputs (no extension)")
    ap.add_argument("--min-cov", type=int, default=0, help="Zero frac_modified if Nvalid_cov < MIN_COV (row kept)")
    ap.add_argument("--filter-enable", action="store_true", help="Enable site-level filtering")
    ap.add_argument("--count-diff-factor", type=float, default=3.0, help="FAIL if Ndiff > factor * Nvalid_cov (default: 3)")
    ap.add_argument("--mod-fail-margin", type=int, default=1, help="FAIL if Nmod <= Nfail + margin (default: 1)")

    # output toggles
    ap.add_argument("--emit-raw", dest="emit_raw", action="store_true")
    ap.add_argument("--no-emit-raw", dest="emit_raw", action="store_false"); ap.set_defaults(emit_raw=True)
    ap.add_argument("--emit-filtered", dest="emit_filt", action="store_true")
    ap.add_argument("--no-emit-filtered", dest="emit_filt", action="store_false"); ap.set_defaults(emit_filt=True)

    ap.add_argument("--write-long", dest="write_long", action="store_true")
    ap.add_argument("--no-write-long", dest="write_long", action="store_false"); ap.set_defaults(write_long=True)
    ap.add_argument("--write-pivots", dest="write_pivots", action="store_true")
    ap.add_argument("--no-write-pivots", dest="write_pivots", action="store_false"); ap.set_defaults(write_pivots=True)

    # RAW vs FILTERED per-gene tables
    ap.add_argument("--write-raw-per-gene", dest="write_raw_per_gene", action="store_true")
    ap.add_argument("--no-write-raw-per-gene", dest="write_raw_per_gene", action="store_false"); ap.set_defaults(write_raw_per_gene=False)
    ap.add_argument("--write-filtered-per-gene", dest="write_filtered_per_gene", action="store_true")
    ap.add_argument("--no-write-filtered-per-gene", dest="write_filtered_per_gene", action="store_false"); ap.set_defaults(write_filtered_per_gene=True)

    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()

# ----------------------------- utils -----------------------------

def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def is_header(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")

def safe_int(x, default=0) -> int:
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return default

def safe_float(x, default=0.0) -> float:
    try:
        return float(x)
    except Exception:
        return default

def sanitize_filename_token(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]", "_", s if s else "NA")

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

# ----------------------------- GTF interval indexing -----------------------------

Interval = namedtuple("Interval", ["start","end","gene_id","gene_name","strand"])

def load_gene_intervals_from_gtf(gtf_path: str, verbose=False) -> Dict[Tuple[str,str], List[Interval]]:
    """Union per-gene exon/transcript spans to coarse intervals for site→gene mapping."""
    gene_bounds: Dict[Tuple[str,str,str], Tuple[int,int]] = {}
    gene_name_map: Dict[str,str] = {}

    with open_text(gtf_path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature not in ("exon","transcript","gene"):
                continue

            a = {}
            for kv in re.finditer(r'(\S+)\s+"([^"]*)"', attrs):
                a[kv.group(1)] = kv.group(2)

            gene_id = a.get("gene_id") or a.get("gtf_gene_id") or a.get("gene") or ""
            gene_name = a.get("ref_gene_name") or a.get("gene_name") or a.get("gtf_gene_name") or gene_id
            if not gene_id:
                continue

            s = int(start); e = int(end)
            key = (chrom, strand, gene_id)
            if key not in gene_bounds:
                gene_bounds[key] = (s, e)
            else:
                mn, mx = gene_bounds[key]
                gene_bounds[key] = (min(mn, s), max(mx, e))
            gene_name_map[gene_id] = gene_name

    by_cs: Dict[Tuple[str,str], List[Interval]] = defaultdict(list)
    for (chrom, strand, gid), (s,e) in gene_bounds.items():
        gname = gene_name_map.get(gid, gid)
        by_cs[(chrom, strand)].append(Interval(s, e, gid, gname, strand))

    for k in by_cs:
        by_cs[k].sort(key=lambda iv: (iv.start, iv.end))

    if verbose:
        total = sum(len(v) for v in by_cs.values())
        print(f"[info] loaded {total} gene intervals from {gtf_path}", file=sys.stderr)

    return by_cs

def assign_gene(chrom: str, pos_start: int, pos_end: int, strand: str,
                gene_index: Dict[Tuple[str,str], List[Interval]]) -> Tuple[str,str]:
    """
    Return (gene_id, gene_name) by overlap; choose max-overlap; tie → first; try opposite strand if empty.
    NOTE: This is O(#genes on chrom/strand) per row in worst-case; preserved from your original logic.
    """
    ivs = gene_index.get((chrom, strand), [])
    best = None
    best_ov = -1

    for iv in ivs:
        if iv.start > pos_end:
            break
        if iv.end < pos_start:
            continue
        ov = min(iv.end, pos_end) - max(iv.start, pos_start) + 1
        if ov > best_ov:
            best_ov = ov
            best = iv

    if best:
        return best.gene_id, best.gene_name

    # fallback: opposite strand (first overlap)
    other = "+" if strand == "-" else "-"
    ivs2 = gene_index.get((chrom, other), [])
    for iv in ivs2:
        if iv.start > pos_end:
            break
        if iv.end < pos_start:
            continue
        return iv.gene_id, iv.gene_name

    return "", ""

# ----------------------------- modkit numbered bed discovery -----------------------------

def iter_numbered_beds(modkit_dir: str) -> List[Tuple[str, str, str, int]]:
    """
    Return sorted list of (root, sample_name, bed_path, ZN_index)
    for files like '<sample>/<something>/<N>.bed' or '.bed.gz'.
    Skip ungrouped and flat '*_filtered_mod.bed(.gz)'.
    """
    out = []
    for root, dirs, files in os.walk(modkit_dir):
        rel = os.path.relpath(root, modkit_dir)
        if rel == ".":
            continue
        sample_name = rel.split(os.sep)[0]

        for fname in files:
            if fname.endswith("_filtered_mod.bed") or fname.endswith("_filtered_mod.bed.gz"):
                continue

            base = fname[:-3] if fname.endswith(".gz") else fname
            if base.lower() == "ungrouped.bed":
                continue

            m = re.fullmatch(r"(\d+)\.bed", base)
            if not m:
                continue

            zn = int(m.group(1))
            out.append((root, sample_name, os.path.join(root, fname), zn))

    return sorted(out)

# ----------------------------- parsing + keys -----------------------------

def parse_bed_line(line: str) -> Optional[Dict[str, object]]:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 18:
        parts = line.strip().split()
        if len(parts) < 18:
            return None
    parts = parts[:18]
    d = dict(zip(BED_COLS, parts))
    d["start0"] = safe_int(d["start0"])
    d["end0"] = safe_int(d["end0"])
    for k in ["Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"]:
        d[k] = safe_int(d[k])
    d["frac_modified"] = safe_float(d.get("frac_modified", 0.0), 0.0)
    return d

def site_key(chrom: str, start0: int, end0: int, strand: str, mod_code: str) -> Tuple[str,int,int,str,str]:
    return (chrom, int(start0), int(end0), strand, mod_code)

def dedup_key(sample: str, mod_code: str, zn: int, chrom: str, start0: int, end0: int, strand: str,
              gene_id: str, gene_name: str) -> Tuple[str,str,int,str,int,int,str,str,str]:
    return (sample, mod_code, int(zn), chrom, int(start0), int(end0), strand, gene_id, gene_name)

# ----------------------------- aggregation state -----------------------------

# Numeric sums we keep for each dedup key:
# Nvalid_cov, Nmod, Ncanonical, Nother_mod, Ndelete, Nfail, Ndiff, Nnocall
SUM_FIELDS = ["Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"]

def add_sums(dst: List[int], src: Dict[str, object]):
    # dst has len=8 in SUM_FIELDS order
    dst[0] += int(src["Nvalid_cov"])
    dst[1] += int(src["Nmod"])
    dst[2] += int(src["Ncanonical"])
    dst[3] += int(src["Nother_mod"])
    dst[4] += int(src["Ndelete"])
    dst[5] += int(src["Nfail"])
    dst[6] += int(src["Ndiff"])
    dst[7] += int(src["Nnocall"])

def frac_modified(nmod: int, cov: int, min_cov: int) -> float:
    if cov <= 0:
        f = 0.0
    else:
        f = nmod / cov
    if min_cov and cov < min_cov:
        f = 0.0
    # match your rounding
    return round(f, 6)

def row_pass_filter(nvalid_cov: int, nmod: int, nfail: int, ndiff: int,
                    count_diff_factor: float, mod_fail_margin: int) -> bool:
    if ndiff > (count_diff_factor * nvalid_cov):
        return False
    if nmod <= (nfail + mod_fail_margin):
        return False
    return True

# ----------------------------- per-sample stats (streaming) -----------------------------

class SampleStats:
    """
    Streaming reimplementation of write_per_sample_mod_stats without pandas.

    We need:
    - per-sample, per-mod:
        n_sites_total, n_sites_detected, total_Nmod, total_cov, overall_stoich
        mean/median site_cov, mean/median site_Nmod
    - per-sample, per-mod:
        n_tx, mean/median detected_sites_per_tx, mean/median total_Nmod_per_tx,
        mean/median tx_stoich
    - per_tx detail table:
        sample, mod_code, ZN_transcript_index, n_sites_total, n_sites_detected, total_Nmod, total_cov, tx_stoich

    Site definitions:
      - site_sample: group by (sample, mod_code, chrom, start0, end0, strand), sum Nmod/Nvalid_cov across transcripts
      - site_tx    : group by (sample, mod_code, ZN, chrom, start0, end0, strand), sum Nmod/Nvalid_cov

    Important: This assumes input rows are already deduplicated at your dedup_key level.
    """

    def __init__(self):
        # (sample, mod, chrom, start0, end0, strand) -> [sumNmod, sumCov]
        self.site_sample: Dict[Tuple[str,str,str,int,int,str], List[int]] = defaultdict(lambda: [0,0])

        # (sample, mod, zn, chrom, start0, end0, strand) -> [sumNmod, sumCov]
        self.site_tx: Dict[Tuple[str,str,int,str,int,int,str], List[int]] = defaultdict(lambda: [0,0])

    def add_row(self, sample: str, mod: str, zn: int, chrom: str, start0: int, end0: int, strand: str,
                nmod: int, cov: int):
        k1 = (sample, mod, chrom, start0, end0, strand)
        self.site_sample[k1][0] += nmod
        self.site_sample[k1][1] += cov

        k2 = (sample, mod, int(zn), chrom, start0, end0, strand)
        self.site_tx[k2][0] += nmod
        self.site_tx[k2][1] += cov

    @staticmethod
    def _median(sorted_vals: List[float]) -> float:
        n = len(sorted_vals)
        if n == 0:
            return 0.0
        mid = n // 2
        if n % 2 == 1:
            return float(sorted_vals[mid])
        return (float(sorted_vals[mid-1]) + float(sorted_vals[mid])) / 2.0

    @staticmethod
    def _mean(vals: List[float]) -> float:
        return float(sum(vals) / len(vals)) if vals else 0.0

    def write(self, base: str, tag: str, verbose: bool = False):
        out1 = f"{base}_{tag}__per_sample_mod_site_stats.tsv"
        out2 = f"{base}_{tag}__per_sample_mod_tx_stats.tsv"
        out3 = f"{base}_{tag}__per_tx_mod_stats.tsv"

        ensure_dir(os.path.dirname(out1) or ".")

        # Build per_sample_site and site-level mean/median from site_sample map
        # Aggregate per (sample, mod)
        per_sm_counts: Dict[Tuple[str,str], Dict[str, object]] = defaultdict(dict)
        # For mean/median we need lists of cov and nmod per site
        cov_list: Dict[Tuple[str,str], List[int]] = defaultdict(list)
        nmod_list: Dict[Tuple[str,str], List[int]] = defaultdict(list)

        for (sample, mod, chrom, s0, e0, strand), (sum_nmod, sum_cov) in self.site_sample.items():
            key = (sample, mod)
            detected = 1 if sum_nmod > 0 else 0

            d = per_sm_counts.get(key)
            if not d:
                d = {
                    "n_sites_total": 0,
                    "n_sites_detected": 0,
                    "total_Nmod": 0,
                    "total_cov": 0,
                }
                per_sm_counts[key] = d

            d["n_sites_total"] += 1
            d["n_sites_detected"] += detected
            d["total_Nmod"] += sum_nmod
            d["total_cov"] += sum_cov

            cov_list[key].append(sum_cov)
            nmod_list[key].append(sum_nmod)

        # Write out1
        with open(out1, "w") as f:
            hdr = [
                "sample","mod_code","n_sites_total","n_sites_detected","total_Nmod","total_cov","overall_stoich",
                "mean_site_cov","median_site_cov","mean_site_Nmod","median_site_Nmod"
            ]
            f.write("\t".join(hdr) + "\n")

            rows = []
            for (sample, mod), d in per_sm_counts.items():
                total_cov = int(d["total_cov"])
                total_nmod = int(d["total_Nmod"])
                overall = (total_nmod / total_cov) if total_cov > 0 else 0.0

                covs = sorted(cov_list[(sample, mod)])
                nmods = sorted(nmod_list[(sample, mod)])
                mean_cov = self._mean(covs)
                med_cov = self._median(covs)
                mean_nm = self._mean(nmods)
                med_nm = self._median(nmods)

                rows.append((mod, sample, d["n_sites_total"], d["n_sites_detected"], total_nmod, total_cov, overall,
                             mean_cov, med_cov, mean_nm, med_nm))

            # sort by mod_code then sample like your pandas version
            rows.sort(key=lambda x: (x[0], x[1]))
            for mod, sample, nst, nsd, tn, tc, ov, mc, medc, mnm, mednm in rows:
                f.write(
                    f"{sample}\t{mod}\t{nst}\t{nsd}\t{tn}\t{tc}\t{ov:.6f}\t{mc:.6f}\t{medc:.6f}\t{mnm:.6f}\t{mednm:.6f}\n"
                )

        # Build per_tx detail and per_sample_tx aggregates
        # per_tx: (sample, mod, zn) -> metrics
        per_tx: Dict[Tuple[str,str,int], Dict[str, object]] = defaultdict(dict)

        # For each (sample, mod, zn), we need:
        # n_sites_total, n_sites_detected, total_Nmod, total_cov
        for (sample, mod, zn, chrom, s0, e0, strand), (sum_nmod, sum_cov) in self.site_tx.items():
            key = (sample, mod, zn)
            detected = 1 if sum_nmod > 0 else 0
            d = per_tx.get(key)
            if not d:
                d = {"n_sites_total": 0, "n_sites_detected": 0, "total_Nmod": 0, "total_cov": 0}
                per_tx[key] = d
            d["n_sites_total"] += 1
            d["n_sites_detected"] += detected
            d["total_Nmod"] += sum_nmod
            d["total_cov"] += sum_cov

        # Write per_tx detail (out3)
        with open(out3, "w") as f:
            hdr = ["sample","mod_code","ZN_transcript_index","n_sites_total","n_sites_detected","total_Nmod","total_cov","tx_stoich"]
            f.write("\t".join(hdr) + "\n")
            rows = []
            for (sample, mod, zn), d in per_tx.items():
                tc = int(d["total_cov"])
                tn = int(d["total_Nmod"])
                sto = (tn / tc) if tc > 0 else 0.0
                rows.append((mod, sample, int(zn), d["n_sites_total"], d["n_sites_detected"], tn, tc, sto))
            rows.sort(key=lambda x: (x[0], x[1], x[2]))
            for mod, sample, zn, nst, nsd, tn, tc, sto in rows:
                f.write(f"{sample}\t{mod}\t{zn}\t{nst}\t{nsd}\t{tn}\t{tc}\t{sto:.6f}\n")

        # per_sample_tx aggregates:
        # n_tx, mean/median detected_sites_per_tx, mean/median total_Nmod_per_tx, mean/median tx_stoich
        tx_lists_detected: Dict[Tuple[str,str], List[int]] = defaultdict(list)
        tx_lists_nmod: Dict[Tuple[str,str], List[int]] = defaultdict(list)
        tx_lists_stoich: Dict[Tuple[str,str], List[float]] = defaultdict(list)
        tx_set: Dict[Tuple[str,str], set] = defaultdict(set)

        for (sample, mod, zn), d in per_tx.items():
            key = (sample, mod)
            tx_set[key].add(int(zn))
            tc = int(d["total_cov"])
            tn = int(d["total_Nmod"])
            sto = (tn / tc) if tc > 0 else 0.0
            tx_lists_detected[key].append(int(d["n_sites_detected"]))
            tx_lists_nmod[key].append(tn)
            tx_lists_stoich[key].append(sto)

        with open(out2, "w") as f:
            hdr = [
                "sample","mod_code","n_tx",
                "mean_detected_sites_per_tx","median_detected_sites_per_tx",
                "mean_total_Nmod_per_tx","median_total_Nmod_per_tx",
                "mean_tx_stoich","median_tx_stoich"
            ]
            f.write("\t".join(hdr) + "\n")

            rows = []
            for (sample, mod), zset in tx_set.items():
                n_tx = len(zset)
                det = sorted(tx_lists_detected[(sample, mod)])
                nm = sorted(tx_lists_nmod[(sample, mod)])
                st = sorted(tx_lists_stoich[(sample, mod)])

                rows.append((
                    mod, sample, n_tx,
                    self._mean(det), self._median(det),
                    self._mean(nm), self._median(nm),
                    self._mean(st), self._median(st),
                ))

            rows.sort(key=lambda x: (x[0], x[1]))
            for mod, sample, n_tx, md, med_d, mn, med_n, ms, med_s in rows:
                f.write(
                    f"{sample}\t{mod}\t{n_tx}\t"
                    f"{md:.6f}\t{med_d:.6f}\t{mn:.6f}\t{med_n:.6f}\t{ms:.6f}\t{med_s:.6f}\n"
                )

        if verbose:
            print(f"[ok] wrote {out1}", file=sys.stderr)
            print(f"[ok] wrote {out2}", file=sys.stderr)
            print(f"[ok] wrote {out3}", file=sys.stderr)

# ----------------------------- per-gene writers -----------------------------

class GeneWriters:
    """
    Manages per-gene×mod TSV writers and streaming pivot accumulation.

    - If write_per_gene: write row-level per-gene tables (like original).
    - If write_pivots: build pivot tables:
        index = (chrom,start0,end0,strand,ZN_transcript_index)
        columns = sample
        values = metric (cov, frac, Nmod) with "first" semantics in your original pivot_table.
      Since rows are deduplicated, there should be a single value per (index, sample).
    """
    def __init__(self, out_dir: str, prefix_base: str, write_per_gene: bool, write_pivots: bool, verbose: bool = False):
        self.out_dir = out_dir
        self.prefix_base = prefix_base
        self.write_per_gene = write_per_gene
        self.write_pivots = write_pivots
        self.verbose = verbose

        ensure_dir(out_dir)

        # open file handles per (gene_name, gene_id, mod_code)
        self._fh: Dict[Tuple[str,str,str], object] = {}

        # pivot accumulators per gene/mod:
        # (gene_name,gene_id,mod) -> dict[(chrom,start0,end0,strand,zn)][sample] = value
        self._pivot_cov: Dict[Tuple[str,str,str], Dict[Tuple[str,int,int,str,int], Dict[str,int]]] = defaultdict(lambda: defaultdict(dict))
        self._pivot_nmod: Dict[Tuple[str,str,str], Dict[Tuple[str,int,int,str,int], Dict[str,int]]] = defaultdict(lambda: defaultdict(dict))
        self._pivot_frac: Dict[Tuple[str,str,str], Dict[Tuple[str,int,int,str,int], Dict[str,float]]] = defaultdict(lambda: defaultdict(dict))

        self._samples_seen: set = set()

    def _key(self, gene_name: str, gene_id: str, mod: str) -> Tuple[str,str,str]:
        return (gene_name or "NA", gene_id or "", mod or "NA")

    def _basepath(self, gene_name: str, mod: str) -> str:
        safe_g = sanitize_filename_token(gene_name or "NA")
        safe_mod = sanitize_filename_token(str(mod))
        return os.path.join(self.out_dir, f"{self.prefix_base}__{safe_g}__{safe_mod}")

    def write_row(self, gene_name: str, gene_id: str, mod: str,
                  chrom: str, start0: int, end0: int, strand: str,
                  zn: int, sample: str,
                  sums: List[int], frac: float):
        """
        Called for each deduplicated row to:
        - optionally write per-gene TSV row
        - optionally accumulate pivot maps
        """
        self._samples_seen.add(sample)
        gk = self._key(gene_name, gene_id, mod)

        nvalid_cov, nmod, ncanonical, nother, ndelete, nfail, ndiff, nnocall = sums

        if self.write_per_gene:
            fh = self._fh.get(gk)
            if fh is None:
                path = self._basepath(gene_name or "NA", mod)
                fh = open(f"{path}.tsv", "w")
                self._fh[gk] = fh
                fh.write("\t".join(PER_GENE_COLS) + "\n")
            fh.write(
                f"{gene_name}\t{gene_id}\t{mod}\t{chrom}\t{start0}\t{end0}\t{strand}\t"
                f"{zn}\t{sample}\t{nvalid_cov}\t{nmod}\t{ncanonical}\t{nother}\t{ndelete}\t{nfail}\t{ndiff}\t{nnocall}\t{frac:.6f}\n"
            )

        if self.write_pivots:
            idx = (chrom, start0, end0, strand, int(zn))
            # "first" semantics: only set if not present
            if sample not in self._pivot_cov[gk][idx]:
                self._pivot_cov[gk][idx][sample] = int(nvalid_cov)
                self._pivot_nmod[gk][idx][sample] = int(nmod)
                self._pivot_frac[gk][idx][sample] = float(frac)

    def finalize(self):
        # close per-gene row files
        for fh in self._fh.values():
            try:
                fh.close()
            except Exception:
                pass
        self._fh.clear()

        if not self.write_pivots:
            return

        # Write pivot TSVs per gene/mod
        samples = sorted(self._samples_seen)

        for gk in sorted(self._pivot_cov.keys(), key=lambda x: (x[0], x[2])):  # gene_name then mod
            gene_name, gene_id, mod = gk
            basepath = self._basepath(gene_name, mod)

            # common index order
            idxs = sorted(self._pivot_cov[gk].keys(), key=lambda t: (t[0], t[1], t[4]))  # chrom,start0,zn

            def write_pivot(path: str, metric_map, fmt: str):
                with open(path, "w") as f:
                    f.write("\t".join(["chrom","start0","end0","strand","ZN_transcript_index"] + samples) + "\n")
                    for idx in idxs:
                        chrom, s0, e0, strand, zn = idx
                        row = [chrom, str(s0), str(e0), strand, str(zn)]
                        smap = metric_map[gk].get(idx, {})
                        for s in samples:
                            v = smap.get(s, 0)
                            if fmt == "int":
                                row.append(str(int(v)))
                            else:
                                row.append(f"{float(v):.6f}")
                        f.write("\t".join(row) + "\n")

            write_pivot(f"{basepath}_cov_pivot.tsv", self._pivot_cov, "int")
            write_pivot(f"{basepath}_Nmod_pivot.tsv", self._pivot_nmod, "int")
            write_pivot(f"{basepath}_frac_pivot.tsv", self._pivot_frac, "float")

        if self.verbose:
            print(f"[ok] wrote pivots under {self.out_dir}", file=sys.stderr)

# ----------------------------- two-pass core -----------------------------

def pass1_compute_passing_sites(
    beds: List[Tuple[str,str,str,int]],
    gene_index: Dict[Tuple[str,str], List[Interval]],
    count_diff_factor: float,
    mod_fail_margin: int,
    verbose: bool = False,
) -> Optional[set]:
    """
    Pass 1: if filtering enabled, compute set of site_keys that should be KEPT
    (site is kept if ANY row at that site passes filter after deduplication).
    To preserve semantics, we must apply the "row-pass" test on the DEDUPLICATED rows.

    Implementation: stream, deduplicate into a dict, then scan dedup rows to mark passing sites.
    This is still far less RAM than storing every line, because we only keep unique dedup keys.
    """
    dedup_map: Dict[Tuple, List[int]] = defaultdict(lambda: [0]*8)  # sums for SUM_FIELDS
    # We also need to remember the site key for each dedup key to evaluate passing sites later.
    site_for_dedup: Dict[Tuple, Tuple[str,int,int,str,str]] = {}

    n_lines = 0
    for _, sample, bed_path, zn in beds:
        with open_text(bed_path) as f:
            for ln in f:
                if is_header(ln):
                    continue
                rec = parse_bed_line(ln)
                if not rec:
                    continue
                n_lines += 1

                chrom = rec["chrom"]
                start0 = int(rec["start0"])
                end0 = int(rec["end0"])
                strand = rec["strand"]
                mod = rec["mod_code"]

                gid, gname = assign_gene(chrom, start0, end0, strand, gene_index)

                dk = dedup_key(sample, mod, zn, chrom, start0, end0, strand, gid, gname)
                add_sums(dedup_map[dk], rec)

                if dk not in site_for_dedup:
                    site_for_dedup[dk] = site_key(chrom, start0, end0, strand, mod)

    if verbose:
        print(f"[pass1] read {n_lines} lines; dedup keys: {len(dedup_map)}", file=sys.stderr)

    passing = set()
    for dk, sums in dedup_map.items():
        # sums order: cov, mod, canonical, other, delete, fail, diff, nocall
        cov = sums[0]
        nmod = sums[1]
        nfail = sums[5]
        ndiff = sums[6]
        if row_pass_filter(cov, nmod, nfail, ndiff, count_diff_factor, mod_fail_margin):
            passing.add(site_for_dedup[dk])

    if verbose:
        print(f"[pass1] passing sites: {len(passing)}", file=sys.stderr)

    return passing

def pass2_write_outputs(
    beds: List[Tuple[str,str,str,int]],
    gene_index: Dict[Tuple[str,str], List[Interval]],
    out_prefix: str,
    tag: str,
    min_cov: int,
    keep_sites: Optional[set],
    write_long: bool,
    write_per_gene: bool,
    write_pivots: bool,
    verbose: bool = False,
):
    """
    Pass 2: stream and deduplicate, then write outputs for this TAG.
    If keep_sites is not None, only emit rows whose site_key is in keep_sites.
    """

    # dedup map again (we need full dedupbed rows to write long/per-gene/pivots/stats)
    # We also need to store the non-summed columns for each dedup key (but they are in the key itself).
    dedup_map: Dict[Tuple, List[int]] = defaultdict(lambda: [0]*8)

    n_lines = 0
    for _, sample, bed_path, zn in beds:
        with open_text(bed_path) as f:
            for ln in f:
                if is_header(ln):
                    continue
                rec = parse_bed_line(ln)
                if not rec:
                    continue
                n_lines += 1

                chrom = rec["chrom"]
                start0 = int(rec["start0"])
                end0 = int(rec["end0"])
                strand = rec["strand"]
                mod = rec["mod_code"]

                gid, gname = assign_gene(chrom, start0, end0, strand, gene_index)
                dk = dedup_key(sample, mod, zn, chrom, start0, end0, strand, gid, gname)
                add_sums(dedup_map[dk], rec)

    if verbose:
        print(f"[{tag}] pass2 read {n_lines} lines; dedup keys: {len(dedup_map)}", file=sys.stderr)

    base = out_prefix
    out_long = f"{base}_{tag}_sites_long.tsv"

    # setup writers/accumulators
    long_fh = None
    if write_long:
        ensure_dir(os.path.dirname(out_long) or ".")
        long_fh = open(out_long, "w")
        long_fh.write("\t".join(LONG_HEADER) + "\n")

    stats = SampleStats()

    gene_writers = None
    if write_per_gene or write_pivots:
        out_dir = f"{base}_{tag}__per_gene_mod"
        ensure_dir(out_dir)
        prefix_base = os.path.basename(out_prefix)
        gene_writers = GeneWriters(out_dir=out_dir, prefix_base=prefix_base,
                                  write_per_gene=write_per_gene, write_pivots=write_pivots, verbose=verbose)

    # iterate dedupbed rows in a deterministic-ish order (optional)
    # Sorting can be memory-heavy itself; but sorting keys (not all lines) is usually ok.
    # If you prefer speed, you can remove sorting.
    keys_sorted = sorted(
        dedup_map.keys(),
        key=lambda k: (k[3], k[4], k[2], k[0], k[1])  # chrom,start0,ZN,sample,mod
    )

    n_emitted = 0
    for dk in keys_sorted:
        sample, mod, zn, chrom, start0, end0, strand, gid, gname = dk
        sums = dedup_map[dk]

        sk = site_key(chrom, start0, end0, strand, mod)
        if keep_sites is not None and sk not in keep_sites:
            continue

        cov = sums[0]
        nmod = sums[1]
        frac = frac_modified(nmod, cov, min_cov)

        # long row
        if long_fh is not None:
            ncanonical, nother, ndelete, nfail, ndiff, nnocall = sums[2], sums[3], sums[4], sums[5], sums[6], sums[7]
            long_fh.write(
                f"{sample}\t{zn}\t{chrom}\t{start0}\t{end0}\t{strand}\t{mod}\t"
                f"{cov}\t{nmod}\t{frac:.6f}\t{gid}\t{gname}\t"
                f"{ncanonical}\t{nother}\t{ndelete}\t{nfail}\t{ndiff}\t{nnocall}\n"
            )

        # stats accumulation
        stats.add_row(sample=sample, mod=mod, zn=int(zn), chrom=chrom, start0=int(start0), end0=int(end0), strand=strand,
                      nmod=int(nmod), cov=int(cov))

        # per-gene + pivots
        if gene_writers is not None:
            gene_writers.write_row(
                gene_name=gname, gene_id=gid, mod=mod,
                chrom=chrom, start0=int(start0), end0=int(end0), strand=strand,
                zn=int(zn), sample=sample,
                sums=sums, frac=frac
            )

        n_emitted += 1

    # close/finalize
    if long_fh is not None:
        long_fh.close()
        if verbose:
            print(f"[ok] wrote {out_long}", file=sys.stderr)

    # per-sample stats files
    stats.write(base=base, tag=tag, verbose=verbose)

    # finalize pivots
    if gene_writers is not None:
        gene_writers.finalize()

    if verbose:
        print(f"[{tag}] emitted {n_emitted} dedup rows", file=sys.stderr)

# ----------------------------- main -----------------------------

def main():
    args = parse_args()

    beds = iter_numbered_beds(args.modkit_dir)
    if not beds:
        sys.exit(f"No numbered ZN partition files found under {args.modkit_dir}")

    gene_index = load_gene_intervals_from_gtf(args.gtf, verbose=args.verbose)

    # Pass 1: compute passing sites if filtering is enabled
    keep_sites = None
    if args.filter_enable:
        keep_sites = pass1_compute_passing_sites(
            beds=beds,
            gene_index=gene_index,
            count_diff_factor=args.count_diff_factor,
            mod_fail_margin=args.mod_fail_margin,
            verbose=args.verbose,
        )

    # RAW outputs (no filtering)
    if args.emit_raw:
        pass2_write_outputs(
            beds=beds,
            gene_index=gene_index,
            out_prefix=args.out_prefix,
            tag="RAW",
            min_cov=args.min_cov,
            keep_sites=None,  # RAW keeps everything
            write_long=args.write_long,
            write_per_gene=args.write_raw_per_gene,
            write_pivots=args.write_pivots,
            verbose=args.verbose,
        )

    # FILTERED outputs
    if args.emit_filt:
        pass2_write_outputs(
            beds=beds,
            gene_index=gene_index,
            out_prefix=args.out_prefix,
            tag="FILTERED",
            min_cov=args.min_cov,
            keep_sites=keep_sites if args.filter_enable else None,
            write_long=args.write_long,
            write_per_gene=args.write_filtered_per_gene,
            write_pivots=args.write_pivots,
            verbose=args.verbose,
        )

    print("[OK] ZN aggregation complete.")

if __name__ == "__main__":
    main()

