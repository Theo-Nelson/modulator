#!/usr/bin/env python3
"""
External-sort + streaming-group implementation of aggregate_by_gene.py (OOM-safe).

This is a drop-in replacement for your current script:
- Same CLI flags (no required new args).
- Same outputs and filenames as your current implementation.

Why this works:
- We NEVER build a global in-memory dict of dedup keys.
- We stream → write a normalized TSV → GNU sort (disk-backed) → single-pass reduce.
- Filtering uses a sorted list of passing site keys + a second sorted view for streaming join.

High-level pipeline
1) Stream all numbered ZN bed files:
   - parse bedMethyl row (18 cols)
   - assign gene_id/gene_name from GTF intervals
   - emit normalized TSV (one line per input row)

2) External sort normalized TSV by the dedup key:
   (sample, mod_code, ZN, chrom, start0, end0, strand, gene_id, gene_name)

3) Reduce sorted file to deduplicated rows:
   - sum numeric fields per key
   - compute frac_modified (+ min_cov display rule)
   - write:
      a) <out_prefix>_RAW_sites_long.tsv    (if enabled)
      b) <out_prefix>__DEDUP_RAW.tsv        (internal temp for downstream steps)
      c) passing_sites.tsv                 (if filter enabled; internal temp)
   - (Stats + per-gene outputs are computed from dedup TSV via additional sort/reduce steps)

4) If filtering enabled:
   - sort+unique passing_sites.tsv by site_key (chrom,start0,end0,strand,mod_code)
   - external sort dedup rows by site_key (+ rest)
   - streaming merge-join to write:
      a) <out_prefix>_FILTERED_sites_long.tsv (if enabled)
      b) <out_prefix>__DEDUP_FILTERED.tsv     (internal temp)

5) For each TAG in {RAW, FILTERED} (as enabled):
   - compute per-sample site stats TSVs (with exact medians) using sort/reduce
   - compute per-tx stats TSVs (with exact medians) using sort/reduce
   - compute per-gene tables + pivots using sort/reduce (bounded memory per gene/mod group)

Notes / assumptions
- Requires GNU sort available on the cluster.
- Uses $TMPDIR if set, else /tmp, for sort temp files and intermediates.
- The per-gene pivot generation holds only ONE gene/mod group’s pivot maps in memory at a time.
  Some genes can still be huge; if you hit memory there, we can shard pivots further.

"""

import os, sys, re, argparse, gzip, subprocess, tempfile, shutil
from collections import defaultdict, namedtuple
from typing import Dict, Tuple, List, Optional, Iterable

# ----------------------------- Constants -----------------------------

BED_COLS = [
    "chrom","start0","end0","mod_code","score","strand",
    "start0_compat","end0_compat","rgb",
    "Nvalid_cov","frac_modified",
    "Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall",
]

# Normalized (pre-dedup) TSV columns we emit:
# NOTE: We do NOT emit "score/rgb/compat" because you don't use them downstream.
NORM_COLS = [
    "sample","ZN_transcript_index","chrom","start0","end0","strand","mod_code",
    "gene_id","gene_name",
    "Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall",
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

SUM_FIELDS = ["Nvalid_cov","Nmod","Ncanonical","Nother_mod","Ndelete","Nfail","Ndiff","Nnocall"]

# ----------------------------- CLI -----------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Aggregate ZN-partitioned modkit outputs per gene/mod with site-level filtering (external sort)")
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

    # optional tuning (NOT required; defaults are safe)
    ap.add_argument("--sort-parallel", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")),
                    help="GNU sort --parallel (default: SLURM_CPUS_PER_TASK or 4)")
    ap.add_argument("--sort-mem", type=str, default=os.environ.get("AGG_SORT_MEM", "2G"),
                    help="GNU sort -S memory limit (default: 2G or env AGG_SORT_MEM)")
    ap.add_argument("--tmpdir", type=str, default=os.environ.get("TMPDIR", "/tmp"),
                    help="Temp directory for sort/intermediates (default: $TMPDIR or /tmp)")
    ap.add_argument("--keep-intermediates", action="store_true", help="Do not delete intermediate TSVs (debugging)")

    return ap.parse_args()

# ----------------------------- Utils -----------------------------

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

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def sanitize_filename_token(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]", "_", s if s else "NA")

def frac_modified(nmod: int, cov: int, min_cov: int) -> float:
    if cov <= 0:
        f = 0.0
    else:
        f = nmod / cov
    if min_cov and cov < min_cov:
        f = 0.0
    return round(f, 6)

def row_pass_filter(cov: int, nmod: int, nfail: int, ndiff: int, count_diff_factor: float, mod_fail_margin: int) -> bool:
    if ndiff > (count_diff_factor * cov):
        return False
    if nmod <= (nfail + mod_fail_margin):
        return False
    return True

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

def assign_gene(chrom: str, pos_start: int, pos_end: int, strand: str, gene_index: Dict[Tuple[str,str], List[Interval]]) -> Tuple[str,str]:
    """Return (gene_id, gene_name) by overlap; choose max-overlap; tie → first; try opposite strand if empty."""
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

    other = "+" if strand == "-" else "-"
    ivs2 = gene_index.get((chrom, other), [])
    for iv in ivs2:
        if iv.start > pos_end:
            break
        if iv.end < pos_start:
            continue
        return iv.gene_id, iv.gene_name

    return "", ""

# ----------------------------- Bed discovery -----------------------------

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

# ----------------------------- bedMethyl parsing -----------------------------

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

# ----------------------------- GNU sort helpers -----------------------------

def run_sort(in_path: str, out_path: str, keys: List[Tuple[int,int,bool]], tmpdir: str, mem: str, parallel: int, unique: bool = False):
    """
    keys: list of (start_col_1based, end_col_1based, numeric)
    """
    # Ensure temp and output dirs exist
    ensure_dir(os.path.dirname(out_path) or ".")

    cmd = ["sort", "--parallel", str(max(1, parallel)), "-S", mem, "-T", tmpdir, "-t", "\t"]
    for (a, b, numeric) in keys:
        k = f"{a},{b}"
        cmd += ["-k", k + ("n" if numeric else "")]
    if unique:
        cmd.append("-u")

    # Enforce bytewise stable locale for speed/consistency
    env = os.environ.copy()
    env["LC_ALL"] = "C"

    with open(in_path, "r") as fin, open(out_path, "w") as fout:
        subprocess.check_call(cmd, stdin=fin, stdout=fout, env=env)

# ----------------------------- Stage 1: normalize -----------------------------

def normalize_to_tsv(beds: List[Tuple[str,str,str,int]], gene_index, out_tsv: str, verbose: bool = False) -> int:
    """
    Stream all numbered bed files and write normalized TSV with NORM_COLS.
    Returns number of rows written.
    """
    ensure_dir(os.path.dirname(out_tsv) or ".")
    n = 0
    with open(out_tsv, "w") as out:
        # no header (for sort speed); we know column order
        for _, sample, bed_path, zn in beds:
            with open_text(bed_path) as f:
                for ln in f:
                    if is_header(ln):
                        continue
                    rec = parse_bed_line(ln)
                    if not rec:
                        continue
                    chrom = rec["chrom"]
                    start0 = int(rec["start0"])
                    end0 = int(rec["end0"])
                    strand = rec["strand"]
                    mod = rec["mod_code"]

                    gid, gname = assign_gene(chrom, start0, end0, strand, gene_index)

                    vals = [
                        sample, str(int(zn)), chrom, str(start0), str(end0), strand, mod,
                        gid, gname,
                        str(int(rec["Nvalid_cov"])),
                        str(int(rec["Nmod"])),
                        str(int(rec["Ncanonical"])),
                        str(int(rec["Nother_mod"])),
                        str(int(rec["Ndelete"])),
                        str(int(rec["Nfail"])),
                        str(int(rec["Ndiff"])),
                        str(int(rec["Nnocall"])),
                    ]
                    out.write("\t".join(vals) + "\n")
                    n += 1
    if verbose:
        print(f"[norm] wrote {n} rows to {out_tsv}", file=sys.stderr)
    return n

# ----------------------------- Stage 2: dedup reduce -----------------------------

def dedup_reduce_sorted(
    sorted_norm_tsv: str,
    out_dedup_tsv: str,
    out_long_tsv: Optional[str],
    out_passing_sites_tsv: Optional[str],
    min_cov: int,
    filter_enable: bool,
    count_diff_factor: float,
    mod_fail_margin: int,
    verbose: bool = False,
) -> int:
    """
    Reduce a sorted normalized TSV to deduplicated rows.
    Writes:
      - out_dedup_tsv: dedup rows with frac_modified as last column (tab-delimited, no header)
      - out_long_tsv : if provided, long format with header
      - out_passing_sites_tsv: if provided and filter_enable, site_key rows for passing rows (no header)
    Returns number of dedup rows written.
    """
    ensure_dir(os.path.dirname(out_dedup_tsv) or ".")
    if out_long_tsv:
        ensure_dir(os.path.dirname(out_long_tsv) or ".")

    # Dedup TSV columns (no header):
    # sample, zn, chrom, start0, end0, strand, mod_code, gene_id, gene_name,
    # Nvalid_cov, Nmod, Ncanonical, Nother_mod, Ndelete, Nfail, Ndiff, Nnocall, frac_modified
    dedup_written = 0

    long_fh = None
    if out_long_tsv is not None:
        long_fh = open(out_long_tsv, "w")
        long_fh.write("\t".join(LONG_HEADER) + "\n")

    pass_fh = None
    if filter_enable and out_passing_sites_tsv is not None:
        ensure_dir(os.path.dirname(out_passing_sites_tsv) or ".")
        pass_fh = open(out_passing_sites_tsv, "w")

    def flush(curr_key, sums):
        nonlocal dedup_written
        if curr_key is None:
            return
        (sample, zn, chrom, start0, end0, strand, mod, gid, gname) = curr_key
        cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall = sums
        frac = frac_modified(nmod, cov, min_cov)

        # dedup line (no header)
        with open(out_dedup_tsv, "a") as d:
            d.write(
                "\t".join([
                    sample, zn, chrom, start0, end0, strand, mod, gid, gname,
                    str(cov), str(nmod), str(ncan), str(nother), str(ndel),
                    str(nfail), str(ndiff), str(nnocall),
                    f"{frac:.6f}"
                ]) + "\n"
            )

        if long_fh is not None:
            long_fh.write(
                f"{sample}\t{zn}\t{chrom}\t{start0}\t{end0}\t{strand}\t{mod}\t"
                f"{cov}\t{nmod}\t{frac:.6f}\t{gid}\t{gname}\t"
                f"{ncan}\t{nother}\t{ndel}\t{nfail}\t{ndiff}\t{nnocall}\n"
            )

        if pass_fh is not None:
            if row_pass_filter(cov, nmod, nfail, ndiff, count_diff_factor, mod_fail_margin):
                pass_fh.write(f"{chrom}\t{start0}\t{end0}\t{strand}\t{mod}\n")

        dedup_written += 1

    # IMPORTANT: we must open out_dedup_tsv in write mode once, not per flush.
    # We'll do that by truncating it first and then appending in flush().
    with open(out_dedup_tsv, "w"):
        pass

    curr_key = None
    sums = [0]*8

    with open(sorted_norm_tsv, "r") as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            parts = ln.split("\t")
            # NORM_COLS order:
            # 0 sample,1 zn,2 chrom,3 start0,4 end0,5 strand,6 mod,7 gid,8 gname,
            # 9 cov,10 nmod,11 ncan,12 nother,13 ndel,14 nfail,15 ndiff,16 nnocall
            key = (parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], parts[7], parts[8])
            vals = [safe_int(parts[9]), safe_int(parts[10]), safe_int(parts[11]), safe_int(parts[12]),
                    safe_int(parts[13]), safe_int(parts[14]), safe_int(parts[15]), safe_int(parts[16])]

            if curr_key is None:
                curr_key = key
                sums = vals
                continue

            if key == curr_key:
                for i in range(8):
                    sums[i] += vals[i]
            else:
                flush(curr_key, sums)
                curr_key = key
                sums = vals

    flush(curr_key, sums)

    if long_fh is not None:
        long_fh.close()
        if verbose:
            print(f"[dedup] wrote long: {out_long_tsv}", file=sys.stderr)

    if pass_fh is not None:
        pass_fh.close()
        if verbose:
            print(f"[dedup] wrote passing-sites: {out_passing_sites_tsv}", file=sys.stderr)

    if verbose:
        print(f"[dedup] wrote {dedup_written} dedup rows to {out_dedup_tsv}", file=sys.stderr)

    return dedup_written

# ----------------------------- Filtering join -----------------------------

def filter_dedup_by_passing_sites(
    dedup_by_site_sorted: str,
    passing_sites_sorted: str,
    out_dedup_filtered: str,
    out_long_filtered: Optional[str],
    verbose: bool = False,
) -> int:
    """
    Streaming merge-join:
    - dedup_by_site_sorted is sorted by (chrom,start0,end0,strand,mod_code, then rest)
    - passing_sites_sorted is sorted unique by (chrom,start0,end0,strand,mod_code)

    Writes filtered dedup TSV (no header) + optional long TSV with header.
    Returns number of filtered dedup rows.
    """
    ensure_dir(os.path.dirname(out_dedup_filtered) or ".")
    if out_long_filtered:
        ensure_dir(os.path.dirname(out_long_filtered) or ".")

    # Truncate filtered dedup
    with open(out_dedup_filtered, "w"):
        pass

    long_fh = None
    if out_long_filtered is not None:
        long_fh = open(out_long_filtered, "w")
        long_fh.write("\t".join(LONG_HEADER) + "\n")

    def site_of_dedup(parts: List[str]) -> Tuple[str,str,str,str,str]:
        # dedup columns:
        # 0 sample,1 zn,2 chrom,3 start0,4 end0,5 strand,6 mod,...
        return (parts[2], parts[3], parts[4], parts[5], parts[6])

    # Read passing sites iterator
    ps_f = open(passing_sites_sorted, "r")
    ps_line = ps_f.readline()
    ps_key = None
    if ps_line:
        p = ps_line.rstrip("\n").split("\t")
        ps_key = (p[0], p[1], p[2], p[3], p[4])

    kept = 0
    with open(dedup_by_site_sorted, "r") as d:
        for ln in d:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            parts = ln.split("\t")
            dk_site = site_of_dedup(parts)

            # advance passing-sites until >= dk_site
            while ps_key is not None and ps_key < dk_site:
                ps_line = ps_f.readline()
                if not ps_line:
                    ps_key = None
                    break
                p = ps_line.rstrip("\n").split("\t")
                ps_key = (p[0], p[1], p[2], p[3], p[4])

            if ps_key is None:
                break

            if dk_site == ps_key:
                # keep this row
                with open(out_dedup_filtered, "a") as out:
                    out.write("\t".join(parts) + "\n")

                if long_fh is not None:
                    # long expects: sample zn chrom start0 end0 strand mod cov nmod frac gid gname ...
                    sample, zn, chrom, start0, end0, strand, mod, gid, gname = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], parts[7], parts[8]
                    cov, nmod = parts[9], parts[10]
                    ncan, nother, ndel, nfail, ndiff, nnocall = parts[11], parts[12], parts[13], parts[14], parts[15], parts[16]
                    frac = parts[17]
                    long_fh.write(
                        f"{sample}\t{zn}\t{chrom}\t{start0}\t{end0}\t{strand}\t{mod}\t"
                        f"{cov}\t{nmod}\t{frac}\t{gid}\t{gname}\t"
                        f"{ncan}\t{nother}\t{ndel}\t{nfail}\t{ndiff}\t{nnocall}\n"
                    )
                kept += 1

    ps_f.close()
    if long_fh is not None:
        long_fh.close()

    if verbose:
        print(f"[filter] kept {kept} dedup rows into {out_dedup_filtered}", file=sys.stderr)
        if out_long_filtered:
            print(f"[filter] wrote long: {out_long_filtered}", file=sys.stderr)

    return kept

# ----------------------------- Stats computation (exact medians) -----------------------------

def compute_per_sample_mod_stats_from_dedup(
    dedup_tsv: str,
    out_prefix: str,
    tag: str,
    tmpdir: str,
    mem: str,
    parallel: int,
    verbose: bool = False
):
    """
    Reproduces your three stats outputs for this TAG:
      1) {base}_{tag}__per_sample_mod_site_stats.tsv
      2) {base}_{tag}__per_sample_mod_tx_stats.tsv
      3) {base}_{tag}__per_tx_mod_stats.tsv

    Uses external sorts + streaming reducers to avoid large RAM.
    """

    base = out_prefix
    out1 = f"{base}_{tag}__per_sample_mod_site_stats.tsv"
    out2 = f"{base}_{tag}__per_sample_mod_tx_stats.tsv"
    out3 = f"{base}_{tag}__per_tx_mod_stats.tsv"

    ensure_dir(os.path.dirname(out1) or ".")

    # --- Build site_sample.tsv: one row per (sample,mod,chrom,start0,end0,strand) with summed Nmod,cov ---
    # dedup columns: sample zn chrom start0 end0 strand mod gid gname cov nmod ...
    # We want group key: sample(1), mod(7), chrom(3), start0(4), end0(5), strand(6)
    site_sample_uns = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_site_sample_uns.tsv")
    site_sample_sorted = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_site_sample_sorted.tsv")
    site_sample_red = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_site_sample.tsv")

    with open(site_sample_uns, "w") as out:
        with open(dedup_tsv, "r") as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 18:
                    continue
                sample = p[0]
                chrom, start0, end0, strand, mod = p[2], p[3], p[4], p[5], p[6]
                cov, nmod = p[9], p[10]
                out.write("\t".join([sample, mod, chrom, start0, end0, strand, nmod, cov]) + "\n")

    # sort by group key: sample,mod,chrom,start0,end0,strand
    run_sort(site_sample_uns, site_sample_sorted,
             keys=[(1,1,False),(2,2,False),(3,3,False),(4,4,True),(5,5,True),(6,6,False)],
             tmpdir=tmpdir, mem=mem, parallel=parallel, unique=False)

    # reduce groups
    with open(site_sample_red, "w") as out:
        # columns: sample mod chrom start0 end0 strand site_Nmod site_cov detected_site
        curr = None
        sum_nmod = 0
        sum_cov = 0
        n_sites_total_by_sm = defaultdict(int)
        n_sites_detected_by_sm = defaultdict(int)
        total_nmod_by_sm = defaultdict(int)
        total_cov_by_sm = defaultdict(int)

        # We'll write site rows to file for median computation
        # We'll write two metric files: cov_metric.tsv and nmod_metric.tsv, both per (sample,mod,metric)
        cov_metric = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_site_cov_metric.tsv")
        nmod_metric = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_site_nmod_metric.tsv")
        cov_mh = open(cov_metric, "w")
        nmod_mh = open(nmod_metric, "w")

        def flush_site():
            nonlocal curr, sum_nmod, sum_cov
            if curr is None:
                return
            sample, mod, chrom, start0, end0, strand = curr
            detected = 1 if sum_nmod > 0 else 0
            out.write("\t".join([sample, mod, chrom, start0, end0, strand, str(sum_nmod), str(sum_cov), str(detected)]) + "\n")

            sm = (sample, mod)
            n_sites_total_by_sm[sm] += 1
            n_sites_detected_by_sm[sm] += detected
            total_nmod_by_sm[sm] += sum_nmod
            total_cov_by_sm[sm] += sum_cov

            cov_mh.write(f"{sample}\t{mod}\t{sum_cov}\n")
            nmod_mh.write(f"{sample}\t{mod}\t{sum_nmod}\n")

        with open(site_sample_sorted, "r") as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 8:
                    continue
                key = (p[0], p[1], p[2], p[3], p[4], p[5])
                nmod = safe_int(p[6]); cov = safe_int(p[7])
                if curr is None:
                    curr = key; sum_nmod = nmod; sum_cov = cov
                elif key == curr:
                    sum_nmod += nmod; sum_cov += cov
                else:
                    flush_site()
                    curr = key; sum_nmod = nmod; sum_cov = cov

        flush_site()
        cov_mh.close()
        nmod_mh.close()

    # compute exact medians per (sample,mod) for site_cov and site_nmod via two-pass select
    def compute_median_from_metric(metric_path: str, metric_is_numeric: bool, out_map: Dict[Tuple[str,str], float], label: str):
        # Pass A: counts per group
        counts_path_uns = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_{label}_counts_uns.tsv")
        counts_sorted = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_{label}_counts_sorted.tsv")
        counts_red = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_{label}_counts.tsv")

        # group key = sample,mod
        run_sort(metric_path, counts_sorted,
                 keys=[(1,1,False),(2,2,False)], tmpdir=tmpdir, mem=mem, parallel=parallel, unique=False)

        with open(counts_red, "w") as out:
            curr = None
            c = 0
            with open(counts_sorted, "r") as f:
                for ln in f:
                    p = ln.rstrip("\n").split("\t")
                    key = (p[0], p[1])
                    if curr is None:
                        curr = key; c = 1
                    elif key == curr:
                        c += 1
                    else:
                        out.write(f"{curr[0]}\t{curr[1]}\t{c}\n")
                        curr = key; c = 1
            if curr is not None:
                out.write(f"{curr[0]}\t{curr[1]}\t{c}\n")

        # load counts (small: sample*mod)
        counts = {}
        with open(counts_red, "r") as f:
            for ln in f:
                s, m, c = ln.rstrip("\n").split("\t")
                counts[(s,m)] = int(c)

        # Pass B: sort by sample,mod,metric and select middle(s)
        metric_sorted = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_{label}_metric_sorted.tsv")
        run_sort(metric_path, metric_sorted,
                 keys=[(1,1,False),(2,2,False),(3,3,True)],
                 tmpdir=tmpdir, mem=mem, parallel=parallel, unique=False)

        # streaming select
        # median positions (0-indexed): if n odd -> mid=n//2; if even -> (n//2 -1, n//2)
        idx_in_group = 0
        curr = None
        need1 = need2 = None
        v1 = v2 = None

        with open(metric_sorted, "r") as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                key = (p[0], p[1])
                val = float(p[2])

                if curr is None or key != curr:
                    # flush previous group
                    if curr is not None:
                        if need2 is None:
                            out_map[curr] = float(v1) if v1 is not None else 0.0
                        else:
                            out_map[curr] = ((float(v1) if v1 is not None else 0.0) + (float(v2) if v2 is not None else 0.0)) / 2.0
                    # init new
                    curr = key
                    n = counts.get(curr, 0)
                    if n <= 0:
                        need1 = need2 = None
                    else:
                        if n % 2 == 1:
                            need1 = n // 2
                            need2 = None
                        else:
                            need1 = (n // 2) - 1
                            need2 = (n // 2)
                    idx_in_group = 0
                    v1 = v2 = None

                # capture if needed
                if need1 is not None and idx_in_group == need1:
                    v1 = val
                if need2 is not None and idx_in_group == need2:
                    v2 = val

                idx_in_group += 1

            # flush last
            if curr is not None:
                if need2 is None:
                    out_map[curr] = float(v1) if v1 is not None else 0.0
                else:
                    out_map[curr] = ((float(v1) if v1 is not None else 0.0) + (float(v2) if v2 is not None else 0.0)) / 2.0

    med_site_cov = {}
    med_site_nmod = {}
    compute_median_from_metric(os.path.join(tmpdir, f"__{os.getpid()}_{tag}_site_cov_metric.tsv"), True, med_site_cov, "site_cov")
    compute_median_from_metric(os.path.join(tmpdir, f"__{os.getpid()}_{tag}_site_nmod_metric.tsv"), True, med_site_nmod, "site_nmod")

    # compute means easily from totals / counts
    # write out1
    with open(out1, "w") as f:
        hdr = [
            "sample","mod_code","n_sites_total","n_sites_detected","total_Nmod","total_cov","overall_stoich",
            "mean_site_cov","median_site_cov","mean_site_Nmod","median_site_Nmod"
        ]
        f.write("\t".join(hdr) + "\n")

        rows = []
        for (sample, mod), n_sites in n_sites_total_by_sm.items():
            tn = total_nmod_by_sm[(sample, mod)]
            tc = total_cov_by_sm[(sample, mod)]
            nd = n_sites_detected_by_sm[(sample, mod)]
            overall = (tn / tc) if tc > 0 else 0.0
            mean_cov = (tc / n_sites) if n_sites > 0 else 0.0
            mean_nm = (tn / n_sites) if n_sites > 0 else 0.0
            med_cov = med_site_cov.get((sample, mod), 0.0)
            med_nm = med_site_nmod.get((sample, mod), 0.0)
            rows.append((mod, sample, n_sites, nd, tn, tc, overall, mean_cov, med_cov, mean_nm, med_nm))

        rows.sort(key=lambda x: (x[0], x[1]))
        for mod, sample, nst, nsd, tn, tc, ov, mc, medc, mnm, mednm in rows:
            f.write(
                f"{sample}\t{mod}\t{nst}\t{nsd}\t{tn}\t{tc}\t{ov:.6f}\t{mc:.6f}\t{medc:.6f}\t{mnm:.6f}\t{mednm:.6f}\n"
            )

    # --- Build per-tx stats (out3) ---
    # Create site_tx rows: group by (sample,mod,zn,chrom,start0,end0,strand) summing Nmod,cov
    site_tx_uns = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_site_tx_uns.tsv")
    site_tx_sorted = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_site_tx_sorted.tsv")

    with open(site_tx_uns, "w") as out:
        with open(dedup_tsv, "r") as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 18:
                    continue
                sample = p[0]; zn = p[1]
                chrom, start0, end0, strand, mod = p[2], p[3], p[4], p[5], p[6]
                cov, nmod = p[9], p[10]
                out.write("\t".join([sample, mod, zn, chrom, start0, end0, strand, nmod, cov]) + "\n")

    run_sort(site_tx_uns, site_tx_sorted,
             keys=[(1,1,False),(2,2,False),(3,3,True),(4,4,False),(5,5,True),(6,6,True),(7,7,False)],
             tmpdir=tmpdir, mem=mem, parallel=parallel, unique=False)

    # reduce site_tx -> per_tx aggregates
    # We'll also create metric files for medians over per-tx metrics later.
    per_tx_rows = []  # small enough to stream-write; we will not store; just write directly to out3 and metric files.
    tx_det_metric = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_tx_det_metric.tsv")
    tx_nmod_metric = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_tx_nmod_metric.tsv")
    tx_sto_metric = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_tx_sto_metric.tsv")

    tx_det_mh = open(tx_det_metric, "w")
    tx_nmod_mh = open(tx_nmod_metric, "w")
    tx_sto_mh = open(tx_sto_metric, "w")

    # per sample/mod aggregates over tx:
    tx_set = defaultdict(set)
    sum_det_sites = defaultdict(int)
    sum_total_nmod_per_tx = defaultdict(int)
    sum_tx_sto = defaultdict(float)

    # for means we also need counts of tx per sample/mod
    # but note: out3 includes all tx groups; we compute n_tx from tx_set size.

    with open(out3, "w") as out:
        out.write("\t".join(["sample","mod_code","ZN_transcript_index","n_sites_total","n_sites_detected","total_Nmod","total_cov","tx_stoich"]) + "\n")

        # Stage 1 reduce: site_tx_sorted -> per_tx accumulators
        # We accumulate within a tx group: (sample,mod,zn) and within that count sites etc.
        curr_site = None
        site_sum_nmod = 0
        site_sum_cov = 0

        # per-tx running:
        curr_tx = None
        tx_n_sites_total = 0
        tx_n_sites_detected = 0
        tx_total_nmod = 0
        tx_total_cov = 0

        def flush_site_into_tx():
            nonlocal tx_n_sites_total, tx_n_sites_detected, tx_total_nmod, tx_total_cov
            nonlocal site_sum_nmod, site_sum_cov
            if curr_site is None:
                return
            tx_n_sites_total += 1
            if site_sum_nmod > 0:
                tx_n_sites_detected += 1
            tx_total_nmod += site_sum_nmod
            tx_total_cov += site_sum_cov

        def flush_tx():
            nonlocal curr_tx, tx_n_sites_total, tx_n_sites_detected, tx_total_nmod, tx_total_cov
            if curr_tx is None:
                return
            sample, mod, zn = curr_tx
            sto = (tx_total_nmod / tx_total_cov) if tx_total_cov > 0 else 0.0
            out.write(f"{sample}\t{mod}\t{zn}\t{tx_n_sites_total}\t{tx_n_sites_detected}\t{tx_total_nmod}\t{tx_total_cov}\t{sto:.6f}\n")

            sm = (sample, mod)
            tx_set[sm].add(int(zn))
            sum_det_sites[sm] += tx_n_sites_detected
            sum_total_nmod_per_tx[sm] += tx_total_nmod
            sum_tx_sto[sm] += sto

            # metrics for medians (per sample/mod distributions)
            tx_det_mh.write(f"{sample}\t{mod}\t{tx_n_sites_detected}\n")
            tx_nmod_mh.write(f"{sample}\t{mod}\t{tx_total_nmod}\n")
            tx_sto_mh.write(f"{sample}\t{mod}\t{sto:.6f}\n")

            # reset
            tx_n_sites_total = 0
            tx_n_sites_detected = 0
            tx_total_nmod = 0
            tx_total_cov = 0

        with open(site_tx_sorted, "r") as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                # columns: sample,mod,zn,chrom,start0,end0,strand,nmod,cov
                sample, mod, zn = p[0], p[1], p[2]
                chrom, start0, end0, strand = p[3], p[4], p[5], p[6]
                nmod = safe_int(p[7]); cov = safe_int(p[8])

                site_key = (sample, mod, zn, chrom, start0, end0, strand)
                tx_key = (sample, mod, zn)

                if curr_site is None:
                    curr_site = site_key
                    curr_tx = tx_key
                    site_sum_nmod = nmod
                    site_sum_cov = cov
                    continue

                if site_key == curr_site:
                    site_sum_nmod += nmod
                    site_sum_cov += cov
                    continue

                # new site
                # if tx changed, flush site then tx
                if tx_key != curr_tx:
                    flush_site_into_tx()
                    flush_tx()
                    # reset for new tx
                    curr_tx = tx_key
                    # start new site
                    curr_site = site_key
                    site_sum_nmod = nmod
                    site_sum_cov = cov
                else:
                    # same tx, flush site and continue
                    flush_site_into_tx()
                    curr_site = site_key
                    site_sum_nmod = nmod
                    site_sum_cov = cov

        # flush last site + tx
        if curr_site is not None:
            flush_site_into_tx()
        flush_tx()

    tx_det_mh.close()
    tx_nmod_mh.close()
    tx_sto_mh.close()

    # Build per_sample_tx (out2) with exact medians for:
    # detected_sites_per_tx, total_Nmod_per_tx, tx_stoich
    med_det = {}
    med_nmod = {}
    med_sto = {}

    def compute_median_metric(metric_path: str, out_map: Dict[Tuple[str,str], float], label: str):
        # counts
        metric_sorted = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_{label}_sorted.tsv")
        run_sort(metric_path, metric_sorted,
                 keys=[(1,1,False),(2,2,False),(3,3,True)],
                 tmpdir=tmpdir, mem=mem, parallel=parallel, unique=False)

        # compute counts per group (sample,mod) quickly by streaming
        counts = defaultdict(int)
        with open(metric_sorted, "r") as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                counts[(p[0], p[1])] += 1

        # select medians by second pass on the already sorted file
        curr = None
        idx = 0
        need1 = need2 = None
        v1 = v2 = None

        with open(metric_sorted, "r") as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                key = (p[0], p[1])
                val = float(p[2])

                if curr is None or key != curr:
                    if curr is not None:
                        if need2 is None:
                            out_map[curr] = float(v1) if v1 is not None else 0.0
                        else:
                            out_map[curr] = ((float(v1) if v1 is not None else 0.0) + (float(v2) if v2 is not None else 0.0)) / 2.0

                    curr = key
                    n = counts[curr]
                    if n % 2 == 1:
                        need1 = n // 2; need2 = None
                    else:
                        need1 = (n // 2) - 1; need2 = (n // 2)
                    idx = 0
                    v1 = v2 = None

                if idx == need1:
                    v1 = val
                if need2 is not None and idx == need2:
                    v2 = val
                idx += 1

            if curr is not None:
                if need2 is None:
                    out_map[curr] = float(v1) if v1 is not None else 0.0
                else:
                    out_map[curr] = ((float(v1) if v1 is not None else 0.0) + (float(v2) if v2 is not None else 0.0)) / 2.0

    compute_median_metric(tx_det_metric, med_det, "tx_det")
    compute_median_metric(tx_nmod_metric, med_nmod, "tx_nmod")
    compute_median_metric(tx_sto_metric, med_sto, "tx_sto")

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
            if n_tx <= 0:
                continue
            mean_det = sum_det_sites[(sample, mod)] / n_tx
            mean_nmod = sum_total_nmod_per_tx[(sample, mod)] / n_tx
            mean_sto = sum_tx_sto[(sample, mod)] / n_tx
            rows.append((
                mod, sample, n_tx,
                mean_det, med_det.get((sample, mod), 0.0),
                mean_nmod, med_nmod.get((sample, mod), 0.0),
                mean_sto, med_sto.get((sample, mod), 0.0),
            ))

        rows.sort(key=lambda x: (x[0], x[1]))
        for mod, sample, n_tx, md, med_d, mn, med_n, ms, med_s in rows:
            f.write(
                f"{sample}\t{mod}\t{n_tx}\t"
                f"{md:.6f}\t{med_d:.6f}\t{mn:.6f}\t{med_n:.6f}\t{ms:.6f}\t{med_s:.6f}\n"
            )

    if verbose:
        print(f"[stats {tag}] wrote {out1}", file=sys.stderr)
        print(f"[stats {tag}] wrote {out2}", file=sys.stderr)
        print(f"[stats {tag}] wrote {out3}", file=sys.stderr)

# ----------------------------- Per-gene outputs + pivots -----------------------------

def generate_per_gene_outputs_from_dedup(
    dedup_tsv: str,
    out_prefix: str,
    tag: str,
    write_per_gene: bool,
    write_pivots: bool,
    tmpdir: str,
    mem: str,
    parallel: int,
    verbose: bool = False
):
    """
    Writes per-gene×mod tables and pivots under:
      <out_prefix>_<TAG>__per_gene_mod/

    This is done via:
      - extracting per-gene rows from dedup TSV
      - sorting by (gene_name, mod_code, chrom, start0, end0, strand, ZN, sample)
      - streaming per gene/mod group:
          * optionally write row TSV
          * optionally accumulate pivot dicts for that group and write pivot TSVs
    """
    if not (write_per_gene or write_pivots):
        return

    out_dir = f"{out_prefix}_{tag}__per_gene_mod"
    ensure_dir(out_dir)
    prefix_base = os.path.basename(out_prefix)

    # Build per_gene_uns.tsv: columns in PER_GENE_COLS order (no header)
    per_gene_uns = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_per_gene_uns.tsv")
    per_gene_sorted = os.path.join(tmpdir, f"__{os.getpid()}_{tag}_per_gene_sorted.tsv")

    with open(per_gene_uns, "w") as out:
        with open(dedup_tsv, "r") as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                if len(p) < 18:
                    continue
                sample, zn = p[0], p[1]
                chrom, start0, end0, strand, mod = p[2], p[3], p[4], p[5], p[6]
                gid, gname = p[7], p[8]
                cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall, frac = p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17]
                out.write("\t".join([
                    gname, gid, mod, chrom, start0, end0, strand, zn, sample,
                    cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall, frac
                ]) + "\n")

    # Sort by gene_name, mod_code, chrom, start0, end0, strand, ZN, sample
    run_sort(per_gene_uns, per_gene_sorted,
             keys=[(1,1,False),(3,3,False),(4,4,False),(5,5,True),(6,6,True),(7,7,False),(8,8,True),(9,9,False)],
             tmpdir=tmpdir, mem=mem, parallel=parallel, unique=False)

    def basepath(gene_name: str, mod: str) -> str:
        safe_g = sanitize_filename_token(gene_name if gene_name else "NA")
        safe_mod = sanitize_filename_token(str(mod))
        return os.path.join(out_dir, f"{prefix_base}__{safe_g}__{safe_mod}")

    # Stream groups
    curr_gm = None  # (gene_name, gene_id, mod)
    row_fh = None

    # Pivot maps for current gene/mod group:
    # idx = (chrom,start0,end0,strand,ZN) -> {sample: value}
    piv_cov = defaultdict(dict)
    piv_nmod = defaultdict(dict)
    piv_frac = defaultdict(dict)
    samples_seen = set()

    def flush_group():
        nonlocal curr_gm, row_fh, piv_cov, piv_nmod, piv_frac, samples_seen
        if curr_gm is None:
            return
        gene_name, gene_id, mod = curr_gm
        bp = basepath(gene_name, mod)

        if row_fh is not None:
            row_fh.close()
            row_fh = None

        if write_pivots:
            samples = sorted(samples_seen)
            idxs = sorted(piv_cov.keys(), key=lambda t: (t[0], int(t[1]), int(t[4])))

            def write_pivot(path, m, is_float: bool):
                with open(path, "w") as f:
                    f.write("\t".join(["chrom","start0","end0","strand","ZN_transcript_index"] + samples) + "\n")
                    for idx in idxs:
                        chrom, s0, e0, strand, zn = idx
                        row = [chrom, str(s0), str(e0), strand, str(zn)]
                        smap = m.get(idx, {})
                        for s in samples:
                            v = smap.get(s, 0)
                            if is_float:
                                row.append(f"{float(v):.6f}")
                            else:
                                row.append(str(int(v)))
                        f.write("\t".join(row) + "\n")

            write_pivot(f"{bp}_cov_pivot.tsv", piv_cov, is_float=False)
            write_pivot(f"{bp}_Nmod_pivot.tsv", piv_nmod, is_float=False)
            write_pivot(f"{bp}_frac_pivot.tsv", piv_frac, is_float=True)

        piv_cov = defaultdict(dict)
        piv_nmod = defaultdict(dict)
        piv_frac = defaultdict(dict)
        samples_seen = set()

    with open(per_gene_sorted, "r") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 18:
                continue
            gname, gid, mod = p[0], p[1], p[2]
            chrom, start0, end0, strand, zn, sample = p[3], p[4], p[5], p[6], p[7], p[8]
            cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall, frac = p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17]

            gm = (gname, gid, mod)
            if curr_gm is None:
                curr_gm = gm
                if write_per_gene:
                    bp = basepath(gname, mod)
                    row_fh = open(f"{bp}.tsv", "w")
                    row_fh.write("\t".join(PER_GENE_COLS) + "\n")

            if gm != curr_gm:
                flush_group()
                curr_gm = gm
                if write_per_gene:
                    bp = basepath(gname, mod)
                    row_fh = open(f"{bp}.tsv", "w")
                    row_fh.write("\t".join(PER_GENE_COLS) + "\n")

            if write_per_gene and row_fh is not None:
                row_fh.write("\t".join([gname, gid, mod, chrom, start0, end0, strand, zn, sample,
                                        cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall, frac]) + "\n")

            if write_pivots:
                samples_seen.add(sample)
                idx = (chrom, int(start0), int(end0), strand, int(zn))
                # "first" semantics
                if sample not in piv_cov[idx]:
                    piv_cov[idx][sample] = int(cov)
                    piv_nmod[idx][sample] = int(nmod)
                    piv_frac[idx][sample] = float(frac)

    flush_group()

    if verbose:
        print(f"[per-gene {tag}] wrote outputs under {out_dir}", file=sys.stderr)

# ----------------------------- Main -----------------------------

def main():
    args = parse_args()

    beds = iter_numbered_beds(args.modkit_dir)
    if not beds:
        sys.exit(f"No numbered ZN partition files found under {args.modkit_dir}")

    gene_index = load_gene_intervals_from_gtf(args.gtf, verbose=args.verbose)

    # Create a private working dir inside tmpdir for this run
    workdir = tempfile.mkdtemp(prefix=f"aggregate_by_gene_{os.getpid()}_", dir=args.tmpdir)
    if args.verbose:
        print(f"[tmp] workdir={workdir}", file=sys.stderr)

    try:
        # Stage 1: normalize
        norm_tsv = os.path.join(workdir, "norm.tsv")
        normalize_to_tsv(beds, gene_index, norm_tsv, verbose=args.verbose)

        # Stage 2: sort for dedup
        norm_sorted = os.path.join(workdir, "norm.sorted.tsv")

        # NORM columns (no header):
        # 1 sample,2 zn,3 chrom,4 start0,5 end0,6 strand,7 mod,8 gid,9 gname, ...
        # Dedup key (same as pandas version):
        # sample, mod_code, ZN, chrom, start0, end0, strand, gene_id, gene_name
        run_sort(
            norm_tsv, norm_sorted,
            keys=[
                (1,1,False),  # sample
                (7,7,False),  # mod_code
                (2,2,True),   # zn
                (3,3,False),  # chrom
                (4,4,True),   # start0
                (5,5,True),   # end0
                (6,6,False),  # strand
                (8,8,False),  # gene_id
                (9,9,False),  # gene_name
            ],
            tmpdir=workdir, mem=args.sort_mem, parallel=args.sort_parallel, unique=False
        )

        base = args.out_prefix

        # Stage 3: dedup reduce (RAW) and collect passing-sites
        dedup_raw = os.path.join(workdir, "dedup.RAW.tsv")
        passing_sites = os.path.join(workdir, "passing_sites.unsorted.tsv") if args.filter_enable else None

        out_long_raw = f"{base}_RAW_sites_long.tsv" if (args.emit_raw and args.write_long) else None

        if args.emit_raw or args.filter_enable or args.emit_filt:
            dedup_reduce_sorted(
                sorted_norm_tsv=norm_sorted,
                out_dedup_tsv=dedup_raw,
                out_long_tsv=out_long_raw,
                out_passing_sites_tsv=passing_sites,
                min_cov=args.min_cov,
                filter_enable=args.filter_enable,
                count_diff_factor=args.count_diff_factor,
                mod_fail_margin=args.mod_fail_margin,
                verbose=args.verbose
            )

        # Per-sample stats + per-gene for RAW
        if args.emit_raw:
            compute_per_sample_mod_stats_from_dedup(
                dedup_tsv=dedup_raw,
                out_prefix=base,
                tag="RAW",
                tmpdir=workdir,
                mem=args.sort_mem,
                parallel=args.sort_parallel,
                verbose=args.verbose
            )
            generate_per_gene_outputs_from_dedup(
                dedup_tsv=dedup_raw,
                out_prefix=base,
                tag="RAW",
                write_per_gene=args.write_raw_per_gene,
                write_pivots=args.write_pivots,
                tmpdir=workdir,
                mem=args.sort_mem,
                parallel=args.sort_parallel,
                verbose=args.verbose
            )

        # Stage 4: FILTERED subset via merge-join (if enabled)
        if args.emit_filt:
            if args.filter_enable:
                # sort+unique passing-sites by site key
                passing_sorted = os.path.join(workdir, "passing_sites.sorted.unique.tsv")
                # passing-sites columns: chrom start0 end0 strand mod_code
                run_sort(
                    passing_sites, passing_sorted,
                    keys=[(1,1,False),(2,2,True),(3,3,True),(4,4,False),(5,5,False)],
                    tmpdir=workdir, mem=args.sort_mem, parallel=args.sort_parallel, unique=True
                )

                # sort dedup_raw by site key for streaming join:
                # dedup columns: sample(1) zn(2) chrom(3) start0(4) end0(5) strand(6) mod(7) ...
                dedup_by_site = os.path.join(workdir, "dedup.by_site.sorted.tsv")
                run_sort(
                    dedup_raw, dedup_by_site,
                    keys=[(3,3,False),(4,4,True),(5,5,True),(6,6,False),(7,7,False),
                          (1,1,False),(2,2,True),(8,8,False),(9,9,False)],  # tie-breakers stable-ish
                    tmpdir=workdir, mem=args.sort_mem, parallel=args.sort_parallel, unique=False
                )

                dedup_filt = os.path.join(workdir, "dedup.FILTERED.tsv")
                out_long_filt = f"{base}_FILTERED_sites_long.tsv" if args.write_long else None

                filter_dedup_by_passing_sites(
                    dedup_by_site_sorted=dedup_by_site,
                    passing_sites_sorted=passing_sorted,
                    out_dedup_filtered=dedup_filt,
                    out_long_filtered=out_long_filt,
                    verbose=args.verbose
                )

            else:
                # no filtering requested; FILTERED == RAW semantics in your original code
                dedup_filt = dedup_raw
                out_long_filt = f"{base}_FILTERED_sites_long.tsv" if args.write_long else None
                if out_long_filt is not None:
                    # just copy RAW long if it exists, else write from dedup_raw
                    if out_long_raw and os.path.exists(out_long_raw):
                        ensure_dir(os.path.dirname(out_long_filt) or ".")
                        shutil.copyfile(out_long_raw, out_long_filt)
                    else:
                        ensure_dir(os.path.dirname(out_long_filt) or ".")
                        with open(out_long_filt, "w") as out:
                            out.write("\t".join(LONG_HEADER) + "\n")
                            with open(dedup_raw, "r") as f:
                                for ln in f:
                                    p = ln.rstrip("\n").split("\t")
                                    if len(p) < 18:
                                        continue
                                    sample, zn, chrom, start0, end0, strand, mod, gid, gname = p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]
                                    cov, nmod = p[9], p[10]
                                    ncan, nother, ndel, nfail, ndiff, nnocall, frac = p[11], p[12], p[13], p[14], p[15], p[16], p[17]
                                    out.write(
                                        f"{sample}\t{zn}\t{chrom}\t{start0}\t{end0}\t{strand}\t{mod}\t"
                                        f"{cov}\t{nmod}\t{frac}\t{gid}\t{gname}\t"
                                        f"{ncan}\t{nother}\t{ndel}\t{nfail}\t{ndiff}\t{nnocall}\n"
                                    )

            # Stats + per-gene for FILTERED
            compute_per_sample_mod_stats_from_dedup(
                dedup_tsv=dedup_filt,
                out_prefix=base,
                tag="FILTERED",
                tmpdir=workdir,
                mem=args.sort_mem,
                parallel=args.sort_parallel,
                verbose=args.verbose
            )
            generate_per_gene_outputs_from_dedup(
                dedup_tsv=dedup_filt,
                out_prefix=base,
                tag="FILTERED",
                write_per_gene=args.write_filtered_per_gene,
                write_pivots=args.write_pivots,
                tmpdir=workdir,
                mem=args.sort_mem,
                parallel=args.sort_parallel,
                verbose=args.verbose
            )

        print("[OK] ZN aggregation complete.")

    finally:
        if args.keep_intermediates:
            print(f"[tmp] kept intermediates at {workdir}", file=sys.stderr)
        else:
            try:
                shutil.rmtree(workdir)
            except Exception:
                pass

if __name__ == "__main__":
    main()


