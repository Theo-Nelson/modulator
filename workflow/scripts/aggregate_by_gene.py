#!/usr/bin/env python3
"""
aggregate_by_gene.py (pure-Python external sort; OOM-safe; no GNU sort)

What it does (ZN aggregation):
- Reads numbered ZN partition bedMethyl files under --modkit-dir (per-sample subdirs)
  e.g. results/modkit_zn/<sample>/**/<N>.bed(.gz)
- Assigns each site to a gene (coarse gene spans derived from the assembler GTF).
- Normalizes all rows into a compact TSV (no header).
- Performs disk-backed external sorts in pure Python (chunk sort + k-way merge).
- Deduplicates by:
    (sample, mod_code, ZN_transcript_index, chrom, start0, end0, strand, gene_id, gene_name)
  summing numeric columns.
- Optional site filtering:
    - FAIL if Ndiff > count_diff_factor * Nvalid_cov
    - FAIL if Nmod <= Nfail + mod_fail_margin
  FILTERED outputs are subset of RAW by site_key (chrom,start0,end0,strand,mod_code).

Outputs (same naming convention as your skeleton):
- <out_prefix>_RAW_sites_long.tsv                         (optional)
- <out_prefix>_FILTERED_sites_long.tsv                    (optional; required by your Snakemake)
- <out_prefix>_RAW__per_sample_mod_site_stats.tsv         (if emit_raw)
- <out_prefix>_RAW__per_sample_mod_tx_stats.tsv           (if emit_raw)
- <out_prefix>_RAW__per_tx_mod_stats.tsv                  (if emit_raw)
- <out_prefix>_FILTERED__per_sample_mod_site_stats.tsv    (if emit_filtered)
- <out_prefix>_FILTERED__per_sample_mod_tx_stats.tsv      (if emit_filtered)
- <out_prefix>_FILTERED__per_tx_mod_stats.tsv             (if emit_filtered)
- <out_prefix>_RAW__per_gene_mod/....                     (optional)
- <out_prefix>_FILTERED__per_gene_mod/....                (optional)

No dependency on GNU sort. Uses:
- chunked in-memory sorting (bounded by --chunk-lines)
- k-way merge across sorted chunks (heapq)
"""

import os
import sys
import re
import gzip
import argparse
import tempfile
import shutil
import heapq
from collections import defaultdict, namedtuple
from typing import Dict, Tuple, List, Optional

# ----------------------------- Constants -----------------------------

BED_COLS = [
    "chrom", "start0", "end0", "mod_code", "score", "strand",
    "start0_compat", "end0_compat", "rgb",
    "Nvalid_cov", "frac_modified",
    "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall",
]

# Normalized (pre-dedup) TSV columns we emit (NO header):
# sample, zn, chrom, start0, end0, strand, mod_code, gene_id, gene_name, plus numeric columns
NORM_COLS = [
    "sample", "ZN_transcript_index", "chrom", "start0", "end0", "strand", "mod_code",
    "gene_id", "gene_name",
    "Nvalid_cov", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall",
]

LONG_HEADER = [
    "sample", "ZN_transcript_index", "chrom", "start0", "end0", "strand", "mod_code",
    "Nvalid_cov", "Nmod", "frac_modified", "gene_id", "gene_name",
    "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall",
]

PER_GENE_COLS = [
    "gene_name", "gene_id", "mod_code", "chrom", "start0", "end0", "strand",
    "ZN_transcript_index", "sample", "Nvalid_cov", "Nmod", "Ncanonical",
    "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall", "frac_modified"
]

# ----------------------------- CLI -----------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Aggregate ZN-partitioned modkit outputs per gene/mod with site-level filtering (pure Python external sort)"
    )
    ap.add_argument("--modkit-dir", required=True, help="Parent dir with per-sample subdirs containing numbered ZN .bed files")
    ap.add_argument("--gtf", required=True, help="Assembler GTF (with gene coordinates). exon/transcript/gene features used.")
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

    ap.add_argument("--write-raw-per-gene", dest="write_raw_per_gene", action="store_true")
    ap.add_argument("--no-write-raw-per-gene", dest="write_raw_per_gene", action="store_false"); ap.set_defaults(write_raw_per_gene=False)
    ap.add_argument("--write-filtered-per-gene", dest="write_filtered_per_gene", action="store_true")
    ap.add_argument("--no-write-filtered-per-gene", dest="write_filtered_per_gene", action="store_false"); ap.set_defaults(write_filtered_per_gene=True)

    ap.add_argument("--verbose", action="store_true")

    # Pure-Python external sort tuning
    ap.add_argument("--chunk-lines", type=int, default=2_000_000,
                    help="Lines per in-memory chunk during external sort (default: 2,000,000)")
    ap.add_argument("--tmpdir", type=str, default=os.environ.get("TMPDIR", "/tmp"),
                    help="Temp directory for intermediates (default: $TMPDIR or /tmp)")
    ap.add_argument("--keep-intermediates", action="store_true", help="Do not delete intermediates (debugging)")

    return ap.parse_args()

# ----------------------------- Utils -----------------------------

def ensure_dir(path: str):
    if path:
        os.makedirs(path, exist_ok=True)

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

# ----------------------------- Pure-Python external sort -----------------------------

def external_sort_tsv(
    in_path: str,
    out_path: str,
    key_func,
    tmpdir: str,
    chunk_lines: int,
    verbose: bool = False
):
    """
    External sort for a tab-delimited file (no header).
    Uses chunked in-memory sort + k-way merge across sorted chunks.
    """
    ensure_dir(os.path.dirname(out_path) or ".")
    ensure_dir(tmpdir)

    chunks: List[str] = []
    buf: List[str] = []
    n_in = 0

    def write_chunk(lines: List[str], idx: int) -> str:
        lines.sort(key=key_func)
        cpath = os.path.join(tmpdir, f"chunk_{idx:05d}.tsv")
        with open(cpath, "w") as f:
            f.writelines(lines)
        return cpath

    with open(in_path, "r") as fin:
        for ln in fin:
            if not ln:
                continue
            buf.append(ln)
            n_in += 1
            if len(buf) >= chunk_lines:
                chunks.append(write_chunk(buf, len(chunks)))
                buf = []
    if buf:
        chunks.append(write_chunk(buf, len(chunks)))

    if verbose:
        print(f"[pysort] {os.path.basename(in_path)}: {n_in} lines -> {len(chunks)} chunks", file=sys.stderr)

    # Single chunk: copy
    if len(chunks) == 1:
        shutil.copyfile(chunks[0], out_path)
        os.remove(chunks[0])
        return

    # Merge chunks
    fps = [open(p, "r") for p in chunks]
    try:
        heap = []
        for i, f in enumerate(fps):
            ln = f.readline()
            if ln:
                heap.append((key_func(ln), ln, i))
        heapq.heapify(heap)

        with open(out_path, "w") as out:
            while heap:
                _, ln, i = heapq.heappop(heap)
                out.write(ln)
                nxt = fps[i].readline()
                if nxt:
                    heapq.heappush(heap, (key_func(nxt), nxt, i))
    finally:
        for f in fps:
            try:
                f.close()
            except Exception:
                pass
        for p in chunks:
            try:
                os.remove(p)
            except Exception:
                pass

def uniq_sorted_file(in_sorted: str, out_path: str, key_func):
    """
    Input must already be sorted by key_func.
    Writes unique lines by key.
    """
    ensure_dir(os.path.dirname(out_path) or ".")
    prev = None
    with open(in_sorted, "r") as fin, open(out_path, "w") as out:
        for ln in fin:
            k = key_func(ln)
            if prev is None or k != prev:
                out.write(ln)
                prev = k

# ----------------------------- GTF interval indexing -----------------------------

Interval = namedtuple("Interval", ["start", "end", "gene_id", "gene_name", "strand"])

def load_gene_intervals_from_gtf(gtf_path: str, verbose=False) -> Dict[Tuple[str, str], List[Interval]]:
    """
    Coarse per-gene spans for site->gene mapping.
    Uses exon/transcript/gene features; takes union span per (chrom,strand,gene_id).
    """
    gene_bounds: Dict[Tuple[str, str, str], Tuple[int, int]] = {}
    gene_name_map: Dict[str, str] = {}

    with open_text(gtf_path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = parts
            if feature not in ("exon", "transcript", "gene"):
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

    by_cs: Dict[Tuple[str, str], List[Interval]] = defaultdict(list)
    for (chrom, strand, gid), (s, e) in gene_bounds.items():
        gname = gene_name_map.get(gid, gid)
        by_cs[(chrom, strand)].append(Interval(s, e, gid, gname, strand))

    for k in by_cs:
        by_cs[k].sort(key=lambda iv: (iv.start, iv.end))

    if verbose:
        total = sum(len(v) for v in by_cs.values())
        print(f"[info] loaded {total} gene spans from {gtf_path}", file=sys.stderr)

    return by_cs

def assign_gene(chrom: str, pos_start: int, pos_end: int, strand: str,
                gene_index: Dict[Tuple[str, str], List[Interval]]) -> Tuple[str, str]:
    """
    Return (gene_id, gene_name) by overlap on same strand; choose max-overlap.
    If none, try opposite strand first overlap.
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
    for files like '<sample>/<...>/<N>.bed' or '<N>.bed.gz'.
    Skips:
      - ungrouped.bed(.gz)
      - *_filtered_mod.bed(.gz)
    """
    out = []
    for root, _, files in os.walk(modkit_dir):
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
    for k in ["Nvalid_cov", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall"]:
        d[k] = safe_int(d[k])
    d["frac_modified"] = safe_float(d.get("frac_modified", 0.0), 0.0)
    return d

# ----------------------------- Stage 1: normalize -----------------------------

def normalize_to_tsv(beds: List[Tuple[str, str, str, int]], gene_index, out_tsv: str, verbose: bool = False) -> int:
    """
    Stream all numbered bed files and write normalized TSV (no header).
    Returns number of rows written.
    """
    ensure_dir(os.path.dirname(out_tsv) or ".")
    n = 0
    with open(out_tsv, "w") as out:
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
        print(f"[norm] wrote {n} rows -> {out_tsv}", file=sys.stderr)
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
    Reduce sorted normalized TSV to deduplicated rows.
    Writes:
      - out_dedup_tsv (no header):
        sample zn chrom start0 end0 strand mod gid gname cov nmod ncan nother ndel nfail ndiff nnocall frac
      - out_long_tsv (if not None) with LONG_HEADER
      - out_passing_sites_tsv (if filter enabled) with site_key lines: chrom start0 end0 strand mod
    """
    ensure_dir(os.path.dirname(out_dedup_tsv) or ".")
    if out_long_tsv:
        ensure_dir(os.path.dirname(out_long_tsv) or ".")
    if out_passing_sites_tsv and filter_enable:
        ensure_dir(os.path.dirname(out_passing_sites_tsv) or ".")

    # truncate outputs
    with open(out_dedup_tsv, "w"):
        pass
    dedup_fh = open(out_dedup_tsv, "a")

    long_fh = None
    if out_long_tsv is not None:
        long_fh = open(out_long_tsv, "w")
        long_fh.write("\t".join(LONG_HEADER) + "\n")

    pass_fh = None
    if filter_enable and out_passing_sites_tsv is not None:
        pass_fh = open(out_passing_sites_tsv, "w")

    dedup_written = 0

    def flush(curr_key, sums):
        nonlocal dedup_written
        if curr_key is None:
            return
        (sample, zn, chrom, start0, end0, strand, mod, gid, gname) = curr_key
        cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall = sums
        frac = frac_modified(nmod, cov, min_cov)

        dedup_fh.write(
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

    curr_key = None
    sums = [0] * 8

    with open(sorted_norm_tsv, "r") as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            p = ln.split("\t")
            # NORM:
            # 0 sample,1 zn,2 chrom,3 start0,4 end0,5 strand,6 mod,7 gid,8 gname
            # 9 cov,10 nmod,11 ncan,12 nother,13 ndel,14 nfail,15 ndiff,16 nnocall
            key = (p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8])
            vals = [
                safe_int(p[9]), safe_int(p[10]), safe_int(p[11]), safe_int(p[12]),
                safe_int(p[13]), safe_int(p[14]), safe_int(p[15]), safe_int(p[16]),
            ]

            if curr_key is None:
                curr_key = key
                sums = vals
            elif key == curr_key:
                for i in range(8):
                    sums[i] += vals[i]
            else:
                flush(curr_key, sums)
                curr_key = key
                sums = vals

    flush(curr_key, sums)

    dedup_fh.close()
    if long_fh:
        long_fh.close()
        if verbose:
            print(f"[dedup] wrote long: {out_long_tsv}", file=sys.stderr)
    if pass_fh:
        pass_fh.close()
        if verbose:
            print(f"[dedup] wrote passing sites: {out_passing_sites_tsv}", file=sys.stderr)

    if verbose:
        print(f"[dedup] wrote {dedup_written} rows -> {out_dedup_tsv}", file=sys.stderr)
    return dedup_written

# ----------------------------- Filtering join -----------------------------

def filter_dedup_by_passing_sites(
    dedup_by_site_sorted: str,
    passing_sites_sorted_unique: str,
    out_dedup_filtered: str,
    out_long_filtered: Optional[str],
    verbose: bool = False,
) -> int:
    """
    Streaming merge-join:
    - dedup_by_site_sorted sorted by site_key then tie-breakers
    - passing_sites_sorted_unique sorted by site_key unique

    Writes filtered dedup TSV (no header) + optional long TSV with header.
    Returns number of filtered dedup rows.
    """
    ensure_dir(os.path.dirname(out_dedup_filtered) or ".")
    with open(out_dedup_filtered, "w"):
        pass
    out_fh = open(out_dedup_filtered, "a")

    long_fh = None
    if out_long_filtered is not None:
        ensure_dir(os.path.dirname(out_long_filtered) or ".")
        long_fh = open(out_long_filtered, "w")
        long_fh.write("\t".join(LONG_HEADER) + "\n")

    def site_of_dedup(parts: List[str]) -> Tuple[str, int, int, str, str]:
        # dedup columns:
        # 0 sample,1 zn,2 chrom,3 start0,4 end0,5 strand,6 mod,...
        return (parts[2], int(parts[3]), int(parts[4]), parts[5], parts[6])

    ps = open(passing_sites_sorted_unique, "r")
    ps_line = ps.readline()
    ps_key = None
    if ps_line:
        p = ps_line.rstrip("\n").split("\t")
        ps_key = (p[0], int(p[1]), int(p[2]), p[3], p[4])

    kept = 0
    with open(dedup_by_site_sorted, "r") as d:
        for ln in d:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            parts = ln.split("\t")
            dk = site_of_dedup(parts)

            while ps_key is not None and ps_key < dk:
                ps_line = ps.readline()
                if not ps_line:
                    ps_key = None
                    break
                p = ps_line.rstrip("\n").split("\t")
                ps_key = (p[0], int(p[1]), int(p[2]), p[3], p[4])

            if ps_key is None:
                break

            if dk == ps_key:
                out_fh.write("\t".join(parts) + "\n")
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

    ps.close()
    out_fh.close()
    if long_fh:
        long_fh.close()

    if verbose:
        print(f"[filter] kept {kept} rows -> {out_dedup_filtered}", file=sys.stderr)
    return kept

# ----------------------------- Stats computation (exact medians via external sort) -----------------------------

def _median_from_sorted_metric(metric_sorted: str, verbose: bool = False) -> Dict[Tuple[str, str], float]:
    """
    metric_sorted must be sorted by (sample, mod, value_numeric) with no header.
    Each line: sample<TAB>mod<TAB>value
    Returns dict (sample,mod)->median (exact).
    """
    # First pass: counts per group
    counts = defaultdict(int)
    with open(metric_sorted, "r") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            counts[(p[0], p[1])] += 1

    # Second pass: select median positions per group
    out = {}
    curr = None
    idx = 0
    need1 = need2 = None
    v1 = v2 = None

    with open(metric_sorted, "r") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            key = (p[0], p[1])
            val = float(p[2])

            if curr is None or key != curr:
                if curr is not None:
                    if need2 is None:
                        out[curr] = float(v1) if v1 is not None else 0.0
                    else:
                        out[curr] = ((float(v1) if v1 is not None else 0.0) + (float(v2) if v2 is not None else 0.0)) / 2.0

                curr = key
                n = counts[curr]
                if n % 2 == 1:
                    need1 = n // 2
                    need2 = None
                else:
                    need1 = (n // 2) - 1
                    need2 = (n // 2)
                idx = 0
                v1 = v2 = None

            if idx == need1:
                v1 = val
            if need2 is not None and idx == need2:
                v2 = val
            idx += 1

        if curr is not None:
            if need2 is None:
                out[curr] = float(v1) if v1 is not None else 0.0
            else:
                out[curr] = ((float(v1) if v1 is not None else 0.0) + (float(v2) if v2 is not None else 0.0)) / 2.0

    return out

def compute_per_sample_mod_stats_from_dedup(
    dedup_tsv: str,
    out_prefix: str,
    tag: str,
    workdir: str,
    chunk_lines: int,
    verbose: bool = False
):
    """
    Writes:
      1) {out_prefix}_{tag}__per_sample_mod_site_stats.tsv
      2) {out_prefix}_{tag}__per_sample_mod_tx_stats.tsv
      3) {out_prefix}_{tag}__per_tx_mod_stats.tsv
    """
    base = out_prefix
    out1 = f"{base}_{tag}__per_sample_mod_site_stats.tsv"
    out2 = f"{base}_{tag}__per_sample_mod_tx_stats.tsv"
    out3 = f"{base}_{tag}__per_tx_mod_stats.tsv"
    ensure_dir(os.path.dirname(out1) or ".")

    # ---------- SITE STATS ----------
    # Build site_sample_uns: sample mod chrom start0 end0 strand nmod cov
    site_sample_uns = os.path.join(workdir, f"{tag}.site_sample.uns.tsv")
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

    def key_site_sample(line: str):
        p = line.rstrip("\n").split("\t")
        return (p[0], p[1], p[2], int(p[3]), int(p[4]), p[5])

    site_sample_sorted = os.path.join(workdir, f"{tag}.site_sample.sorted.tsv")
    external_sort_tsv(site_sample_uns, site_sample_sorted, key_site_sample, tmpdir=workdir, chunk_lines=chunk_lines, verbose=verbose)

    # Reduce per-site and create metric files for medians
    cov_metric = os.path.join(workdir, f"{tag}.site_cov.metric.tsv")
    nmod_metric = os.path.join(workdir, f"{tag}.site_nmod.metric.tsv")

    n_sites_total_by_sm = defaultdict(int)
    n_sites_detected_by_sm = defaultdict(int)
    total_nmod_by_sm = defaultdict(int)
    total_cov_by_sm = defaultdict(int)

    with open(cov_metric, "w") as cov_mh, open(nmod_metric, "w") as nmod_mh:
        curr = None
        sum_nmod = 0
        sum_cov = 0

        def flush_site():
            nonlocal curr, sum_nmod, sum_cov
            if curr is None:
                return
            sample, mod, *_ = curr
            detected = 1 if sum_nmod > 0 else 0
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
                key = (p[0], p[1], p[2], p[3], p[4], p[5])  # sample,mod,chrom,start0,end0,strand
                nmod = safe_int(p[6]); cov = safe_int(p[7])
                if curr is None:
                    curr = key
                    sum_nmod = nmod
                    sum_cov = cov
                elif key == curr:
                    sum_nmod += nmod
                    sum_cov += cov
                else:
                    flush_site()
                    curr = key
                    sum_nmod = nmod
                    sum_cov = cov
        flush_site()

    # Sort metrics by (sample,mod,value) and compute medians
    def key_metric(line: str):
        p = line.rstrip("\n").split("\t")
        return (p[0], p[1], float(p[2]))

    cov_metric_sorted = os.path.join(workdir, f"{tag}.site_cov.metric.sorted.tsv")
    nmod_metric_sorted = os.path.join(workdir, f"{tag}.site_nmod.metric.sorted.tsv")
    external_sort_tsv(cov_metric, cov_metric_sorted, key_metric, tmpdir=workdir, chunk_lines=chunk_lines, verbose=verbose)
    external_sort_tsv(nmod_metric, nmod_metric_sorted, key_metric, tmpdir=workdir, chunk_lines=chunk_lines, verbose=verbose)

    med_site_cov = _median_from_sorted_metric(cov_metric_sorted)
    med_site_nmod = _median_from_sorted_metric(nmod_metric_sorted)

    # Write out1
    with open(out1, "w") as f:
        hdr = [
            "sample", "mod_code", "n_sites_total", "n_sites_detected", "total_Nmod", "total_cov", "overall_stoich",
            "mean_site_cov", "median_site_cov", "mean_site_Nmod", "median_site_Nmod"
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

    # ---------- TX STATS ----------
    # Build site_tx_uns: sample mod zn chrom start0 end0 strand nmod cov
    site_tx_uns = os.path.join(workdir, f"{tag}.site_tx.uns.tsv")
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

    def key_site_tx(line: str):
        p = line.rstrip("\n").split("\t")
        return (p[0], p[1], int(p[2]), p[3], int(p[4]), int(p[5]), p[6])

    site_tx_sorted = os.path.join(workdir, f"{tag}.site_tx.sorted.tsv")
    external_sort_tsv(site_tx_uns, site_tx_sorted, key_site_tx, tmpdir=workdir, chunk_lines=chunk_lines, verbose=verbose)

    # Reduce site_tx_sorted into per-tx stats (out3) + metrics for medians in out2
    tx_det_metric = os.path.join(workdir, f"{tag}.tx_det.metric.tsv")
    tx_nmod_metric = os.path.join(workdir, f"{tag}.tx_nmod.metric.tsv")
    tx_sto_metric = os.path.join(workdir, f"{tag}.tx_sto.metric.tsv")

    # per sample/mod accumulators
    tx_set = defaultdict(set)
    sum_det_sites = defaultdict(int)
    sum_total_nmod_per_tx = defaultdict(int)
    sum_tx_sto = defaultdict(float)

    with open(tx_det_metric, "w") as det_mh, open(tx_nmod_metric, "w") as nm_mh, open(tx_sto_metric, "w") as sto_mh:
        with open(out3, "w") as out:
            out.write("\t".join(["sample", "mod_code", "ZN_transcript_index", "n_sites_total", "n_sites_detected", "total_Nmod", "total_cov", "tx_stoich"]) + "\n")

            curr_site = None
            site_sum_nmod = 0
            site_sum_cov = 0

            curr_tx = None
            tx_n_sites_total = 0
            tx_n_sites_detected = 0
            tx_total_nmod = 0
            tx_total_cov = 0

            def flush_site_into_tx():
                nonlocal tx_n_sites_total, tx_n_sites_detected, tx_total_nmod, tx_total_cov
                nonlocal site_sum_nmod, site_sum_cov, curr_site
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

                det_mh.write(f"{sample}\t{mod}\t{tx_n_sites_detected}\n")
                nm_mh.write(f"{sample}\t{mod}\t{tx_total_nmod}\n")
                sto_mh.write(f"{sample}\t{mod}\t{sto:.6f}\n")

                tx_n_sites_total = 0
                tx_n_sites_detected = 0
                tx_total_nmod = 0
                tx_total_cov = 0

            with open(site_tx_sorted, "r") as f:
                for ln in f:
                    p = ln.rstrip("\n").split("\t")
                    if len(p) < 9:
                        continue
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
                    if tx_key != curr_tx:
                        flush_site_into_tx()
                        flush_tx()
                        curr_tx = tx_key
                        curr_site = site_key
                        site_sum_nmod = nmod
                        site_sum_cov = cov
                    else:
                        flush_site_into_tx()
                        curr_site = site_key
                        site_sum_nmod = nmod
                        site_sum_cov = cov

            if curr_site is not None:
                flush_site_into_tx()
            flush_tx()

    # Sort tx metrics + compute medians for out2
    def key_metric3(line: str):
        p = line.rstrip("\n").split("\t")
        return (p[0], p[1], float(p[2]))

    det_sorted = os.path.join(workdir, f"{tag}.tx_det.metric.sorted.tsv")
    nm_sorted = os.path.join(workdir, f"{tag}.tx_nmod.metric.sorted.tsv")
    sto_sorted = os.path.join(workdir, f"{tag}.tx_sto.metric.sorted.tsv")

    external_sort_tsv(tx_det_metric, det_sorted, key_metric3, tmpdir=workdir, chunk_lines=chunk_lines, verbose=verbose)
    external_sort_tsv(tx_nmod_metric, nm_sorted, key_metric3, tmpdir=workdir, chunk_lines=chunk_lines, verbose=verbose)
    external_sort_tsv(tx_sto_metric, sto_sorted, key_metric3, tmpdir=workdir, chunk_lines=chunk_lines, verbose=verbose)

    med_det = _median_from_sorted_metric(det_sorted)
    med_nmod = _median_from_sorted_metric(nm_sorted)
    med_sto = _median_from_sorted_metric(sto_sorted)

    with open(out2, "w") as f:
        hdr = [
            "sample", "mod_code", "n_tx",
            "mean_detected_sites_per_tx", "median_detected_sites_per_tx",
            "mean_total_Nmod_per_tx", "median_total_Nmod_per_tx",
            "mean_tx_stoich", "median_tx_stoich"
        ]
        f.write("\t".join(hdr) + "\n")

        rows = []
        for (sample, mod), zset in tx_set.items():
            n_tx = len(zset)
            if n_tx <= 0:
                continue
            mean_det = sum_det_sites[(sample, mod)] / n_tx
            mean_nm = sum_total_nmod_per_tx[(sample, mod)] / n_tx
            mean_st = sum_tx_sto[(sample, mod)] / n_tx
            rows.append((
                mod, sample, n_tx,
                mean_det, med_det.get((sample, mod), 0.0),
                mean_nm, med_nmod.get((sample, mod), 0.0),
                mean_st, med_sto.get((sample, mod), 0.0),
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
    workdir: str,
    chunk_lines: int,
    verbose: bool = False
):
    """
    Writes per-gene×mod tables and pivots under:
      <out_prefix>_<TAG>__per_gene_mod/

    - row TSV per gene/mod: <prefix_base>__<gene>__<mod>.tsv
    - pivots: *_cov_pivot.tsv, *_Nmod_pivot.tsv, *_frac_pivot.tsv
    """
    if not (write_per_gene or write_pivots):
        return

    out_dir = f"{out_prefix}_{tag}__per_gene_mod"
    ensure_dir(out_dir)
    prefix_base = os.path.basename(out_prefix)

    per_gene_uns = os.path.join(workdir, f"{tag}.per_gene.uns.tsv")
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

    def key_per_gene(line: str):
        p = line.rstrip("\n").split("\t")
        # 0 gname,1 gid,2 mod,3 chrom,4 start0,5 end0,6 strand,7 zn,8 sample,...
        return (p[0], p[2], p[3], int(p[4]), int(p[5]), p[6], int(p[7]), p[8])

    per_gene_sorted = os.path.join(workdir, f"{tag}.per_gene.sorted.tsv")
    external_sort_tsv(per_gene_uns, per_gene_sorted, key_per_gene, tmpdir=workdir, chunk_lines=chunk_lines, verbose=verbose)

    def basepath(gene_name: str, mod: str) -> str:
        safe_g = sanitize_filename_token(gene_name if gene_name else "NA")
        safe_mod = sanitize_filename_token(str(mod))
        return os.path.join(out_dir, f"{prefix_base}__{safe_g}__{safe_mod}")

    curr_gm = None  # (gene_name, gene_id, mod)
    row_fh = None

    piv_cov = defaultdict(dict)   # idx -> {sample: cov}
    piv_nmod = defaultdict(dict)  # idx -> {sample: nmod}
    piv_frac = defaultdict(dict)  # idx -> {sample: frac}
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
            # idx: (chrom,start0,end0,strand,zn)
            idxs = sorted(piv_cov.keys(), key=lambda t: (t[0], t[1], t[2], t[3], t[4]))

            def write_pivot(path, m, is_float: bool):
                with open(path, "w") as f:
                    f.write("\t".join(["chrom", "start0", "end0", "strand", "ZN_transcript_index"] + samples) + "\n")
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
            cov, nmod, frac = p[9], p[10], p[17]
            ncan, nother, ndel, nfail, ndiff, nnocall = p[11], p[12], p[13], p[14], p[15], p[16]

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
                row_fh.write("\t".join([
                    gname, gid, mod, chrom, start0, end0, strand, zn, sample,
                    cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall, frac
                ]) + "\n")

            if write_pivots:
                samples_seen.add(sample)
                idx = (chrom, int(start0), int(end0), strand, int(zn))
                # first-seen semantics per (idx,sample)
                if sample not in piv_cov[idx]:
                    piv_cov[idx][sample] = int(cov)
                    piv_nmod[idx][sample] = int(nmod)
                    piv_frac[idx][sample] = float(frac)

    flush_group()

    if verbose:
        print(f"[per-gene {tag}] wrote outputs under {out_dir}", file=sys.stderr)

# ----------------------------- Key functions for sorts -----------------------------

def key_norm_for_dedup(line: str):
    p = line.rstrip("\n").split("\t")
    # NORM: 0 sample,1 zn,2 chrom,3 start0,4 end0,5 strand,6 mod,7 gid,8 gname,...
    # Dedup sort key: sample, mod, zn, chrom, start0, end0, strand, gid, gname
    return (p[0], p[6], int(p[1]), p[2], int(p[3]), int(p[4]), p[5], p[7], p[8])

def key_passing_site(line: str):
    p = line.rstrip("\n").split("\t")
    return (p[0], int(p[1]), int(p[2]), p[3], p[4])

def key_dedup_by_site(line: str):
    p = line.rstrip("\n").split("\t")
    # dedup: 0 sample,1 zn,2 chrom,3 start0,4 end0,5 strand,6 mod,7 gid,8 gname,...
    return (p[2], int(p[3]), int(p[4]), p[5], p[6], p[0], int(p[1]), p[7], p[8])

# ----------------------------- Main -----------------------------

def main():
    args = parse_args()

    beds = iter_numbered_beds(args.modkit_dir)
    if not beds:
        sys.exit(f"No numbered ZN partition files found under {args.modkit_dir}")

    gene_index = load_gene_intervals_from_gtf(args.gtf, verbose=args.verbose)

    workdir = tempfile.mkdtemp(prefix=f"aggregate_by_gene_{os.getpid()}_", dir=args.tmpdir)
    if args.verbose:
        print(f"[tmp] workdir={workdir}", file=sys.stderr)

    try:
        # Stage 1: normalize
        norm_tsv = os.path.join(workdir, "norm.tsv")
        normalize_to_tsv(beds, gene_index, norm_tsv, verbose=args.verbose)

        # Stage 2: external sort for dedup reduce
        norm_sorted = os.path.join(workdir, "norm.sorted.tsv")
        external_sort_tsv(
            norm_tsv, norm_sorted,
            key_norm_for_dedup,
            tmpdir=workdir,
            chunk_lines=args.chunk_lines,
            verbose=args.verbose
        )

        base = args.out_prefix

        # Stage 3: dedup reduce (RAW) and collect passing sites (if filter enabled)
        dedup_raw = os.path.join(workdir, "dedup.RAW.tsv")
        passing_sites_uns = os.path.join(workdir, "passing_sites.unsorted.tsv") if args.filter_enable else None

        out_long_raw = f"{base}_RAW_sites_long.tsv" if (args.emit_raw and args.write_long) else None

        if args.emit_raw or args.filter_enable or args.emit_filt:
            dedup_reduce_sorted(
                sorted_norm_tsv=norm_sorted,
                out_dedup_tsv=dedup_raw,
                out_long_tsv=out_long_raw,
                out_passing_sites_tsv=passing_sites_uns,
                min_cov=args.min_cov,
                filter_enable=args.filter_enable,
                count_diff_factor=args.count_diff_factor,
                mod_fail_margin=args.mod_fail_margin,
                verbose=args.verbose
            )

        # RAW stats + per-gene
        if args.emit_raw:
            compute_per_sample_mod_stats_from_dedup(
                dedup_tsv=dedup_raw,
                out_prefix=base,
                tag="RAW",
                workdir=workdir,
                chunk_lines=args.chunk_lines,
                verbose=args.verbose
            )
            generate_per_gene_outputs_from_dedup(
                dedup_tsv=dedup_raw,
                out_prefix=base,
                tag="RAW",
                write_per_gene=args.write_raw_per_gene,
                write_pivots=args.write_pivots,
                workdir=workdir,
                chunk_lines=args.chunk_lines,
                verbose=args.verbose
            )

        # Stage 4: FILTERED subset + stats + per-gene
        if args.emit_filt:
            if args.filter_enable:
                # Sort + unique passing sites by site key
                passing_sorted_tmp = os.path.join(workdir, "passing_sites.sorted.tsv")
                external_sort_tsv(
                    passing_sites_uns,
                    passing_sorted_tmp,
                    key_passing_site,
                    tmpdir=workdir,
                    chunk_lines=args.chunk_lines,
                    verbose=args.verbose
                )
                passing_sorted_unique = os.path.join(workdir, "passing_sites.sorted.unique.tsv")
                uniq_sorted_file(passing_sorted_tmp, passing_sorted_unique, key_passing_site)

                # Sort dedup_raw by site key for streaming join
                dedup_by_site = os.path.join(workdir, "dedup.by_site.sorted.tsv")
                external_sort_tsv(
                    dedup_raw, dedup_by_site,
                    key_dedup_by_site,
                    tmpdir=workdir,
                    chunk_lines=args.chunk_lines,
                    verbose=args.verbose
                )

                dedup_filt = os.path.join(workdir, "dedup.FILTERED.tsv")
                out_long_filt = f"{base}_FILTERED_sites_long.tsv" if args.write_long else None

                filter_dedup_by_passing_sites(
                    dedup_by_site_sorted=dedup_by_site,
                    passing_sites_sorted_unique=passing_sorted_unique,
                    out_dedup_filtered=dedup_filt,
                    out_long_filtered=out_long_filt,
                    verbose=args.verbose
                )
            else:
                # no filtering requested => FILTERED == RAW semantics
                dedup_filt = dedup_raw
                out_long_filt = f"{base}_FILTERED_sites_long.tsv" if args.write_long else None
                if out_long_filt is not None:
                    ensure_dir(os.path.dirname(out_long_filt) or ".")
                    if out_long_raw and os.path.exists(out_long_raw):
                        shutil.copyfile(out_long_raw, out_long_filt)
                    else:
                        # write long directly from dedup
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

            # FILTERED stats + per-gene
            compute_per_sample_mod_stats_from_dedup(
                dedup_tsv=dedup_filt,
                out_prefix=base,
                tag="FILTERED",
                workdir=workdir,
                chunk_lines=args.chunk_lines,
                verbose=args.verbose
            )
            generate_per_gene_outputs_from_dedup(
                dedup_tsv=dedup_filt,
                out_prefix=base,
                tag="FILTERED",
                write_per_gene=args.write_filtered_per_gene,
                write_pivots=args.write_pivots,
                workdir=workdir,
                chunk_lines=args.chunk_lines,
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

