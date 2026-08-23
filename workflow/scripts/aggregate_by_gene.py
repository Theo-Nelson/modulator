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
    - FAIL if Nmod < nfail_score_k * (Nfail + 1)   (NFail-SCORE k-ratio; k may be set per mod_code)
  FILTERED outputs are subset of RAW by site_key (chrom,start0,end0,strand,mod_code).

Stats output (MEANS ONLY; no medians; avoids huge metric sorts):
- <out_prefix>_{RAW|FILTERED}__per_sample_mod_site_stats.tsv
- <out_prefix>_{RAW|FILTERED}__per_sample_mod_tx_stats.tsv
- <out_prefix>_{RAW|FILTERED}__per_tx_mod_stats.tsv

No dependency on GNU sort. Uses:
- chunked in-memory sorting (bounded by --chunk-lines)
- k-way merge across sorted chunks (heapq)

Important fixes included in this version:
- Robust boolean parsing for CLI flags that may be passed as strings (Snakemake sometimes does this).
- Explicit on/off flags supported: --emit-raw / --no-emit-raw, etc.
- Avoids “truthy string” surprises (e.g., "false" treated as True).
- Stats are MEANS ONLY (no median metric sorts).
"""

import os
import sys
import re
import gzip
import argparse
import tempfile
import shutil
import heapq
from concurrent.futures import ProcessPoolExecutor
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


def exon_overlap_len(ex1, ex2):
    total = 0
    i = 0
    j = 0
    ex1 = sorted(ex1)
    ex2 = sorted(ex2)
    while i < len(ex1) and j < len(ex2):
        s1, e1 = ex1[i]
        s2, e2 = ex2[j]
        lo = max(s1, s2)
        hi = min(e1, e2)
        if hi >= lo:
            total += (hi - lo + 1)
        if e1 < e2:
            i += 1
        else:
            j += 1
    return total


def parse_nfail_score_k(spec) -> "tuple[float, dict]":
    """Parse the --nfail-score-k spec into (default_k, per_mod_k).

    Accepts either:
      - a bare number, applied to every mod code:      "1.0"        -> (1.0, {})
      - a per-mod-code map (comma-separated k=v):       "a=0.4,17802=1.0"
        keys are mod codes exactly as they appear in the data (single letters like 'a','m'
        or numeric ChEBI codes like '17802'). Mods NOT listed inherit the fallback, which is the
        standard k=1 guard unless overridden with an explicit 'default=' key -> (1.0, {'a':0.4,'17802':1.0}).
        (Using default_k=1.0 here -- NOT 0.0 -- is critical: a map like 'a=0.4' must still filter the
        other mod codes in the data, not silently disable the confident-call guard for them. To turn the
        guard OFF for unlisted mods, opt in explicitly with 'default=0'.)
      - empty / None:                                   -> (0.0, {})   (filter disabled)

    k is the NFail-SCORE k-ratio calibrated per modification (basecaller + model + version); see
    resources/nfail_score_k_calibration.tsv.
    """
    if spec is None:
        return (0.0, {})
    s = str(spec).strip()
    if s == "":
        return (0.0, {})
    if ("=" not in s) and ("," not in s):
        return (float(s), {})                      # bare scalar -> all mods
    default_k = 1.0                                 # unlisted mods keep the standard guard, not k=0
    per_mod = {}
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        key, sep, val = tok.partition("=")
        if not sep:
            raise ValueError(f"--nfail-score-k token '{tok}' is not of the form mod=k or default=k")
        key = key.strip()
        k = float(val.strip())
        if key.lower() == "default":
            default_k = k
        else:
            per_mod[key] = k
    return (default_k, per_mod)


def resolve_nfail_score_k(mod_code, default_k: float, per_mod: dict) -> float:
    """Look up the k-ratio threshold for one mod code, falling back to default_k."""
    if per_mod:
        return per_mod.get(str(mod_code), default_k)
    return default_k


def row_pass_filter(
    cov: int,
    nmod: int,
    nfail: int,
    ndiff: int,
    count_diff_factor: float,
    nfail_score_k: float = 1.0,
) -> bool:
    # (1) variant/misalignment guard: drop positions dominated by a different base, not a modification
    if ndiff > (count_diff_factor * cov):
        return False
    # (2) NFail-SCORE k-ratio (Nelson et al., "NFail-SCORE") -- the confident-call guard: error-prone
    # false-positive sites carry a large NFail; true sites carry a small one. FAIL if
    # k = Nmod / (NFail + 1) < nfail_score_k. This SUPERSEDES the old "Nmod > Nfail + mod_fail_margin"
    # rule -- k=1 is equivalent to margin=0 (Nmod must strictly exceed Nfail). k is calibrated per
    # modification (basecaller + model + Dorado version; see resources/nfail_score_k_calibration.tsv)
    # and can be set per mod_code (parse_nfail_score_k / resolve_nfail_score_k). Disabled when k <= 0.
    if nfail_score_k > 0.0 and nmod < nfail_score_k * (nfail + 1):
        return False
    return True


def parse_bool(x, default: bool) -> bool:
    """
    Robust bool parsing:
    - If x is already bool => return it
    - If x is None => default
    - If x is int/float => bool(x)
    - If x is str => parse common true/false tokens
    """
    if isinstance(x, bool):
        return x
    if x is None:
        return default
    if isinstance(x, (int, float)):
        return bool(x)
    if isinstance(x, str):
        s = x.strip().lower()
        if s in ("1", "true", "t", "yes", "y", "on"):
            return True
        if s in ("0", "false", "f", "no", "n", "off", ""):
            return False
    # fallback: python truthiness
    return bool(x)


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
    ap.add_argument("--mod-fail-margin", type=int, default=1,
                    help="DEPRECATED / no-op: the confident-call guard is now the NFail-SCORE k-ratio "
                         "(--nfail-score-k); k=1 reproduces the old margin=0 behaviour. Accepted but ignored.")
    ap.add_argument("--nfail-score-k", type=str, default="1.0",
                    help="NFail-SCORE k-ratio confident-call filter: FAIL if Nmod < k*(Nfail+1), i.e. "
                         "k = Nmod/(Nfail+1) < this. Either a single value for all mods (e.g. '1.0') or a "
                         "per-mod-code map (e.g. 'a=0.4,17802=1.0,default=1.0'). Calibrate k per modification "
                         "basecaller + version (see resources/nfail_score_k_calibration.tsv). 0 disables.")

    # output toggles (explicit on/off)
    ap.add_argument("--emit-raw", dest="emit_raw", action="store_true")
    ap.add_argument("--no-emit-raw", dest="emit_raw", action="store_false")
    ap.set_defaults(emit_raw=True)

    ap.add_argument("--emit-filtered", dest="emit_filt", action="store_true")
    ap.add_argument("--no-emit-filtered", dest="emit_filt", action="store_false")
    ap.set_defaults(emit_filt=True)

    ap.add_argument("--write-long", dest="write_long", action="store_true")
    ap.add_argument("--no-write-long", dest="write_long", action="store_false")
    ap.set_defaults(write_long=True)

    # Per-gene pivots are optional inspection outputs (3 dense files per gene x mod group);
    # nothing downstream reads them. 'auto' writes them unless the run would explode into too
    # many tiny files (see --pivot-max-groups); 'on' always writes them (even at scale); 'off'
    # never does. The legacy --write-pivots/--no-write-pivots map onto 'on'/'off'.
    ap.add_argument("--pivot-mode", dest="pivot_mode", choices=["auto", "on", "off"],
                    default="auto",
                    help="auto (default): write per-gene pivots unless (gene x mod) groups exceed "
                         "--pivot-max-groups; on: always write; off: never write.")
    ap.add_argument("--pivot-max-groups", dest="pivot_max_groups", type=int, default=2000,
                    help="auto-mode ceiling on the number of (gene x mod) pivot groups; each group "
                         "writes 3 files, so this bounds the small-file count at ~3x this value.")
    ap.add_argument("--write-pivots", dest="pivot_mode", action="store_const", const="on",
                    help="(legacy alias for --pivot-mode on)")
    ap.add_argument("--no-write-pivots", dest="pivot_mode", action="store_const", const="off",
                    help="(legacy alias for --pivot-mode off)")

    ap.add_argument("--write-raw-per-gene", dest="write_raw_per_gene", action="store_true")
    ap.add_argument("--no-write-raw-per-gene", dest="write_raw_per_gene", action="store_false")
    ap.set_defaults(write_raw_per_gene=False)

    ap.add_argument("--write-filtered-per-gene", dest="write_filtered_per_gene", action="store_true")
    ap.add_argument("--no-write-filtered-per-gene", dest="write_filtered_per_gene", action="store_false")
    ap.set_defaults(write_filtered_per_gene=True)

    ap.add_argument("--verbose", action="store_true")

    # Pure-Python external sort tuning
    ap.add_argument(
        "--chunk-lines",
        type=int,
        default=2_000_000,
        help="Lines per in-memory chunk during external sort (default: 2,000,000)"
    )
    ap.add_argument(
        "--tmpdir",
        type=str,
        default=os.environ.get("TMPDIR", "/tmp"),
        help="Temp directory for intermediates (default: $TMPDIR or /tmp)"
    )
    ap.add_argument("--keep-intermediates", action="store_true", help="Do not delete intermediates (debugging)")

    args = ap.parse_args()

    # Defensive: in case something upstream passes string-ish values into args
    # (argparse itself normally produces bools here, but this prevents surprises).
    args.emit_raw = parse_bool(args.emit_raw, default=True)
    args.emit_filt = parse_bool(args.emit_filt, default=True)
    args.write_long = parse_bool(args.write_long, default=True)
    args.write_raw_per_gene = parse_bool(args.write_raw_per_gene, default=False)
    args.write_filtered_per_gene = parse_bool(args.write_filtered_per_gene, default=True)

    return args


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

    # Isolate this sort's chunk files in a private subdir. Two external sorts that share the
    # same --tmpdir -- concurrent aggregate_by_gene invocations (jobs>1), or the RAW and
    # FILTERED per-gene sorts within one run -- would otherwise both write chunk_00000.tsv,
    # chunk_00001.tsv, ... into tmpdir and silently corrupt each other's intermediates.
    work = tempfile.mkdtemp(prefix="agg_sort_", dir=tmpdir)
    try:
        _external_sort_into(in_path, out_path, key_func, work, chunk_lines, verbose)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _external_sort_into(in_path, out_path, key_func, work, chunk_lines, verbose):
    chunks: List[str] = []
    buf: List[str] = []
    n_in = 0

    def write_chunk(lines: List[str], idx: int) -> str:
        lines.sort(key=key_func)
        cpath = os.path.join(work, f"chunk_{idx:05d}.tsv")
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

    if not chunks:
        # empty input
        with open(out_path, "w"):
            pass
        return

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

TxInterval = namedtuple("TxInterval", ["start", "end", "gene_id", "gene_name", "strand", "zn", "exons"])
GeneInterval = namedtuple("GeneInterval", ["start", "end", "gene_id", "gene_name", "strand", "exons"])


def load_gene_intervals_from_gtf(gtf_path: str, verbose=False):
    """
    Build transcript-aware interval indices from the assembler GTF.
    Uses transcript-level zn_index to map numbered ZN partitions back to gene labels.
    """
    tx_meta = {}
    gene_exons = defaultdict(list)
    gene_name_map = {}

    with open_text(gtf_path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = parts
            a = {}
            for kv in re.finditer(r'(\S+)\s+"([^"]*)"', attrs):
                a[kv.group(1)] = kv.group(2)

            gene_id = a.get("gene_id") or a.get("gtf_gene_id") or a.get("gene") or ""
            gene_name = a.get("ref_gene_name") or a.get("gene_name") or a.get("gtf_gene_name") or gene_id
            if not gene_id:
                continue

            s = int(start)
            e = int(end)
            if feature == "transcript":
                tid = a.get("transcript_id")
                if not tid:
                    continue
                zn = safe_int(a.get("zn_index", a.get("transcript_index", 0)))
                tx_meta[tid] = dict(
                    chrom=chrom,
                    strand=strand,
                    start=s,
                    end=e,
                    gene_id=gene_id,
                    gene_name=gene_name,
                    zn=zn,
                    exons=[],
                )
            elif feature == "exon":
                tid = a.get("transcript_id")
                if tid and tid in tx_meta:
                    tx_meta[tid]["exons"].append((s, e))
                gene_exons[(chrom, strand, gene_id)].append((s, e))
                gene_name_map[gene_id] = gene_name

    tx_by_cs = defaultdict(list)
    for meta in tx_meta.values():
        exons = sorted(meta["exons"]) if meta["exons"] else [(meta["start"], meta["end"])]
        tx_by_cs[(meta["chrom"], meta["strand"])].append(
            TxInterval(meta["start"], meta["end"], meta["gene_id"], meta["gene_name"], meta["strand"], meta["zn"], exons)
        )

    gene_by_cs = defaultdict(list)
    for (chrom, strand, gid), exons in gene_exons.items():
        exons = sorted(exons)
        merged = []
        for s, e in exons:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        merged_t = [(s, e) for s, e in merged]
        gene_by_cs[(chrom, strand)].append(
            GeneInterval(merged_t[0][0], merged_t[-1][1], gid, gene_name_map.get(gid, gid), strand, merged_t)
        )

    for k in tx_by_cs:
        tx_by_cs[k].sort(key=lambda iv: (iv.start, iv.end, iv.zn, iv.gene_id))
    for k in gene_by_cs:
        gene_by_cs[k].sort(key=lambda iv: (iv.start, iv.end, iv.gene_id))

    if verbose:
        total_tx = sum(len(v) for v in tx_by_cs.values())
        total_genes = sum(len(v) for v in gene_by_cs.values())
        print(f"[info] loaded {total_tx} transcript partitions and {total_genes} gene spans from {gtf_path}", file=sys.stderr)

    return tx_by_cs, gene_by_cs


def site_interval_1based(pos_start: int, pos_end: int):
    s = int(pos_start) + 1
    e = max(s, int(pos_end))
    return [(s, e)]


def assign_gene(
    chrom: str,
    pos_start: int,
    pos_end: int,
    strand: str,
    zn: int,
    tx_index,
    gene_index,
) -> Tuple[str, str]:
    """
    Return (gene_id, gene_name) using same-strand transcript partitions first,
    keyed by zn_index and exonic overlap. Falls back to coarse gene overlap.
    """
    site_exon = site_interval_1based(pos_start, pos_end)

    best = None
    best_ov = -1
    for iv in tx_index.get((chrom, strand), []):
        if iv.start > pos_end:
            break
        if iv.end < (pos_start + 1):
            continue
        if int(iv.zn) != int(zn):
            continue
        ov = exon_overlap_len(site_exon, iv.exons)
        if ov > best_ov:
            best_ov = ov
            best = iv
    if best and best_ov > 0:
        return best.gene_id, best.gene_name

    best = None
    best_ov = -1
    for iv in gene_index.get((chrom, strand), []):
        if iv.start > pos_end:
            break
        if iv.end < (pos_start + 1):
            continue
        ov = exon_overlap_len(site_exon, iv.exons)
        if ov > best_ov:
            best_ov = ov
            best = iv
    if best and best_ov > 0:
        return best.gene_id, best.gene_name

    other = "+" if strand == "-" else "-"
    for iv in gene_index.get((chrom, other), []):
        if iv.start > pos_end:
            break
        if iv.end < (pos_start + 1):
            continue
        ov = exon_overlap_len(site_exon, iv.exons)
        if ov > 0:
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
    # dedup by (root, sample, ZN): if BOTH N.bed and N.bed.gz exist in a dir, the partition would
    # otherwise be read TWICE (silently doubling Nvalid_cov/Nmod/Nfail). Keep one, preferring .gz.
    picked = {}
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
            key = (root, sample_name, zn)
            prev = picked.get(key)
            if prev is None or (fname.endswith(".gz") and not prev.endswith(".gz")):
                picked[key] = os.path.join(root, fname)
    out = [(root, sample_name, path, zn) for (root, sample_name, zn), path in picked.items()]
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


def normalize_to_tsv(
    beds: List[Tuple[str, str, str, int]],
    tx_index,
    gene_index,
    out_tsv: str,
    verbose: bool = False
) -> int:
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

                    gid, gname = assign_gene(chrom, start0, end0, strand, zn, tx_index, gene_index)

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
    k_default: float = 1.0,
    k_per_mod: Optional[dict] = None,
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
                f"{frac:.6f}",
            ]) + "\n"
        )

        if long_fh is not None:
            long_fh.write(
                f"{sample}\t{zn}\t{chrom}\t{start0}\t{end0}\t{strand}\t{mod}\t"
                f"{cov}\t{nmod}\t{frac:.6f}\t{gid}\t{gname}\t"
                f"{ncan}\t{nother}\t{ndel}\t{nfail}\t{ndiff}\t{nnocall}\n"
            )

        if pass_fh is not None:
            k = resolve_nfail_score_k(mod, k_default, k_per_mod or {})
            if row_pass_filter(cov, nmod, nfail, ndiff, count_diff_factor, k):
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
            if len(p) < 17:
                continue
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
        return (parts[2], int(parts[3]), int(parts[4]), parts[5], parts[6])

    ps = open(passing_sites_sorted_unique, "r")
    try:
        ps_line = ps.readline()
        ps_key = None
        if ps_line:
            p = ps_line.rstrip("\n").split("\t")
            if len(p) >= 5:
                ps_key = (p[0], int(p[1]), int(p[2]), p[3], p[4])

        kept = 0
        with open(dedup_by_site_sorted, "r") as d:
            for ln in d:
                ln = ln.rstrip("\n")
                if not ln:
                    continue
                parts = ln.split("\t")
                if len(parts) < 18:
                    continue
                dk = site_of_dedup(parts)

                while ps_key is not None and ps_key < dk:
                    ps_line = ps.readline()
                    if not ps_line:
                        ps_key = None
                        break
                    p = ps_line.rstrip("\n").split("\t")
                    if len(p) < 5:
                        continue
                    ps_key = (p[0], int(p[1]), int(p[2]), p[3], p[4])

                if ps_key is None:
                    break

                if dk == ps_key:
                    out_fh.write("\t".join(parts) + "\n")
                    if long_fh is not None:
                        sample, zn, chrom, start0, end0, strand, mod, gid, gname = (
                            parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], parts[7], parts[8]
                        )
                        cov, nmod = parts[9], parts[10]
                        ncan, nother, ndel, nfail, ndiff, nnocall = (
                            parts[11], parts[12], parts[13], parts[14], parts[15], parts[16]
                        )
                        frac = parts[17]
                        long_fh.write(
                            f"{sample}\t{zn}\t{chrom}\t{start0}\t{end0}\t{strand}\t{mod}\t"
                            f"{cov}\t{nmod}\t{frac}\t{gid}\t{gname}\t"
                            f"{ncan}\t{nother}\t{ndel}\t{nfail}\t{ndiff}\t{nnocall}\n"
                        )
                    kept += 1

    finally:
        ps.close()
        out_fh.close()
        if long_fh:
            long_fh.close()

    if verbose:
        print(f"[filter] kept {kept} rows -> {out_dedup_filtered}", file=sys.stderr)
    return kept


# ----------------------------- Stats computation (MEANS ONLY) -----------------------------


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

    Notes:
    - Reports MEANS only (no medians).
    - Avoids generating gigantic metric files and sorting them.
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
    external_sort_tsv(
        site_sample_uns,
        site_sample_sorted,
        key_site_sample,
        tmpdir=workdir,
        chunk_lines=chunk_lines,
        verbose=verbose,
    )

    # Reduce per-site -> per (sample,mod) aggregates needed for means
    n_sites_total_by_sm = defaultdict(int)
    n_sites_detected_by_sm = defaultdict(int)
    total_nmod_by_sm = defaultdict(int)
    total_cov_by_sm = defaultdict(int)

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

    with open(site_sample_sorted, "r") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            key = (p[0], p[1], p[2], p[3], p[4], p[5])
            nmod = safe_int(p[6])
            cov = safe_int(p[7])

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

    # Write out1 (means only)
    with open(out1, "w") as f:
        hdr = [
            "sample", "mod_code",
            "n_sites_total", "n_sites_detected",
            "total_Nmod", "total_cov", "overall_stoich",
            "mean_site_cov", "mean_site_Nmod",
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
            rows.append((mod, sample, n_sites, nd, tn, tc, overall, mean_cov, mean_nm))

        rows.sort(key=lambda x: (x[0], x[1]))
        for mod, sample, nst, nsd, tn, tc, ov, mc, mnm in rows:
            f.write(
                f"{sample}\t{mod}\t{nst}\t{nsd}\t{tn}\t{tc}\t{ov:.6f}\t{mc:.6f}\t{mnm:.6f}\n"
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
                sample = p[0]
                zn = p[1]
                chrom, start0, end0, strand, mod = p[2], p[3], p[4], p[5], p[6]
                cov, nmod = p[9], p[10]
                out.write("\t".join([sample, mod, zn, chrom, start0, end0, strand, nmod, cov]) + "\n")

    def key_site_tx(line: str):
        p = line.rstrip("\n").split("\t")
        return (p[0], p[1], int(p[2]), p[3], int(p[4]), int(p[5]), p[6])

    site_tx_sorted = os.path.join(workdir, f"{tag}.site_tx.sorted.tsv")
    external_sort_tsv(
        site_tx_uns,
        site_tx_sorted,
        key_site_tx,
        tmpdir=workdir,
        chunk_lines=chunk_lines,
        verbose=verbose,
    )

    tx_set = defaultdict(set)
    sum_det_sites = defaultdict(int)
    sum_total_nmod_per_tx = defaultdict(int)
    sum_tx_sto = defaultdict(float)

    with open(out3, "w") as out:
        out.write("\t".join([
            "sample", "mod_code", "ZN_transcript_index",
            "n_sites_total", "n_sites_detected",
            "total_Nmod", "total_cov", "tx_stoich"
        ]) + "\n")

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
            out.write(
                f"{sample}\t{mod}\t{zn}\t{tx_n_sites_total}\t{tx_n_sites_detected}\t"
                f"{tx_total_nmod}\t{tx_total_cov}\t{sto:.6f}\n"
            )

            sm = (sample, mod)
            tx_set[sm].add(int(zn))
            sum_det_sites[sm] += tx_n_sites_detected
            sum_total_nmod_per_tx[sm] += tx_total_nmod
            sum_tx_sto[sm] += sto

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
                nmod = safe_int(p[7])
                cov = safe_int(p[8])

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

    with open(out2, "w") as f:
        hdr = [
            "sample", "mod_code", "n_tx",
            "mean_detected_sites_per_tx",
            "mean_total_Nmod_per_tx",
            "mean_tx_stoich",
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
            rows.append((mod, sample, n_tx, mean_det, mean_nm, mean_st))

        rows.sort(key=lambda x: (x[0], x[1]))
        for mod, sample, n_tx, md, mn, ms in rows:
            f.write(f"{sample}\t{mod}\t{n_tx}\t{md:.6f}\t{mn:.6f}\t{ms:.6f}\n")

    if verbose:
        print(f"[stats {tag}] wrote {out1}", file=sys.stderr)
        print(f"[stats {tag}] wrote {out2}", file=sys.stderr)
        print(f"[stats {tag}] wrote {out3}", file=sys.stderr)


# ----------------------------- Per-gene outputs + pivots -----------------------------


def _write_one_gene_group(lines, gene_name, gene_id, mod, out_dir, prefix_base,
                          write_per_gene, write_pivots):
    """Write the per-gene row TSV and/or the 3 pivots for ONE (gene_name, gene_id, mod) group.
    `lines` are the already-grouped rows (18 tab-cols each) from the sorted per_gene file. Each
    group writes to its own basepath files, so this is safe to run in parallel across groups."""
    # The filename must be INJECTIVE over the (gene_name, gene_id, mod) group key. Including gene_id
    # is not sufficient on its own: sanitize_filename_token maps every char outside [A-Za-z0-9._+-] to
    # '_', so two genuinely distinct keys that differ only in a sanitized character collide on one path
    # -- and "w" mode then truncates (serial) or two workers tear the file (jobs>1). Append a short hash
    # of the RAW key so the path is unique even when the human-readable part collides.
    import hashlib
    _keyhash = hashlib.sha1(f"{gene_name}\x00{gene_id}\x00{mod}".encode("utf-8", "replace")).hexdigest()[:8]
    safe_g = sanitize_filename_token(gene_name if gene_name else "NA")
    safe_gid = sanitize_filename_token(gene_id if gene_id else "NA")
    safe_mod = sanitize_filename_token(str(mod))
    bp = os.path.join(out_dir, f"{prefix_base}__{safe_g}__{safe_gid}__{safe_mod}__{_keyhash}")

    row_fh = open(f"{bp}.tsv", "w") if write_per_gene else None
    if row_fh is not None:
        row_fh.write("\t".join(PER_GENE_COLS) + "\n")

    piv_cov = defaultdict(dict)
    piv_nmod = defaultdict(dict)
    piv_frac = defaultdict(dict)
    samples_seen = set()

    for ln in lines:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 18:
            continue
        gname, gid, m = p[0], p[1], p[2]
        chrom, start0, end0, strand, zn, sample = p[3], p[4], p[5], p[6], p[7], p[8]
        cov, nmod, frac = p[9], p[10], p[17]
        ncan, nother, ndel, nfail, ndiff, nnocall = p[11], p[12], p[13], p[14], p[15], p[16]
        if row_fh is not None:
            row_fh.write("\t".join([
                gname, gid, m, chrom, start0, end0, strand, zn, sample,
                cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall, frac
            ]) + "\n")
        if write_pivots:
            samples_seen.add(sample)
            idx = (chrom, int(start0), int(end0), strand, int(zn))
            if sample not in piv_cov[idx]:
                piv_cov[idx][sample] = int(cov)
                piv_nmod[idx][sample] = int(nmod)
                piv_frac[idx][sample] = float(frac)

    if row_fh is not None:
        row_fh.close()

    if write_pivots:
        samples = sorted(samples_seen)
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
                        row.append(f"{float(v):.6f}" if is_float else str(int(v)))
                    f.write("\t".join(row) + "\n")

        write_pivot(f"{bp}_cov_pivot.tsv", piv_cov, is_float=False)
        write_pivot(f"{bp}_Nmod_pivot.tsv", piv_nmod, is_float=False)
        write_pivot(f"{bp}_frac_pivot.tsv", piv_frac, is_float=True)


def _emit_gene_group_worker(task):
    """Picklable ProcessPool entry point: read one group's byte range from the sorted file and
    write its outputs. Reading one gene's rows keeps per-worker memory to a single gene."""
    (sorted_path, start_off, end_off, gene_name, gene_id, mod,
     out_dir, prefix_base, write_per_gene, write_pivots) = task
    lines = []
    with open(sorted_path, "rb") as f:
        f.seek(start_off)
        while f.tell() < end_off:
            raw = f.readline()
            if not raw:
                break
            lines.append(raw.decode("utf-8", "replace"))
    _write_one_gene_group(lines, gene_name, gene_id, mod, out_dir, prefix_base,
                          write_per_gene, write_pivots)
    return 1


def generate_per_gene_outputs_from_dedup(
    dedup_tsv: str,
    out_prefix: str,
    tag: str,
    write_per_gene: bool,
    pivot_mode="auto",
    workdir: str = None,
    chunk_lines: int = 2_000_000,
    verbose: bool = False,
    jobs: int = 1,
    pivot_max_groups: int = 2000,
):
    """
    Writes per-gene×mod tables and pivots under:
      <out_prefix>_<TAG>__per_gene_mod/

    - row TSV per gene/mod: <prefix_base>__<gene>__<mod>.tsv
    - pivots: *_cov_pivot.tsv, *_Nmod_pivot.tsv, *_frac_pivot.tsv

    Per-gene groups are independent (distinct output files), so with jobs>1 they are written
    concurrently across a ProcessPool (each worker reads one gene's byte range from the sorted
    file -> bounded memory). jobs<=1 keeps a plain serial pass (single-core fallback). Output is
    identical to the serial path.
    """
    write_per_gene = parse_bool(write_per_gene, default=False)
    # Normalize the pivot mode. Accept the tri-state strings plus legacy bools/"true"/"false".
    mode = str(pivot_mode).strip().lower()
    if mode not in ("auto", "on", "off"):
        mode = "on" if parse_bool(pivot_mode, default=True) else "off"
    # If pivots are explicitly off and no per-gene tables are requested, there is nothing to do --
    # skip the (potentially large) external sort entirely. In 'auto' we cannot decide yet (the
    # group count is only known after the sort), so we fall through.
    if mode == "off" and not write_per_gene:
        return

    out_dir = f"{out_prefix}_{tag}__per_gene_mod"
    ensure_dir(out_dir)
    # Prune stale per-gene files from a previous run of the same prefix (e.g. a looser filter): the
    # long table is authoritative and is rewritten, but orphaned per-gene files would survive and make
    # a per-gene-view Nmod total exceed the long table. Clear the dir so it only ever holds THIS run.
    for _stale in os.listdir(out_dir):
        _sp = os.path.join(out_dir, _stale)
        if os.path.isfile(_sp):
            try:
                os.remove(_sp)
            except OSError:
                pass
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
                cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall, frac = (
                    p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17]
                )
                out.write("\t".join([
                    gname, gid, mod, chrom, start0, end0, strand, zn, sample,
                    cov, nmod, ncan, nother, ndel, nfail, ndiff, nnocall, frac
                ]) + "\n")

    def key_per_gene(line: str):
        # MUST include gene_id (p[1]): the group boundaries below split on (gene_name, gene_id, mod),
        # so if two gene_ids share a gene_name and the sort ignored gene_id, their rows interleave and
        # form multiple non-contiguous groups that then collide on one gene_name-only filename.
        p = line.rstrip("\n").split("\t")
        return (p[0], p[1], p[2], p[3], int(p[4]), int(p[5]), p[6], int(p[7]), p[8])

    per_gene_sorted = os.path.join(workdir, f"{tag}.per_gene.sorted.tsv")
    external_sort_tsv(per_gene_uns, per_gene_sorted, key_per_gene, tmpdir=workdir, chunk_lines=chunk_lines, verbose=verbose)

    # Index the (gene_name, gene_id, mod) group byte-ranges in the sorted file (one cheap serial
    # pass), then emit each group -- serially (jobs<=1) or across a ProcessPool (jobs>1).
    groups = []  # (gname, gid, mod, start_off, end_off)
    with open(per_gene_sorted, "rb") as f:
        cur = None
        start = 0
        while True:
            off = f.tell()
            raw = f.readline()
            if not raw:
                if cur is not None:
                    groups.append((cur[0], cur[1], cur[2], start, off))
                break
            p = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
            if len(p) < 18:
                continue
            gm = (p[0], p[1], p[2])
            if cur is None:
                cur = gm
                start = off
            elif gm != cur:
                groups.append((cur[0], cur[1], cur[2], start, off))
                cur = gm
                start = off

    # Resolve the pivot decision now that the exact (gene x mod) group count is known.
    n_groups = len(groups)
    if mode == "on":
        write_pivots = True
    elif mode == "off":
        write_pivots = False
    else:  # auto
        write_pivots = n_groups <= int(pivot_max_groups)
        if verbose:
            if write_pivots:
                print(f"[per-gene {tag}] pivots: auto -> ON ({n_groups} gene x mod groups "
                      f"<= pivot_max_groups={pivot_max_groups})", file=sys.stderr)
            else:
                print(f"[per-gene {tag}] pivots: auto -> OFF ({n_groups} gene x mod groups "
                      f"> pivot_max_groups={pivot_max_groups}; would write ~{3 * n_groups} small "
                      f"files). Long/per-gene tables are unaffected; set aggregation.zn.write_pivots"
                      f"=on (or raise pivot_max_groups) to force pivots at this scale.",
                      file=sys.stderr)

    if not (write_per_gene or write_pivots):
        # auto disabled pivots and no per-gene tables requested -> nothing left to emit.
        return

    tasks = [
        (per_gene_sorted, s_off, e_off, gname, gid, mod, out_dir, prefix_base,
         write_per_gene, write_pivots)
        for (gname, gid, mod, s_off, e_off) in groups
    ]

    jobs = max(1, int(jobs))
    if jobs <= 1 or len(tasks) <= 1:
        for t in tasks:
            _emit_gene_group_worker(t)
    else:
        try:
            with ProcessPoolExecutor(max_workers=min(jobs, len(tasks))) as ex:
                for _ in ex.map(_emit_gene_group_worker, tasks):
                    pass
        except Exception as exc:
            if verbose:
                print(f"[per-gene {tag}] parallel emit failed ({exc}); serial fallback",
                      file=sys.stderr)
            for t in tasks:
                _emit_gene_group_worker(t)

    if verbose:
        print(f"[per-gene {tag}] wrote outputs for {len(groups)} gene-group(s) under "
              f"{out_dir} (jobs={jobs})", file=sys.stderr)


# ----------------------------- Key functions for sorts -----------------------------


def key_norm_for_dedup(line: str):
    p = line.rstrip("\n").split("\t")
    return (p[0], p[6], int(p[1]), p[2], int(p[3]), int(p[4]), p[5], p[7], p[8])


def key_passing_site(line: str):
    p = line.rstrip("\n").split("\t")
    return (p[0], int(p[1]), int(p[2]), p[3], p[4])


def key_dedup_by_site(line: str):
    p = line.rstrip("\n").split("\t")
    return (p[2], int(p[3]), int(p[4]), p[5], p[6], p[0], int(p[1]), p[7], p[8])


# ----------------------------- Main -----------------------------


def main():
    args = parse_args()
    _k_default, _k_per_mod = parse_nfail_score_k(args.nfail_score_k)

    beds = iter_numbered_beds(args.modkit_dir)
    if not beds:
        sys.exit(f"No numbered ZN partition files found under {args.modkit_dir}")

    tx_index, gene_index = load_gene_intervals_from_gtf(args.gtf, verbose=args.verbose)

    workdir = tempfile.mkdtemp(prefix=f"aggregate_by_gene_{os.getpid()}_", dir=args.tmpdir)
    if args.verbose:
        print(f"[tmp] workdir={workdir}", file=sys.stderr)
        print(
            "[cfg] "
            f"emit_raw={args.emit_raw} emit_filtered={args.emit_filt} "
            f"write_long={args.write_long} pivot_mode={args.pivot_mode} "
            f"write_raw_per_gene={args.write_raw_per_gene} write_filtered_per_gene={args.write_filtered_per_gene} "
            f"filter_enable={args.filter_enable}",
            file=sys.stderr
        )

    try:
        # Stage 1: normalize
        norm_tsv = os.path.join(workdir, "norm.tsv")
        normalize_to_tsv(beds, tx_index, gene_index, norm_tsv, verbose=args.verbose)

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
                k_default=_k_default,
                k_per_mod=_k_per_mod,
                verbose=args.verbose
            )

        # RAW stats + per-gene outputs
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
                pivot_mode=args.pivot_mode,
                pivot_max_groups=args.pivot_max_groups,
                workdir=workdir,
                chunk_lines=args.chunk_lines,
                verbose=args.verbose
            )

        # FILTERED subset + stats + per-gene
        if args.emit_filt:
            if args.filter_enable:
                if passing_sites_uns is None:
                    sys.exit("Internal error: filter_enable true but passing_sites_uns is None")

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
                        with open(out_long_filt, "w") as out:
                            out.write("\t".join(LONG_HEADER) + "\n")
                            with open(dedup_raw, "r") as f:
                                for ln in f:
                                    p = ln.rstrip("\n").split("\t")
                                    if len(p) < 18:
                                        continue
                                    sample, zn, chrom, start0, end0, strand, mod, gid, gname = (
                                        p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]
                                    )
                                    cov, nmod = p[9], p[10]
                                    ncan, nother, ndel, nfail, ndiff, nnocall, frac = (
                                        p[11], p[12], p[13], p[14], p[15], p[16], p[17]
                                    )
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
                pivot_mode=args.pivot_mode,
                pivot_max_groups=args.pivot_max_groups,
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
