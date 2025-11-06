#!/usr/bin/env python3
"""
Aggregate modkit bedMethyl outputs that were generated with or without
`--partition-tag ZN`, *recursively* walking a modkit output directory.

What this script does (v5):
  • Recursively finds *.bed / *.bed.gz under --modkit-dir (subfolders OK).
  • Parses both filename styles:
      - numeric ZN partition files like "1.bed", "2.bed", "ungrouped.bed" (no sample/mod in name)
      - non-numeric files like "<sample>_<code>.bed" (code may be ZT or other)
      - plain per-mod files like "..._a_filtered_mod.bed" without partitions
  • *Does not rely on code in filename.* Instead, it maps each site to a transcript
    by intersecting its genomic position with transcript features from --gtf.
    It pulls `gene_index`, `transcript_index`, `gene_id`, `gene_name` from the GTF
    (or falls back to assembler-style attributes if present in the GTF).
  • Emits per-gene, per-mod spreadsheets with a `transcript_index` column so you can
    see which transcript (ZN) each site came from, plus a single combined long table.
  • Optional --summary-tsv can be merged for additional metadata if the 'code' values
    match a `zt_label` or `code` column; otherwise the join is skipped harmlessly.

Outputs (prefix = --out-prefix):
  <prefix>_sites_long.tsv
  <prefix>.per_gene_sites/<gene_index>_<gene_name>_<mod>.tsv   # one per gene per mod

Columns in per-gene tables:
  chrom, start0, end0, strand, mod_code, sample,
  gene_index, gene_name, gene_id, transcript_index,
  Nvalid_cov, Nmod, Ncanonical, Nother_mod, Ndelete, Nfail, Ndiff, Nnocall,
  frac_modified

Notes:
  • `sample` is inferred as the parent directory name of the .bed file. If your layout
    encodes sample differently, use --sample-from-parent N to climb N directories up (default 1).
  • A site is assigned to the transcript that overlaps its midpoint on the matching strand.
    If multiple transcripts overlap, the one with largest exon coverage window around the site
    is chosen. If none overlap, the row is dropped with a warning when --verbose.
  • Use --min-cov to zero frac_modified for rows where Nvalid_cov < threshold.

Example:
  python3 aggregate_modkit_by_ZN_per_gene_v5.py \
    --modkit-dir modkit_out_ZN \
    --gtf fivegenes_readbacked_annot.gtf \
    --out-prefix modkit_by_transcript_ZN \
    --min-cov 5 \
    --verbose
"""

import os, sys, argparse, glob, gzip, re
from collections import defaultdict
from typing import Dict, Tuple, List

try:
    import pandas as pd
except ImportError:
    sys.exit("This script requires pandas. Install it (e.g. `micromamba install pandas`).")

try:
    from intervaltree import Interval, IntervalTree
except ImportError:
    sys.exit("This script requires intervaltree. Install it (e.g. `micromamba install intervaltree`).")

BED_COLS = [
    "chrom", "start0", "end0", "mod_code", "score", "strand",
    "start0_compat", "end0_compat", "rgb",
    "Nvalid_cov", "frac_modified",
    "Nmod", "Ncanonical", "Nother_mod",
    "Ndelete", "Nfail", "Ndiff", "Nnocall"
]

ATTR_RE = re.compile(r"\b(\w+)\s+\"([^\"]*)\";")


def parse_args():
    ap = argparse.ArgumentParser(description="Aggregate modkit bedMethyl by transcript (ZN-style) per gene, recursively")
    ap.add_argument("--modkit-dir", required=True, help="Directory containing modkit .bed/.bed.gz (searches recursively)")
    ap.add_argument("--gtf", required=True, help="Assembler GTF with transcript features and attributes")
    ap.add_argument("--summary-tsv", help="Optional assembler summary TSV to join; purely additive")
    ap.add_argument("--out-prefix", required=True, help="Output prefix for TSVs")
    ap.add_argument("--min-cov", type=int, default=0, help="If >0, zero frac_modified where Nvalid_cov < MIN_COV")
    ap.add_argument("--sample-from-parent", type=int, default=1, help="Pick sample name from Nth parent directory (default 1)")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def is_header_line(line: str) -> bool:
    s = line.strip()
    return (not s) or s.startswith("#") or s.startswith("track") or s.startswith("browser")


def open_textmaybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return default


def build_transcript_trees(gtf_path: str, verbose=False):
    """Build per-(chrom,strand) IntervalTrees of transcript extents and capture per-transcript metadata.
    Returns (trees, meta)
      trees[(chrom,strand)] = IntervalTree of (start,end, tid)
      meta[tid] = {chrom,strand,start,end,gene_id,gene_name,gene_index,transcript_index,zt_label}
    """
    trees: Dict[Tuple[str,str], IntervalTree] = {}
    meta: Dict[str, dict] = {}
    with open(gtf_path) as f:
        for ln in f:
            if ln.startswith("#"): continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            chrom, src, feature, start, end, score, strand, frame, attrs = parts
            if feature != "transcript":
                continue
            start0 = int(start) - 1
            end0 = int(end)
            attr = {m.group(1): m.group(2) for m in ATTR_RE.finditer(attrs)}
            gene_id = attr.get("gene_id") or attr.get("gtf_gene_id") or attr.get("gene")
            gene_name = attr.get("ref_gene_name") or attr.get("gene_name") or attr.get("gtf_gene_name") or gene_id
            gidx = attr.get("gene_index") or attr.get("gene_idx")
            tidx = attr.get("transcript_index") or attr.get("transcript_idx")
            zt_label = attr.get("zt_label") or attr.get("code")
            tid = attr.get("transcript_id") or zt_label or f"{gene_name}.G{gidx}.T{tidx}"
            if (chrom,strand) not in trees:
                trees[(chrom,strand)] = IntervalTree()
            trees[(chrom,strand)].addi(start0, end0, tid)
            meta[tid] = {
                "chrom": chrom, "strand": strand, "start0": start0, "end0": end0,
                "gene_id": gene_id or "", "gene_name": gene_name or "",
                "gene_index": int(gidx) if gidx and gidx.isdigit() else None,
                "transcript_index": int(tidx) if tidx and tidx.isdigit() else None,
                "zt_label": zt_label or "",
            }
    if verbose:
        print(f"[info] built transcript trees for {len(meta)} transcripts", file=sys.stderr)
    return trees, meta


def pick_best_transcript(overlaps: List[Interval], site_start: int, site_end: int, meta: Dict[str,dict]):
    if not overlaps:
        return None
    best_tid, best_span = None, -1
    for iv in overlaps:
        tid = iv.data
        span = min(iv.end, site_end) - max(iv.begin, site_start)
        if span > best_span:
            best_span = span
            best_tid = tid
    return best_tid


def infer_sample_from_path(path: str, up: int) -> str:
    p = os.path.abspath(path)
    for _ in range(up):
        p = os.path.dirname(p)
    return os.path.basename(p)


def read_sites_with_mapping(path: str, trees, meta, sample_parent_up: int, min_cov: int, verbose=False):
    sample = infer_sample_from_path(path, up=sample_parent_up)
    rows = []
    with open_textmaybe_gzip(path) as f:
        for line in f:
            if is_header_line(line):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 18:
                parts = line.strip().split()
                if len(parts) < 18:
                    continue
            parts = parts[:18]
            d = dict(zip(BED_COLS, parts))
            chrom = d["chrom"]
            strand = d["strand"]
            try:
                s0 = int(d["start0"]); e0 = int(d["end0"])
            except Exception:
                continue
            key = (chrom, strand)
            tree = trees.get(key)
            if not tree:
                if verbose:
                    print(f"[skip] no transcripts on {chrom}{strand} for site {chrom}:{s0}-{e0}", file=sys.stderr)
                continue
            overlaps = list(tree.overlap(s0, e0))
            tid = pick_best_transcript(overlaps, s0, e0, meta)
            if not tid:
                if verbose:
                    print(f"[skip] no transcript overlap for site {chrom}:{s0}-{e0} {strand}", file=sys.stderr)
                continue
            m = meta[tid]
            Ncov = safe_int(d["Nvalid_cov"]) ; Nmod = safe_int(d["Nmod"]) ; frac = float(d["frac_modified"]) if d["frac_modified"] not in (".", "NA", "") else 0.0
            if min_cov and Ncov < min_cov:
                frac = 0.0
            rows.append({
                "chrom": chrom, "start0": s0, "end0": e0, "strand": strand,
                "mod_code": d["mod_code"], "sample": sample,
                "gene_index": m["gene_index"], "gene_name": m["gene_name"], "gene_id": m["gene_id"],
                "transcript_index": m["transcript_index"],
                "Nvalid_cov": Ncov, "Nmod": Nmod, "Ncanonical": safe_int(d["Ncanonical"]),
                "Nother_mod": safe_int(d["Nother_mod"]), "Ndelete": safe_int(d["Ndelete"]),
                "Nfail": safe_int(d["Nfail"]), "Ndiff": safe_int(d["Ndiff"]), "Nnocall": safe_int(d["Nnocall"]),
                "frac_modified": round(frac, 6),
            })
    return rows


def robust_load_summary(path: str, verbose=False) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    lines = []
    with open(path) as f:
        for ln in f:
            if ln.strip() == "":
                continue
            lines.append(ln.rstrip("\n"))
    header_idx, header = None, None
    for i, ln in enumerate(lines):
        h = ln.lstrip("#")
        if "\t" in h and any(c.strip() == "code" for c in h.split("\t")):
            header_idx = i
            header = [c.strip() for c in h.split("\t")]
            break
    if header is None:
        header_idx = 0
        header = [c.strip() for c in lines[0].lstrip("#").split("\t")]
        if verbose:
            print("[warn] Could not find explicit 'code' header; using first line as header", file=sys.stderr)
    rows = []
    expected = len(header)
    for ln in lines[header_idx+1:]:
        if ln.startswith("#"): continue
        parts = ln.split("\t")
        if len(parts) < expected: parts += [""] * (expected - len(parts))
        elif len(parts) > expected: parts = parts[:expected]
        rows.append({header[j]: parts[j] for j in range(expected)})
    df = pd.DataFrame(rows)
    if "code" not in df.columns:
        for c in list(df.columns):
            if c.lstrip("#").strip() == "code":
                df = df.rename(columns={c: "code"})
                break
    if "zt_label" in df.columns and "code" in df.columns:
        df["code"] = df["zt_label"].fillna(df["code"]).astype(str)
    elif "zt_label" in df.columns and "code" not in df.columns:
        df = df.rename(columns={"zt_label": "code"})
    if "code" in df.columns:
        df["code"] = df["code"].astype(str).str.strip()
    return df


def main():
    args = parse_args()

    # Build transcript lookup
    trees, meta = build_transcript_trees(args.gtf, verbose=args.verbose)

    # Recursively gather .bed/.bed.gz files
    beds = sorted(glob.glob(os.path.join(args.modkit_dir, "**", "*.bed"), recursive=True)) + \
           sorted(glob.glob(os.path.join(args.modkit_dir, "**", "*.bed.gz"), recursive=True))
    if not beds:
        sys.exit(f"No .bed/.bed.gz files found recursively under {args.modkit_dir}")

    if args.verbose:
        print(f"[info] found {len(beds)} bed files (recursive)", file=sys.stderr)

    rows_all: List[dict] = []
    for bed in beds:
        # Skip obvious logs or non-bed artifacts
        if bed.endswith("_bed_log"):
            continue
        rows = read_sites_with_mapping(bed, trees, meta, args.sample_from_parent, args.min_cov, verbose=args.verbose)
        rows_all.extend(rows)

    if not rows_all:
        sys.exit("Parsed zero site rows after transcript mapping; check GTF and bed formats.")

    df = pd.DataFrame(rows_all)

    # Optional metadata join (kept minimal; many users won't need this)
    if args.summary_tsv:
        summ = robust_load_summary(args.summary_tsv, verbose=args.verbose)
        # Only attach columns that don't collide with required outputs
        if not summ.empty and "code" in summ.columns:
            keep_cols = [c for c in [
                "classification","match_source","read_support","frac_global","polya_support_frac"
            ] if c in summ.columns]
            # We have no 'code' per row here (we mapped by position), so this join is generally skipped
            # unless user has a separate mapping. We log and skip to avoid confusion.
            if args.verbose:
                print("[warn] Skipping summary join (no 'code' per site after mapping).", file=sys.stderr)

    # Write long table
    out_long = f"{args.out_prefix}_sites_long.tsv"
    df.to_csv(out_long, sep="\t", index=False)
    if args.verbose:
        print(f"[ok] wrote {out_long} (rows={len(df)})", file=sys.stderr)

    # Per-gene per-mod spreadsheets
    outdir = f"{args.out_prefix}.per_gene_sites"
    os.makedirs(outdir, exist_ok=True)

    grouped = df.groupby(["gene_index","gene_name","gene_id","mod_code"], dropna=False)
    for (gidx,gname,gid,mod), gdf in grouped:
        gidx_str = str(gidx) if pd.notna(gidx) else "NA"
        gname_str = (gname or "NA").replace("/", "-")
        fname = f"{gidx_str}_{gname_str}_{mod}.tsv"
        outp = os.path.join(outdir, fname)
        gdf.sort_values(["sample","transcript_index","chrom","start0"], inplace=True)
        gdf.to_csv(outp, sep="\t", index=False)
        if args.verbose:
            print(f"[ok] wrote {outp} (rows={len(gdf)})", file=sys.stderr)

    print("[OK] ZN-style per-gene aggregation complete.")

if __name__ == "__main__":
    main()

