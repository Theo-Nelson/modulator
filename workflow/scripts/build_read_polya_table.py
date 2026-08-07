#!/usr/bin/env python3
"""Per-read poly(A) tail-length table from ZT-tagged BAMs.

Reads the dorado poly(A) estimate (``pt:i`` tag, tail length in bases) for every ASSIGNED read
(one carrying a fragmentform ``ZT`` tag) in one pass over the ZT-tagged BAMs -- the same BAMs the
assembler/read_stats already produce, which retain ``pt:i`` alongside ``ZT/ZG/ZN/ZM``. The tag is
written by dorado basecalling with ``--estimate-poly-a``; reads with no estimate carry ``pt:i:0``.

Output is one row per (sample, qname) assigned read with its tail length pre-joined to the
fragmentform, so downstream steps group tail lengths by fragmentform / gene / metagene / sample and
join to ``molecule_mod_calls`` (same ``(sample, qname)`` key) for tail x modification. Mirrors
build_read_assignment_table.py (per-BAM x chrom sharding, run_process_jobs, atomic write).
"""
import argparse
import os
import sys

import pandas as pd
import pysam

from genotype_utils import (
    normalize_string_series,
    robust_load_summary,
    run_process_jobs,
    sample_name_from_bam,
    safe_int,
)

OUT_COLS = ["sample", "qname", "strand", "tail_len", "tail_estimated",
            "ZT", "ZG", "ZN", "ZM", "gene_name", "metagene_index", "transcript_index", "classification"]


def parse_args():
    ap = argparse.ArgumentParser(description="Per-read poly(A) tail length (dorado pt:i) from ZT-tagged BAMs.")
    ap.add_argument("--bams", nargs="+", required=True, help="ZT-tagged BAMs carrying pt:i + ZT/ZG/ZN/ZM tags")
    ap.add_argument("--summary-tsv", default="", help="Classification summary TSV to join gene/metagene/classification")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--jobs", type=int, default=1, help="Number of BAM x chrom shards to scan in parallel")
    ap.add_argument("--primary-only", action="store_true", help="Skip secondary/supplementary alignments")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def safe_get_tag(aln, tag, default=""):
    try:
        return aln.get_tag(tag)
    except Exception:
        return default


def bam_chroms_with_reads(bam: str):
    try:
        with pysam.AlignmentFile(bam, "rb") as fh:
            return [s.contig for s in fh.get_index_statistics() if s.mapped > 0]
    except Exception:
        return []


def collect_rows_from_bam(bam: str, primary_only: bool, verbose: bool = False, region=None):
    sample = sample_name_from_bam(bam)
    if verbose:
        print(f"[info] polya scan start: {sample} {region or 'all'}", file=sys.stderr, flush=True)
    rows = []
    with pysam.AlignmentFile(bam, "rb") as fh:
        for aln in (fh.fetch(contig=region) if region else fh.fetch()):
            if aln.is_unmapped:
                continue
            if primary_only and (aln.is_secondary or aln.is_supplementary):
                continue
            zt = str(safe_get_tag(aln, "ZT", ""))
            if not zt:
                continue  # assigned reads only -- unassigned reads have no fragmentform
            tail = safe_int(safe_get_tag(aln, "pt", 0))
            rows.append({
                "sample": sample,
                "qname": aln.query_name,
                "strand": "-" if aln.is_reverse else "+",
                "tail_len": int(tail),
                "tail_estimated": bool(tail > 0),  # dorado writes pt:i:0 when it finds no tail
                "ZT": zt,
                "ZG": safe_int(safe_get_tag(aln, "ZG", "")),
                "ZN": safe_int(safe_get_tag(aln, "ZN", "")),
                "ZM": safe_int(safe_get_tag(aln, "ZM", "")),
            })
    if verbose:
        print(f"[info] polya scan done: {sample} assigned_rows={len(rows)}", file=sys.stderr, flush=True)
    return rows


def main():
    args = parse_args()
    task_args = []
    for bam in args.bams:
        chroms = bam_chroms_with_reads(bam)
        if chroms:
            task_args.extend((bam, args.primary_only, args.verbose, c) for c in chroms)
        else:
            task_args.append((bam, args.primary_only, args.verbose, None))
    jobs = max(1, min(int(args.jobs), len(task_args)))
    rows = []
    if jobs == 1:
        for item in task_args:
            rows.extend(collect_rows_from_bam(*item))
    else:
        for result in run_process_jobs(collect_rows_from_bam, task_args, jobs,
                                       verbose=args.verbose, label="build_read_polya_table"):
            rows.extend(result)

    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(columns=["sample", "qname", "strand", "tail_len", "tail_estimated", "ZT", "ZG", "ZN", "ZM"])

    # Join gene_name / metagene_index / transcript_index / fragmentform classification from the summary.
    summ = robust_load_summary(args.summary_tsv) if args.summary_tsv else pd.DataFrame()
    if not summ.empty and "zt_label" in summ.columns:
        keep = [c for c in ["zt_label", "gtf_gene_name", "gene_index", "transcript_index",
                            "metagene_index", "classification"] if c in summ.columns]
        meta = summ[keep].drop_duplicates("zt_label").rename(columns={
            "zt_label": "ZT", "gtf_gene_name": "gene_name"})
        df = df.merge(meta, on="ZT", how="left")
    # Fallbacks derivable from the ZT tag itself if no summary provided.
    if "gene_name" not in df.columns:
        df["gene_name"] = normalize_string_series(df.get("ZT", pd.Series(dtype=str))).str.split(".").str[0]
    if "metagene_index" not in df.columns:
        df["metagene_index"] = df.get("ZM", pd.Series(dtype="Int64"))
    for c in OUT_COLS:
        if c not in df.columns:
            df[c] = pd.NA
    df = df[OUT_COLS]

    if not df.empty:
        df = df.sort_values([c for c in ["sample", "qname"] if c in df.columns]).reset_index(drop=True)

    out_dir = os.path.dirname(args.out_tsv) or "."
    os.makedirs(out_dir, exist_ok=True)
    _tmp = args.out_tsv + ".tmp"
    df.to_csv(_tmp, sep="\t", index=False)
    os.replace(_tmp, args.out_tsv)
    if args.verbose:
        print(f"[info] wrote {len(df)} assigned-read tail rows -> {args.out_tsv}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
