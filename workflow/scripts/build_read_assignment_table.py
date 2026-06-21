#!/usr/bin/env python3

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


def parse_args():
    ap = argparse.ArgumentParser(description="Build a read-level assignment table from ZT/ZN-tagged BAMs.")
    ap.add_argument("--bams", nargs="+", required=True, help="Tagged or cleaned BAMs containing ZT/ZG/ZN/ZM tags")
    ap.add_argument("--summary-tsv", default="", help="Classification summary TSV to join transcript metadata")
    ap.add_argument("--out-tsv", required=True, help="Output TSV")
    ap.add_argument("--jobs", type=int, default=1, help="Number of BAMs to scan in parallel")
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
        print(f"[info] read assignments start: {sample} {region or 'all'}", file=sys.stderr, flush=True)

    rows = []
    with pysam.AlignmentFile(bam, "rb") as fh:
        for aln in (fh.fetch(contig=region) if region else fh.fetch()):
            if aln.is_unmapped:
                continue
            if primary_only and (aln.is_secondary or aln.is_supplementary):
                continue
            rows.append({
                "sample": sample,
                "qname": aln.query_name,
                "chrom": fh.get_reference_name(aln.reference_id),
                "start0": int(aln.reference_start),
                "end0": int(aln.reference_end or aln.reference_start),
                "strand": "-" if aln.is_reverse else "+",
                "mapq": int(aln.mapping_quality),
                "ZT": str(safe_get_tag(aln, "ZT", "")),
                "ZG": safe_int(safe_get_tag(aln, "ZG", "")),
                "ZN": safe_int(safe_get_tag(aln, "ZN", "")),
                "ZM": safe_int(safe_get_tag(aln, "ZM", "")),
            })

    if verbose:
        print(f"[info] read assignments done: {sample} rows={len(rows)}", file=sys.stderr, flush=True)
    return rows


def main():
    args = parse_args()
    # Shard per (BAM x chromosome-with-reads) for parallelism beyond the sample count.
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
        for result in run_process_jobs(
            collect_rows_from_bam,
            task_args,
            jobs,
            verbose=args.verbose,
            label="build_read_assignment_table",
        ):
            rows.extend(result)

    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(columns=["sample", "qname", "chrom", "start0", "end0", "strand", "mapq", "ZT", "ZG", "ZN", "ZM"])
    df["assigned"] = normalize_string_series(df.get("ZT", pd.Series(dtype=str))).ne("")

    summ = robust_load_summary(args.summary_tsv) if args.summary_tsv else pd.DataFrame()
    if not summ.empty and "zt_label" in summ.columns:
        keep = [
            c for c in [
                "zt_label", "gtf_gene_id", "gtf_gene_name", "gene_index", "transcript_index",
                "metagene_index", "zn_index", "metagene_partition_count", "classification",
                "match_source", "assignment_mode", "read_support"
            ] if c in summ.columns
        ]
        meta = summ[keep].drop_duplicates("zt_label").rename(columns={
            "zt_label": "ZT",
            "gtf_gene_id": "gene_id",
            "gtf_gene_name": "gene_name",
            "zn_index": "summary_zn_index",
        })
        df = df.merge(meta, on="ZT", how="left")

    if not df.empty:
        # Deterministic on-disk order across (BAM x chrom-with-reads) shards; also makes
        # the downstream drop_duplicates(["sample","qname"]) selection reproducible.
        sort_cols = [c for c in ["sample", "qname", "chrom", "start0", "end0", "strand"] if c in df.columns]
        if sort_cols:
            df = df.sort_values(sort_cols).reset_index(drop=True)

    out_dir = os.path.dirname(args.out_tsv) or "."
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
