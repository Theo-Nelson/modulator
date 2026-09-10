#!/usr/bin/env python3

import argparse
import os
import shutil
import sys
import tempfile

import pandas as pd
import pysam

from genotype_utils import (
    normalize_string_series,
    robust_load_summary,
    run_process_jobs,
    sample_name_from_bam,
    safe_int,
)

BASE_COLS = ["sample", "qname", "chrom", "start0", "end0", "strand", "mapq", "ZT", "ZG", "ZN", "ZM"]


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


def collect_rows_from_bam(bam: str, primary_only: bool, verbose: bool = False, region=None, shard_dir=None):
    """Scan one (BAM x chromosome) and write its rows to a pickle shard. Returns
    (sample, shard_path_or_None, n_rows, n_unassigned).

    The rows are NOT returned to the parent: at cohort scale every read of every sample as a Python
    dict in one process is >100 GiB. The parent assembles the table one SAMPLE at a time from the
    shards (the output is sorted by sample first, so per-sample blocks reproduce the global order)."""
    sample = sample_name_from_bam(bam)
    if verbose:
        print(f"[info] read assignments start: {sample} {region or 'all'}", file=sys.stderr, flush=True)

    rows = []
    n_unassigned = 0
    with pysam.AlignmentFile(bam, "rb") as fh:
        for aln in (fh.fetch(contig=region) if region else fh.fetch()):
            if aln.is_unmapped:
                continue
            if primary_only and (aln.is_secondary or aln.is_supplementary):
                continue
            zt = str(safe_get_tag(aln, "ZT", ""))
            if zt == "":
                n_unassigned += 1
            rows.append({
                "sample": sample,
                "qname": aln.query_name,
                "chrom": fh.get_reference_name(aln.reference_id),
                "start0": int(aln.reference_start),
                "end0": int(aln.reference_end or aln.reference_start),
                "strand": "-" if aln.is_reverse else "+",
                "mapq": int(aln.mapping_quality),
                "ZT": zt,
                "ZG": safe_int(safe_get_tag(aln, "ZG", "")),
                "ZN": safe_int(safe_get_tag(aln, "ZN", "")),
                "ZM": safe_int(safe_get_tag(aln, "ZM", "")),
            })

    if verbose:
        print(f"[info] read assignments done: {sample} rows={len(rows)}", file=sys.stderr, flush=True)
    if not rows:
        return sample, None, 0, 0
    path = os.path.join(shard_dir, f"{sample}.{region or 'all'}.pkl")
    pd.DataFrame(rows)[BASE_COLS].to_pickle(path)
    return sample, path, len(rows), n_unassigned


def _meta_frame(summary_tsv: str):
    summ = robust_load_summary(summary_tsv) if summary_tsv else pd.DataFrame()
    if summ.empty or "zt_label" not in summ.columns:
        return None
    keep = [
        c for c in [
            "zt_label", "gtf_gene_id", "gtf_gene_name", "gene_index", "transcript_index",
            "metagene_index", "zn_index", "metagene_partition_count", "classification",
            "match_source", "assignment_mode", "read_support"
        ] if c in summ.columns
    ]
    return summ[keep].drop_duplicates("zt_label").rename(columns={
        "zt_label": "ZT",
        "gtf_gene_id": "gene_id",
        "gtf_gene_name": "gene_name",
        "zn_index": "summary_zn_index",
    })


def _finish_frame(df: pd.DataFrame, meta, cast_float_cols):
    df["assigned"] = normalize_string_series(df.get("ZT", pd.Series(dtype=str))).ne("")
    if meta is not None:
        df = df.merge(meta, on="ZT", how="left")
        # A left merge that leaves ANY row unmatched turns the summary's integer columns into float64
        # (NaN) for the WHOLE table. The table is now written per sample, so a sample whose reads are
        # all assigned must still write those columns as floats to keep one dtype across the file.
        for c in cast_float_cols:
            if c in df.columns:
                df[c] = df[c].astype("float64")
    if not df.empty:
        sort_cols = [c for c in ["sample", "qname", "chrom", "start0", "end0", "strand"] if c in df.columns]
        if sort_cols:
            df = df.sort_values(sort_cols).reset_index(drop=True)
    return df


def main():
    args = parse_args()
    shard_dir = tempfile.mkdtemp(prefix=".read_assign_shards.", dir=os.path.dirname(args.out_tsv) or ".")
    try:
        # Shard per (BAM x chromosome-with-reads) for parallelism beyond the sample count.
        task_args = []
        for bam in args.bams:
            chroms = bam_chroms_with_reads(bam)
            if chroms:
                task_args.extend((bam, args.primary_only, args.verbose, c, shard_dir) for c in chroms)
            else:
                task_args.append((bam, args.primary_only, args.verbose, None, shard_dir))
        jobs = max(1, min(int(args.jobs), len(task_args)))
        if jobs == 1:
            results = [collect_rows_from_bam(*item) for item in task_args]
        else:
            results = run_process_jobs(
                collect_rows_from_bam, task_args, jobs, verbose=args.verbose, label="build_read_assignment_table")

        by_sample = {}
        any_unassigned = False
        for sample, path, n, n_un in results:
            if path is not None:
                by_sample.setdefault(sample, []).append(path)
            if n_un:
                any_unassigned = True

        meta = _meta_frame(args.summary_tsv)
        cast_float_cols = []
        if meta is not None and any_unassigned:
            cast_float_cols = [c for c in meta.columns if c != "ZT" and pd.api.types.is_integer_dtype(meta[c])]

        out_dir = os.path.dirname(args.out_tsv) or "."
        os.makedirs(out_dir, exist_ok=True)
        _tmp = args.out_tsv + ".tmp"          # atomic write: a partial file from an OOM/kill
        with open(_tmp, "w") as out:          # is never left named as the final output, so
            wrote_header = False              # --resume can safely reuse a non-empty final file.
            if not by_sample:
                df = _finish_frame(pd.DataFrame(columns=BASE_COLS), meta, cast_float_cols)
                df.to_csv(out, sep="\t", index=False)
                wrote_header = True
            for sample in sorted(by_sample):
                df = pd.concat([pd.read_pickle(p) for p in by_sample[sample]], ignore_index=True)
                df = _finish_frame(df, meta, cast_float_cols)
                df.to_csv(out, sep="\t", index=False, header=not wrote_header)
                wrote_header = True
                del df
        os.replace(_tmp, args.out_tsv)
    finally:
        shutil.rmtree(shard_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
