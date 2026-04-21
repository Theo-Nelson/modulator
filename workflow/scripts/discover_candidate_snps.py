#!/usr/bin/env python3

import argparse
from collections import defaultdict
import os
import re
import sys

import pandas as pd
import pysam

from genotype_utils import run_process_jobs, sample_name_from_bam


DNA_BASES = ("A", "C", "G", "T")
BASE_INDEX = {base: idx for idx, base in enumerate(DNA_BASES)}


def parse_args():
    ap = argparse.ArgumentParser(description="Discover segregating candidate SNPs from tagged BAMs inside assembled transcript loci.")
    ap.add_argument("--bams", nargs="+", required=True, help="Input BAMs")
    ap.add_argument("--reference-fa", required=True, help="Reference FASTA")
    ap.add_argument("--gtf", required=True, help="Assembler GTF used to define transcribed loci")
    ap.add_argument("--out-tsv", required=True, help="Output candidate SNP TSV")
    ap.add_argument("--min-alt-reads", type=int, default=4)
    ap.add_argument("--min-total-cov", type=int, default=8)
    ap.add_argument("--min-alt-frac", type=float, default=0.10)
    ap.add_argument("--max-alt-frac", type=float, default=0.90)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=1, help="Number of BAMs to scan in parallel")
    ap.add_argument("--primary-only", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def parse_attrs(attrs: str):
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attrs):
        out[m.group(1)] = m.group(2)
    return out


def load_gtf_exons(gtf_path: str):
    exon_records = defaultdict(list)
    merged_intervals = defaultdict(list)
    with open(gtf_path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = parts
            if feature != "exon":
                continue
            a = parse_attrs(attrs)
            start = int(start)
            end = int(end)
            exon_records[chrom].append({
                "start": start,
                "end": end,
                "strand": strand,
                "gene_id": a.get("gene_id", ""),
                "gene_name": a.get("ref_gene_name", a.get("gene_name", a.get("gene_id", ""))),
                "metagene_index": a.get("metagene_index", ""),
                "zt_label": a.get("zt_label", a.get("transcript_id", "")),
            })
            merged_intervals[chrom].append((start, end))

    for chrom, intervals in merged_intervals.items():
        intervals = sorted(intervals)
        merged = []
        for s, e in intervals:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        merged_intervals[chrom] = [(s, e) for s, e in merged]
        exon_records[chrom].sort(key=lambda x: (x["start"], x["end"], x["gene_name"], x["zt_label"]))
    return exon_records, merged_intervals


def annotate_site(chrom: str, pos1: int, exon_records):
    genes = set()
    gene_ids = set()
    metagenes = set()
    zts = set()
    for rec in exon_records.get(chrom, []):
        if rec["start"] > pos1:
            break
        if rec["end"] < pos1:
            continue
        gene_ids.add(rec["gene_id"])
        genes.add(rec["gene_name"])
        if rec["metagene_index"]:
            metagenes.add(str(rec["metagene_index"]))
        if rec["zt_label"]:
            zts.add(rec["zt_label"])
    return {
        "gene_ids": ";".join(sorted(g for g in gene_ids if g)),
        "gene_names": ";".join(sorted(g for g in genes if g)),
        "metagene_indices": ";".join(sorted(m for m in metagenes if m)),
        "zt_labels": ";".join(sorted(z for z in zts if z)),
    }


def init_count_record():
    return {
        "counts": {b: 0 for b in DNA_BASES},
        "per_sample": defaultdict(lambda: {b: 0 for b in DNA_BASES}),
    }


def make_read_callback(primary_only: bool, min_mapq: int):
    if (not primary_only) and int(min_mapq) <= 0:
        return "all"

    def callback(read):
        if read.is_unmapped:
            return False
        if primary_only and (read.is_secondary or read.is_supplementary):
            return False
        if read.mapping_quality < int(min_mapq):
            return False
        return True

    return callback


def scan_bam_counts(
    bam: str,
    reference_fa: str,
    merged_intervals,
    min_baseq: int,
    min_mapq: int,
    primary_only: bool,
    verbose: bool = False,
):
    sample = sample_name_from_bam(bam)
    if verbose:
        print(f"[info] SNP scan start: {sample}", file=sys.stderr, flush=True)

    sample_counts = {}
    fasta = pysam.FastaFile(reference_fa)
    callback = make_read_callback(primary_only, min_mapq)
    covered_positions = 0
    try:
        with pysam.AlignmentFile(bam, "rb") as fh:
            for chrom, intervals in merged_intervals.items():
                for start1, end1 in intervals:
                    ref_seq = fasta.fetch(chrom, start1 - 1, end1).upper()
                    counts = fh.count_coverage(
                        chrom,
                        start1 - 1,
                        end1,
                        quality_threshold=int(min_baseq),
                        read_callback=callback,
                    )
                    for offset, ref in enumerate(ref_seq):
                        if ref not in BASE_INDEX:
                            continue
                        base_counts = tuple(int(arr[offset]) for arr in counts)
                        if not any(base_counts):
                            continue
                        covered_positions += 1
                        sample_counts[(chrom, start1 + offset, ref)] = base_counts
    finally:
        fasta.close()

    if verbose:
        print(
            f"[info] SNP scan done: {sample} positions_with_coverage={covered_positions}",
            file=sys.stderr,
            flush=True,
        )
    return sample, sample_counts


def main():
    args = parse_args()
    exon_records, merged_intervals = load_gtf_exons(args.gtf)
    jobs = max(1, min(int(args.jobs), len(args.bams)))

    results = []
    worker_args = [
        (bam, args.reference_fa, merged_intervals, args.min_baseq, args.min_mapq, args.primary_only, args.verbose)
        for bam in args.bams
    ]
    if jobs == 1:
        for item in worker_args:
            results.append(scan_bam_counts(*item))
    else:
        results = run_process_jobs(
            scan_bam_counts,
            worker_args,
            jobs,
            verbose=args.verbose,
            label="discover_candidate_snps",
        )

    site_counts = {}
    for sample, sample_counts in results:
        for key, sample_base_counts in sample_counts.items():
            rec = site_counts.setdefault(key, {"counts": [0, 0, 0, 0], "per_sample": {}})
            rec["per_sample"][sample] = sample_base_counts
            for idx, value in enumerate(sample_base_counts):
                rec["counts"][idx] += int(value)

    rows = []
    for (chrom, pos1, ref), rec in sorted(site_counts.items()):
        total_counts = rec["counts"]
        counts = {base: int(total_counts[idx]) for idx, base in enumerate(DNA_BASES)}
        total_cov = sum(total_counts)
        ref_count = counts.get(ref, 0)
        alts = sorted(((b, c) for b, c in counts.items() if b != ref), key=lambda x: (-x[1], x[0]))
        alt, alt_count = alts[0] if alts else ("", 0)
        second_alt_count = alts[1][1] if len(alts) > 1 else 0
        alt_frac = (alt_count / total_cov) if total_cov > 0 else 0.0
        ref_frac = (ref_count / total_cov) if total_cov > 0 else 0.0
        if total_cov < args.min_total_cov:
            continue
        if alt_count < args.min_alt_reads:
            continue
        if alt_frac < args.min_alt_frac:
            continue
        if alt_frac > args.max_alt_frac:
            continue
        if second_alt_count >= args.min_alt_reads:
            continue

        ann = annotate_site(chrom, pos1, exon_records)
        sample_summaries = []
        samples_with_alt = 0
        alt_idx = BASE_INDEX.get(alt, -1)
        for sample in sorted(rec["per_sample"]):
            sc = rec["per_sample"][sample]
            if alt_idx >= 0 and sc[alt_idx] > 0:
                samples_with_alt += 1
            sample_summaries.append(
                f"{sample}:A={sc[0]},C={sc[1]},G={sc[2]},T={sc[3]}"
            )
        rows.append({
            "snp_id": f"{chrom}:{pos1}:{ref}>{alt}",
            "chrom": chrom,
            "pos1": pos1,
            "start0": pos1 - 1,
            "end0": pos1,
            "ref": ref,
            "alt": alt,
            "total_cov": total_cov,
            "ref_count": ref_count,
            "ref_frac": round(ref_frac, 6),
            "alt_count": alt_count,
            "second_alt_count": second_alt_count,
            "alt_frac": round(alt_frac, 6),
            "site_class": "segregating",
            "samples_with_alt": samples_with_alt,
            "gene_ids": ann["gene_ids"],
            "gene_names": ann["gene_names"],
            "metagene_indices": ann["metagene_indices"],
            "zt_labels": ann["zt_labels"],
            "sample_base_counts": "|".join(sample_summaries),
        })

    df = pd.DataFrame(rows)
    out_dir = os.path.dirname(args.out_tsv) or "."
    os.makedirs(out_dir, exist_ok=True)
    if df.empty:
        df = pd.DataFrame(columns=[
            "snp_id", "chrom", "pos1", "start0", "end0", "ref", "alt", "total_cov", "ref_count",
            "ref_frac", "alt_count", "second_alt_count", "alt_frac", "site_class", "samples_with_alt", "gene_ids",
            "gene_names", "metagene_indices", "zt_labels", "sample_base_counts"
        ])
    df.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
