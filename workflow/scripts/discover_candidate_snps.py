#!/usr/bin/env python3

import argparse
from collections import defaultdict
import os
import re

import pandas as pd
import pysam

from genotype_utils import sample_name_from_bam


DNA_BASES = ("A", "C", "G", "T")


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


def main():
    args = parse_args()
    exon_records, merged_intervals = load_gtf_exons(args.gtf)
    fasta = pysam.FastaFile(args.reference_fa)

    site_counts = {}
    for bam in args.bams:
        sample = sample_name_from_bam(bam)
        if args.verbose:
            print(f"[info] scanning {sample}", file=os.sys.stderr)
        with pysam.AlignmentFile(bam, "rb") as fh:
            for chrom, intervals in merged_intervals.items():
                for start1, end1 in intervals:
                    ref_seq = fasta.fetch(chrom, start1 - 1, end1).upper()
                    for col in fh.pileup(
                        chrom,
                        start1 - 1,
                        end1,
                        truncate=True,
                        stepper="samtools",
                        min_base_quality=args.min_baseq,
                        min_mapping_quality=args.min_mapq,
                    ):
                        pos1 = int(col.reference_pos) + 1
                        ref_idx = pos1 - start1
                        if ref_idx < 0 or ref_idx >= len(ref_seq):
                            continue
                        ref = ref_seq[ref_idx]
                        if ref not in DNA_BASES:
                            continue
                        rec = site_counts.setdefault((chrom, pos1, ref), init_count_record())
                        per_sample = rec["per_sample"][sample]
                        for pr in col.pileups:
                            aln = pr.alignment
                            if pr.is_del or pr.is_refskip or pr.query_position is None:
                                continue
                            if args.primary_only and (aln.is_secondary or aln.is_supplementary):
                                continue
                            if aln.mapping_quality < args.min_mapq:
                                continue
                            qpos = pr.query_position
                            if aln.query_qualities is not None and aln.query_qualities[qpos] < args.min_baseq:
                                continue
                            base = aln.query_sequence[qpos].upper()
                            if base not in DNA_BASES:
                                continue
                            rec["counts"][base] += 1
                            per_sample[base] += 1

    rows = []
    for (chrom, pos1, ref), rec in sorted(site_counts.items()):
        counts = rec["counts"]
        total_cov = sum(counts.values())
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
        for sample in sorted(rec["per_sample"]):
            sc = rec["per_sample"][sample]
            if sc.get(alt, 0) > 0:
                samples_with_alt += 1
            sample_summaries.append(
                f"{sample}:A={sc['A']},C={sc['C']},G={sc['G']},T={sc['T']}"
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
