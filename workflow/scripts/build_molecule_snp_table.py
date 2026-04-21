#!/usr/bin/env python3

import argparse
import os

import pandas as pd
import pysam

from genotype_utils import sample_name_from_bam, safe_int


def parse_args():
    ap = argparse.ArgumentParser(description="Build a per-read candidate SNP table from tagged BAMs.")
    ap.add_argument("--bams", nargs="+", required=True, help="Input BAMs")
    ap.add_argument("--candidate-snps", required=True, help="Candidate SNP TSV from discover_candidate_snps.py")
    ap.add_argument("--out-tsv", required=True, help="Output molecule SNP TSV")
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--primary-only", action="store_true")
    return ap.parse_args()


def safe_get_tag(aln, tag, default=""):
    try:
        return aln.get_tag(tag)
    except Exception:
        return default


def build_windows(positions, max_gap=5000):
    if not positions:
        return []
    positions = sorted(set(int(p) for p in positions))
    windows = []
    start = positions[0]
    prev = positions[0]
    for pos in positions[1:]:
        if pos - prev > max_gap:
            windows.append((start, prev))
            start = pos
        prev = pos
    windows.append((start, prev))
    return windows


def main():
    args = parse_args()
    cand = pd.read_csv(args.candidate_snps, sep="\t", low_memory=False)
    if cand.empty:
        out = pd.DataFrame(columns=[
            "sample", "qname", "snp_id", "chrom", "pos1", "start0", "end0", "ref", "alt",
            "observed_base", "allele_class", "baseq", "mapq", "strand",
            "ZT", "ZG", "ZN", "ZM", "gene_names", "gene_ids", "metagene_indices"
        ])
        os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
        out.to_csv(args.out_tsv, sep="\t", index=False)
        return

    cand_by_chrom = {}
    for chrom, sub in cand.groupby("chrom", sort=False):
        pos_map = {}
        for row in sub.to_dict("records"):
            pos_map[int(row["pos1"])] = row
        cand_by_chrom[chrom] = pos_map

    rows = []
    for bam in args.bams:
        sample = sample_name_from_bam(bam)
        with pysam.AlignmentFile(bam, "rb") as fh:
            for chrom, pos_map in cand_by_chrom.items():
                windows = build_windows(pos_map.keys())
                for win_start1, win_end1 in windows:
                    wanted = pos_map
                    for col in fh.pileup(
                        chrom,
                        win_start1 - 1,
                        win_end1,
                        truncate=True,
                        stepper="samtools",
                        min_base_quality=args.min_baseq,
                        min_mapping_quality=args.min_mapq,
                    ):
                        pos1 = int(col.reference_pos) + 1
                        cand_row = wanted.get(pos1)
                        if cand_row is None:
                            continue
                        ref = str(cand_row["ref"])
                        alt = str(cand_row["alt"])
                        for pr in col.pileups:
                            aln = pr.alignment
                            if pr.is_del or pr.is_refskip or pr.query_position is None:
                                continue
                            if args.primary_only and (aln.is_secondary or aln.is_supplementary):
                                continue
                            if aln.mapping_quality < args.min_mapq:
                                continue
                            qpos = pr.query_position
                            bq = int(aln.query_qualities[qpos]) if aln.query_qualities is not None else 0
                            if bq < args.min_baseq:
                                continue
                            base = aln.query_sequence[qpos].upper()
                            if len(base) != 1:
                                continue
                            allele_class = "other"
                            if base == ref:
                                allele_class = "ref"
                            elif base == alt:
                                allele_class = "alt"
                            rows.append({
                                "sample": sample,
                                "qname": aln.query_name,
                                "snp_id": cand_row["snp_id"],
                                "chrom": chrom,
                                "pos1": pos1,
                                "start0": pos1 - 1,
                                "end0": pos1,
                                "ref": ref,
                                "alt": alt,
                                "observed_base": base,
                                "allele_class": allele_class,
                                "baseq": bq,
                                "mapq": int(aln.mapping_quality),
                                "strand": "-" if aln.is_reverse else "+",
                                "ZT": str(safe_get_tag(aln, "ZT", "")),
                                "ZG": safe_int(safe_get_tag(aln, "ZG", "")),
                                "ZN": safe_int(safe_get_tag(aln, "ZN", "")),
                                "ZM": safe_int(safe_get_tag(aln, "ZM", "")),
                                "gene_names": cand_row.get("gene_names", ""),
                                "gene_ids": cand_row.get("gene_ids", ""),
                                "metagene_indices": cand_row.get("metagene_indices", ""),
                            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates(["sample", "qname", "snp_id"], keep="first").reset_index(drop=True)
    else:
        df = pd.DataFrame(columns=[
            "sample", "qname", "snp_id", "chrom", "pos1", "start0", "end0", "ref", "alt",
            "observed_base", "allele_class", "baseq", "mapq", "strand",
            "ZT", "ZG", "ZN", "ZM", "gene_names", "gene_ids", "metagene_indices"
        ])

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    df.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
