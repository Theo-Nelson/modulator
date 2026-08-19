#!/usr/bin/env python3
"""
build_candidate_regions_bed.py -- union of candidate SNP positions and candidate mod sites,
merged into a sorted, non-overlapping BED.

This BED is what the genotype stage uses to pre-subset each BAM (`samtools view -M -L`) before the
per-molecule scans. Every read that can contribute to snp/mod/haplotype association overlaps at
least one candidate site, so restricting the per-read tables to reads over these regions is lossless
for every genotype output -- while removing the genome-wide off-target reads that otherwise make
build_read_assignment_table materialize an all-reads table (~22 GB on disk / ~70 GB in pandas).

Emits an empty file (0 lines) when there are no candidates; the caller must then skip subsetting.
"""

import argparse
import os
import sys


def read_snp_intervals(path):
    """candidate_snps.tsv -> [(chrom, start0, end0)]"""
    out = []
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header:
            header[0] = header[0].lstrip("#")   # tolerate a '#'-commented header (else all SNPs dropped)
        idx = {c: i for i, c in enumerate(header)}
        if "chrom" not in idx:
            return out
        has_se = "start0" in idx and "end0" in idx
        has_pos1 = "pos1" in idx
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            chrom = p[idx["chrom"]]
            try:
                if has_se:
                    s, e = int(p[idx["start0"]]), int(p[idx["end0"]])
                elif has_pos1:
                    pos1 = int(p[idx["pos1"]])
                    s, e = pos1 - 1, pos1
                else:
                    continue
            except ValueError:
                continue
            if e > s:
                out.append((chrom, s, e))
    return out


def read_bed_intervals(path):
    out = []
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            if e > s:
                out.append((p[0], s, e))
    return out


def merge_intervals(intervals, pad, chrom_lengths=None):
    """Sort + merge overlapping/adjacent intervals (after padding). Pure Python -- no bedtools."""
    padded = []
    for chrom, s, e in intervals:
        s = max(0, s - pad)
        e = e + pad
        if chrom_lengths and chrom in chrom_lengths:
            e = min(e, chrom_lengths[chrom])
        if e > s:
            padded.append((chrom, s, e))
    padded.sort()
    merged = []
    for chrom, s, e in padded:
        if merged and merged[-1][0] == chrom and s <= merged[-1][2]:
            if e > merged[-1][2]:
                merged[-1][2] = e
        else:
            merged.append([chrom, s, e])
    return merged


def load_fai(path):
    lengths = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) >= 2:
                    try:
                        lengths[p[0]] = int(p[1])
                    except ValueError:
                        pass
    return lengths


def parse_args():
    ap = argparse.ArgumentParser(description="Merge candidate SNP + mod-site positions into one BED.")
    ap.add_argument("--candidate-snps", default="", help="candidate_snps.tsv")
    ap.add_argument("--candidate-mod-bed", default="", help="candidate_mod_sites.bed")
    ap.add_argument("--fai", default="", help="Reference .fai, to clamp padded intervals to contig ends")
    ap.add_argument("--pad", type=int, default=0,
                    help="Bases to pad each candidate site by before merging. Reads only need to OVERLAP "
                         "a site, so 0 is lossless; a small pad just widens the kept region.")
    ap.add_argument("--out-bed", required=True)
    return ap.parse_args()


def main():
    args = parse_args()
    ivals = read_snp_intervals(args.candidate_snps) + read_bed_intervals(args.candidate_mod_bed)
    merged = merge_intervals(ivals, max(0, int(args.pad)), load_fai(args.fai))

    os.makedirs(os.path.dirname(args.out_bed) or ".", exist_ok=True)
    with open(args.out_bed, "w") as out:
        for chrom, s, e in merged:
            out.write(f"{chrom}\t{s}\t{e}\n")

    total_bp = sum(e - s for _, s, e in merged)
    print(f"[ok] wrote {args.out_bed}: {len(merged)} merged region(s), {total_bp:,} bp "
          f"(from {len(ivals)} candidate site interval(s), pad={args.pad})")
    if not merged:
        print("[warn] no candidate regions -- caller should skip BAM subsetting", file=sys.stderr)


if __name__ == "__main__":
    main()
