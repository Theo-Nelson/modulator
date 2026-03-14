#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict, Counter

import pysam


ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def parse_gtf_attrs(attr_field):
    d = {}
    for m in ATTR_RE.finditer(attr_field):
        d[m.group(1)] = m.group(2)
    return d


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


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


def aln_exon_blocks_1based(aln):
    exons = []
    ref = aln.reference_start
    cur_start = ref
    for op, ln in (aln.cigartuples or []):
        if op == 3:  # N
            exons.append((cur_start + 1, ref))
            ref += ln
            cur_start = ref
        elif op in (0, 2, 7, 8):  # M/D/=/X
            ref += ln
        elif op in (1, 4, 5, 6):
            # I / S / H / P do not consume reference
            pass
    exons.append((cur_start + 1, ref))
    return [(s, e) for s, e in exons if e >= s]


def get_read_strand(aln):
    return "-" if aln.is_reverse else "+"


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def load_gene_exon_unions(gtf_path):
    """
    Build per-gene exonic union from the assembler GTF.
    Uses exon lines and groups by (chrom, strand, gene_id).
    """
    gene_exons = defaultdict(list)
    gene_name_map = {}

    with open_text(gtf_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "exon":
                continue
            a = parse_gtf_attrs(attrs)
            gene_id = a.get("gene_id")
            gene_name = a.get("ref_gene_name", a.get("gene_name", gene_id))
            if not gene_id:
                continue
            start = int(start)
            end = int(end)
            gene_exons[(chrom, strand, gene_id)].append((start, end))
            gene_name_map[gene_id] = gene_name

    gene_index = defaultdict(list)
    for (chrom, strand, gene_id), exons in gene_exons.items():
        merged = merge_intervals(exons)
        span_start = merged[0][0]
        span_end = merged[-1][1]
        gene_index[(chrom, strand)].append({
            "gene_id": gene_id,
            "gene_name": gene_name_map.get(gene_id, gene_id),
            "exons": merged,
            "span_start": span_start,
            "span_end": span_end,
        })

    for key in gene_index:
        gene_index[key].sort(key=lambda g: (g["span_start"], g["span_end"], g["gene_id"]))

    return gene_index


def find_overlapping_genes(read_exons, chrom, strand, gene_index, same_strand_only=True):
    """
    Return list of overlapping genes with positive exonic overlap.
    """
    if not read_exons:
        return []

    read_start = min(s for s, e in read_exons)
    read_end = max(e for s, e in read_exons)

    candidates = []
    strands = [strand] if same_strand_only else ["+", "-"]
    for strand_key in strands:
        for g in gene_index.get((chrom, strand_key), []):
            if g["span_start"] > read_end:
                break
            if g["span_end"] < read_start:
                continue
            ov = exon_overlap_len(read_exons, g["exons"])
            if ov > 0:
                candidates.append((g, ov))

    candidates.sort(key=lambda x: (-x[1], x[0]["gene_id"]))
    return candidates


def safe_get_tag(aln, tag, default=None):
    try:
        return aln.get_tag(tag)
    except KeyError:
        return default


def main():
    ap = argparse.ArgumentParser(
        description="Remove reads from a ZT-tagged BAM if they overlap exonic regions of multiple genes."
    )
    ap.add_argument("--bam", required=True, help="Input ZT-tagged BAM")
    ap.add_argument("--gtf", required=True, help="Assembler GTF")
    ap.add_argument("--sample", default=None, help="Optional sample name for scrap quantification output")
    ap.add_argument("--out-clean-bam", required=True, help="Output BAM for reads overlapping exactly one gene")
    ap.add_argument("--out-scrap-bam", required=True, help="Output BAM for reads overlapping multiple genes")
    ap.add_argument("--out-summary-tsv", required=True, help="Summary TSV")
    ap.add_argument("--out-removed-tsv", required=True, help="Detailed TSV of removed reads")
    ap.add_argument("--out-scrap-tx-counts-tsv", required=True, help="Per-sample transcript counts for scrapped assigned reads")
    ap.add_argument(
        "--zero-gene-action",
        choices=["keep", "scrap"],
        default="keep",
        help="What to do with reads overlapping zero genes in the GTF exonic union"
    )
    ap.add_argument(
        "--same-strand-only",
        action="store_true",
        default=True,
        help="Retained for compatibility; overlap checks use the read strand by default"
    )
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_clean_bam) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_scrap_bam) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_summary_tsv) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_removed_tsv) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_scrap_tx_counts_tsv) or ".", exist_ok=True)

    gene_index = load_gene_exon_unions(args.gtf)
    sample = args.sample or os.path.basename(args.bam).replace(".bam", "")

    summary_counts = Counter()
    per_gene_removed = Counter()
    per_zt_removed = Counter()
    per_scrap_tx_removed = Counter()
    per_scrap_tx_meta = {}

    with pysam.AlignmentFile(args.bam, "rb") as inp, \
         pysam.AlignmentFile(args.out_clean_bam, "wb", header=inp.header) as clean_out, \
         pysam.AlignmentFile(args.out_scrap_bam, "wb", header=inp.header) as scrap_out, \
         open(args.out_removed_tsv, "w") as removed_fh:

        removed_fh.write(
            "\t".join([
                "read_name",
                "chrom",
                "strand",
                "zt",
                "zg",
                "zn",
                "n_overlapping_genes",
                "overlapping_gene_ids",
                "overlapping_gene_names",
                "overlap_bases_per_gene",
            ]) + "\n"
        )

        for aln in inp.fetch(until_eof=True):
            if aln.is_unmapped:
                clean_out.write(aln)
                summary_counts["unmapped_kept"] += 1
                continue

            chrom = inp.get_reference_name(aln.reference_id)
            strand = get_read_strand(aln)
            read_exons = aln_exon_blocks_1based(aln)

            overlaps = find_overlapping_genes(
                read_exons,
                chrom,
                strand,
                gene_index,
                same_strand_only=args.same_strand_only,
            )
            n_genes = len(overlaps)

            zt = safe_get_tag(aln, "ZT", "")
            zg = safe_get_tag(aln, "ZG", "")
            zn = safe_get_tag(aln, "ZN", "")

            if n_genes == 1:
                clean_out.write(aln)
                summary_counts["single_gene_kept"] += 1
            elif n_genes == 0:
                if args.zero_gene_action == "keep":
                    clean_out.write(aln)
                    summary_counts["zero_gene_kept"] += 1
                else:
                    scrap_out.write(aln)
                    summary_counts["zero_gene_scrapped"] += 1
            else:
                scrap_out.write(aln)
                summary_counts["multi_gene_scrapped"] += 1

                genes = [g["gene_id"] for g, ov in overlaps]
                names = [g["gene_name"] for g, ov in overlaps]
                ovs = [str(ov) for g, ov in overlaps]

                removed_fh.write(
                    "\t".join([
                        aln.query_name,
                        chrom,
                        strand,
                        str(zt),
                        str(zg),
                        str(zn),
                        str(n_genes),
                        ",".join(genes),
                        ",".join(names),
                        ",".join(ovs),
                    ]) + "\n"
                )

                for gid in genes:
                    per_gene_removed[gid] += 1
                if zt:
                    per_zt_removed[zt] += 1
                    per_scrap_tx_removed[zt] += 1
                    per_scrap_tx_meta[zt] = {
                        "sample": sample,
                        "code": zt,
                        "zg": "" if zg is None else zg,
                        "zn": "" if zn is None else zn,
                    }

    try:
        pysam.index(args.out_clean_bam)
    except Exception:
        pass

    try:
        pysam.index(args.out_scrap_bam)
    except Exception:
        pass

    with open(args.out_summary_tsv, "w") as out:
        out.write("metric\tvalue\n")
        for k in sorted(summary_counts):
            out.write(f"{k}\t{summary_counts[k]}\n")

        out.write("\nremoved_reads_per_gene_id\tcount\n")
        for gid, c in sorted(per_gene_removed.items(), key=lambda x: (-x[1], x[0])):
            out.write(f"{gid}\t{c}\n")

        out.write("\nremoved_reads_per_zt_label\tcount\n")
        for zt, c in sorted(per_zt_removed.items(), key=lambda x: (-x[1], x[0])):
            out.write(f"{zt}\t{c}\n")

    with open(args.out_scrap_tx_counts_tsv, "w") as out:
        out.write("sample\tcode\tzg\tzn\tscrapped_assigned_reads\n")
        for code, count in sorted(per_scrap_tx_removed.items(), key=lambda x: (-x[1], x[0])):
            meta = per_scrap_tx_meta.get(code, {})
            out.write(
                "\t".join([
                    str(meta.get("sample", sample)),
                    str(meta.get("code", code)),
                    str(meta.get("zg", "")),
                    str(meta.get("zn", "")),
                    str(count),
                ]) + "\n"
            )

    print(f"[OK] Wrote clean BAM: {args.out_clean_bam}", file=sys.stderr)
    print(f"[OK] Wrote scrap BAM: {args.out_scrap_bam}", file=sys.stderr)
    print(f"[OK] Wrote summary TSV: {args.out_summary_tsv}", file=sys.stderr)
    print(f"[OK] Wrote removed-read detail TSV: {args.out_removed_tsv}", file=sys.stderr)
    print(f"[OK] Wrote scrap transcript counts TSV: {args.out_scrap_tx_counts_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
