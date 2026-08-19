#!/usr/bin/env python3
"""
classify_splice_junctions.py -- canonical vs non-canonical splice junctions of the assembled
read-backed fragmentforms.

Every intron of every fragmentform in the assembled GTF is looked up in the reference genome and
classified by its donor/acceptor dinucleotides (in TRANSCRIPT orientation, i.e. reverse-complemented
on the minus strand):

    GT..AG   CANONICAL_GT_AG      major (U2) spliceosome, ~98-99% of human introns
    GC..AG   SEMI_CANONICAL_GC_AG major spliceosome, non-canonical donor (~0.5-1%)
    AT..AC   MINOR_AT_AC          minor (U12) spliceosome
    other    NONCANONICAL         everything else (mis-assembly, alignment artifact, or real oddity)

Because the fragmentform intron chains are read-derived, these are the junctions the reads actually
support -- so this is "canonical vs non-canonical junction usage within the samples", not an
annotation lookup.

Outputs
-------
--out-junctions : one row per (fragmentform, intron)  -- donor/acceptor, motif, class, support
--out-genes     : one row per gene -- junction counts by class, frac_canonical, has_noncanonical,
                  and an intron_category summarising the gene's junction repertoire.

Intron categories (per gene):
    ALL_CANONICAL          every distinct junction is GT-AG
    CANONICAL_WITH_GC_AG   only GT-AG + GC-AG (both major spliceosome)
    HAS_MINOR_U12          at least one AT-AC junction
    HAS_NONCANONICAL       at least one junction outside {GT-AG, GC-AG, AT-AC}
    NO_JUNCTIONS           single-exon fragmentforms only
(HAS_NONCANONICAL takes precedence over HAS_MINOR_U12 when both apply.)
"""

import argparse
import os
import sys
from collections import defaultdict

import pysam


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

MOTIF_CLASS = {
    ("GT", "AG"): "CANONICAL_GT_AG",
    ("GC", "AG"): "SEMI_CANONICAL_GC_AG",
    ("AT", "AC"): "MINOR_AT_AC",
}

JUNCTION_CLASSES = ["CANONICAL_GT_AG", "SEMI_CANONICAL_GC_AG", "MINOR_AT_AC", "NONCANONICAL"]


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def _attr(attr: str, key: str) -> str:
    needle = f'{key} "'
    i = attr.find(needle)
    if i < 0:
        return ""
    j = attr.find('"', i + len(needle))
    return attr[i + len(needle):j] if j > 0 else ""


def load_fragmentforms(gtf_path):
    """(gene_name, zn) -> dict(chrom, strand, exons=[(s1,e1)...], zt_label, read_support, gene_id).

    Exons are 1-based inclusive, as written by assemble_transcripts.py.
    """
    exons = defaultdict(list)
    meta = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ, attr = f[2], f[8]
            gname = _attr(attr, "ref_gene_name") or _attr(attr, "gene_id")
            zn = _attr(attr, "zn_index") or _attr(attr, "transcript_index")
            if not gname or not zn:
                continue
            key = (gname, zn)
            if typ == "exon":
                exons[key].append((int(f[3]), int(f[4])))
            elif typ == "transcript":
                rs = _attr(attr, "read_support")
                meta[key] = dict(
                    chrom=f[0], strand=f[6],
                    zt_label=_attr(attr, "zt_label"),
                    gene_id=_attr(attr, "gene_id"),
                    read_support=int(rs) if rs.isdigit() else 0,
                )
    out = {}
    for key, exl in exons.items():
        exl.sort()
        m = meta.get(key, {})
        out[key] = dict(chrom=m.get("chrom", ""), strand=m.get("strand", "+"), exons=exl,
                        zt_label=m.get("zt_label", ""), gene_id=m.get("gene_id", ""),
                        read_support=m.get("read_support", 0))
    return out


def introns_of(exons):
    """Genomic introns as 1-based inclusive (start, end) between consecutive exons."""
    return [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)
            if exons[i + 1][0] - 1 >= exons[i][1] + 1]


def junction_motif(fa, chrom, i_start1, i_end1, strand):
    """Return (donor, acceptor) dinucleotides in TRANSCRIPT orientation.

    Intron 1-based inclusive [i_start1, i_end1] -> 0-based half-open [i_start1-1, i_end1).
    On '+', donor = first 2 intron bases, acceptor = last 2.
    On '-', the transcript reads the genome right-to-left, so the donor is at the HIGHER genomic
    coordinate: donor = revcomp(last 2), acceptor = revcomp(first 2).
    """
    s0 = i_start1 - 1
    e0 = i_end1                      # exclusive
    if e0 - s0 < 4:                  # too short to have distinct donor+acceptor dinucleotides
        return "", ""
    try:
        first2 = fa.fetch(chrom, s0, s0 + 2).upper()
        last2 = fa.fetch(chrom, e0 - 2, e0).upper()
    except (KeyError, ValueError):
        return "", ""
    if len(first2) < 2 or len(last2) < 2:
        return "", ""
    if strand == "-":
        return revcomp(last2), revcomp(first2)
    return first2, last2


def classify_motif(donor, acceptor):
    if not donor or not acceptor:
        return "NONCANONICAL"
    return MOTIF_CLASS.get((donor, acceptor), "NONCANONICAL")


def gene_intron_category(counts, n_junctions):
    if n_junctions == 0:
        return "NO_JUNCTIONS"
    if counts["NONCANONICAL"] > 0:
        return "HAS_NONCANONICAL"
    if counts["MINOR_AT_AC"] > 0:
        return "HAS_MINOR_U12"
    if counts["SEMI_CANONICAL_GC_AG"] > 0:
        return "CANONICAL_WITH_GC_AG"
    return "ALL_CANONICAL"


def parse_args():
    ap = argparse.ArgumentParser(description="Classify canonical vs non-canonical splice junctions of assembled fragmentforms.")
    ap.add_argument("--gtf", required=True, help="Assembled fragmentform GTF")
    ap.add_argument("--reference-fa", required=True, help="Reference FASTA (indexed)")
    ap.add_argument("--out-junctions", required=True)
    ap.add_argument("--out-genes", required=True)
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def main():
    args = parse_args()
    iso = load_fragmentforms(args.gtf)
    if args.verbose:
        print(f"[sj] fragmentforms: {len(iso)}", file=sys.stderr)

    fa = pysam.FastaFile(args.reference_fa)

    # Preflight: if the assembled-GTF contigs are not in the FASTA index (e.g. GTF "1" vs FASTA
    # "chr1"), every fa.fetch below fails silently and EVERY junction is mislabelled NONCANONICAL
    # with exit 0. Warn loudly rather than emit a corrupt QC table.
    fa_contigs = set(fa.references)
    gtf_contigs = {d["chrom"] for d in iso.values()}
    unknown = sorted(gtf_contigs - fa_contigs)
    if unknown:
        print(f"[sj][WARNING] {len(unknown)}/{len(gtf_contigs)} GTF contig(s) are NOT in the "
              f"reference FASTA index (e.g. {unknown[:3]}); their splice junctions will all be "
              f"classified NONCANONICAL. Check that the reference FASTA matches the alignment "
              f"(contig naming).", file=sys.stderr)

    jrows = []
    # gene -> distinct junction (chrom,start,end,strand) -> class
    gene_junctions = defaultdict(dict)
    gene_forms = defaultdict(set)

    for (gname, zn), d in sorted(iso.items()):
        gene_forms[gname].add(zn)
        chrom, strand = d["chrom"], d["strand"]
        for (i_s, i_e) in introns_of(d["exons"]):
            donor, acceptor = junction_motif(fa, chrom, i_s, i_e, strand)
            cls = classify_motif(donor, acceptor)
            motif = f"{donor}-{acceptor}" if donor and acceptor else "NA"
            jrows.append([
                gname, d["gene_id"], d["zt_label"], zn, chrom, strand,
                i_s, i_e, i_e - i_s + 1, donor, acceptor, motif, cls, d["read_support"],
            ])
            gene_junctions[gname][(chrom, i_s, i_e, strand)] = cls

    fa.close()

    os.makedirs(os.path.dirname(args.out_junctions) or ".", exist_ok=True)
    with open(args.out_junctions, "w") as out:
        out.write("\t".join([
            "gene_name", "gene_id", "zt_label", "zn_index", "chrom", "strand",
            "intron_start1", "intron_end1", "intron_length", "donor", "acceptor",
            "motif", "junction_class", "read_support",
        ]) + "\n")
        for r in jrows:
            out.write("\t".join(str(x) for x in r) + "\n")

    # Per-gene summary over DISTINCT junctions (a junction shared by several fragmentforms counts once)
    os.makedirs(os.path.dirname(args.out_genes) or ".", exist_ok=True)
    grows = []
    for gname in sorted(set(list(gene_junctions.keys()) + list(gene_forms.keys()))):
        juncs = gene_junctions.get(gname, {})
        counts = {c: 0 for c in JUNCTION_CLASSES}
        for cls in juncs.values():
            counts[cls] += 1
        n_j = len(juncs)
        n_canon = counts["CANONICAL_GT_AG"]
        frac_canon = round(n_canon / n_j, 6) if n_j else 0.0
        n_noncanon = counts["NONCANONICAL"]
        grows.append([
            gname, len(gene_forms.get(gname, ())), n_j,
            counts["CANONICAL_GT_AG"], counts["SEMI_CANONICAL_GC_AG"],
            counts["MINOR_AT_AC"], n_noncanon,
            frac_canon, int(n_noncanon > 0), gene_intron_category(counts, n_j),
        ])
    with open(args.out_genes, "w") as out:
        out.write("\t".join([
            "gene_name", "n_fragmentforms", "n_distinct_junctions",
            "n_canonical_GT_AG", "n_semi_canonical_GC_AG", "n_minor_AT_AC", "n_noncanonical",
            "frac_canonical", "has_noncanonical", "intron_category",
        ]) + "\n")
        for r in grows:
            out.write("\t".join(str(x) for x in r) + "\n")

    tot = len(jrows)
    by_cls = defaultdict(int)
    for r in jrows:
        by_cls[r[12]] += 1
    print(f"[ok] wrote {args.out_junctions}: {tot} fragmentform intron(s) over {len(grows)} gene(s)")
    for c in JUNCTION_CLASSES:
        n = by_cls.get(c, 0)
        if n:
            print(f"    {c:<24} {n:>7}  ({100.0 * n / tot:.2f}%)")
    print(f"[ok] wrote {args.out_genes}")


if __name__ == "__main__":
    main()
