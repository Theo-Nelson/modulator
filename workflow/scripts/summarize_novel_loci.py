#!/usr/bin/env python3
"""
summarize_novel_loci.py -- stats for read-backed NOVEL loci (fragmentforms that match no gene in
the reference GTF).

A fragmentform whose annotation found no overlapping reference transcript is tagged
classification="NOVEL_LOCUS" by assemble_transcripts.py and given a unique, coordinate-anchored
gene name (NOVEL_<chrom>_<strand>_<NNN>_<start>_<end>). This script rolls those up into:

  --out-loci           one row per novel locus: span, #fragmentforms, read support, per-sample
                       support, modification sites detected, and (optionally) its splice-junction
                       category.
  --out-fragmentforms  one row per novel fragmentform: exon count, read support, TES, per-sample counts.

Both are empty (header-only) when a run has no novel loci, which is the normal case for a
well-annotated reference restricted to known genes.
"""

import argparse
import os
import sys
from collections import defaultdict


def _attr(attr: str, key: str) -> str:
    needle = f'{key} "'
    i = attr.find(needle)
    if i < 0:
        return ""
    j = attr.find('"', i + len(needle))
    return attr[i + len(needle):j] if j > 0 else ""


def load_novel_from_gtf(gtf_path):
    """zt_label -> dict(gene_name, gene_id, chrom, strand, exons, tes, read_support, zn)"""
    tx = {}
    exons = defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ, attr = f[2], f[8]
            if _attr(attr, "classification") != "NOVEL_LOCUS":
                continue
            zt = _attr(attr, "zt_label")
            if not zt:
                continue
            if typ == "transcript":
                rs = _attr(attr, "read_support")
                tes = _attr(attr, "tes")
                tx[zt] = dict(
                    gene_name=_attr(attr, "ref_gene_name"),
                    gene_id=_attr(attr, "gene_id"),
                    chrom=f[0], strand=f[6],
                    read_support=int(rs) if rs.isdigit() else 0,
                    tes=int(tes) if tes.lstrip("-").isdigit() else 0,
                    zn=_attr(attr, "zn_index"),
                )
            elif typ == "exon":
                exons[zt].append((int(f[3]), int(f[4])))
    for zt, ex in exons.items():
        if zt in tx:
            ex.sort()
            tx[zt]["exons"] = ex
    return {zt: d for zt, d in tx.items() if d.get("exons")}


def load_sample_counts(path):
    """zt_label -> sample_counts string, from *_classification_summary.tsv (optional)."""
    out = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        try:
            i_zt = header.index("zt_label")
        except ValueError:
            return out
        i_sc = header.index("sample_counts") if "sample_counts" in header else -1
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= i_zt:
                continue
            out[p[i_zt]] = p[i_sc] if (0 <= i_sc < len(p)) else ""
    return out


def load_mod_sites_by_gene(zn_long):
    """gene_name -> (set of (chrom,start0,mod_code), set of mod_codes) from the FILTERED long table."""
    sites = defaultdict(set)
    codes = defaultdict(set)
    if not zn_long or not os.path.exists(zn_long) or os.path.getsize(zn_long) == 0:
        return sites, codes
    with open(zn_long) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        need = ("gene_name", "chrom", "start0", "mod_code")
        if not all(c in idx for c in need):
            return sites, codes
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            g = p[idx["gene_name"]]
            sites[g].add((p[idx["chrom"]], p[idx["start0"]], p[idx["mod_code"]]))
            codes[g].add(p[idx["mod_code"]])
    return sites, codes


def load_splice_gene_summary(path):
    """gene_name -> (intron_category, frac_canonical, n_noncanonical)"""
    out = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            out[p[idx["gene_name"]]] = (
                p[idx.get("intron_category", 0)],
                p[idx.get("frac_canonical", 0)],
                p[idx.get("n_noncanonical", 0)],
            )
    return out


def parse_args():
    ap = argparse.ArgumentParser(description="Summarize read-backed novel loci (NOVEL_LOCUS fragmentforms).")
    ap.add_argument("--gtf", required=True, help="Assembled fragmentform GTF")
    ap.add_argument("--classification", default="", help="*_classification_summary.tsv (for per-sample counts)")
    ap.add_argument("--zn-long", default="", help="*_FILTERED_sites_long.tsv (for modification sites)")
    ap.add_argument("--splice-genes", default="", help="*_gene_splice_summary.tsv (for intron category)")
    ap.add_argument("--out-loci", required=True)
    ap.add_argument("--out-fragmentforms", required=True)
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


LOCI_COLS = [
    "locus_name", "chrom", "strand", "start1", "end1", "span_bp",
    "n_fragmentforms", "total_read_support", "n_exons_max",
    "n_mod_sites", "mod_codes", "intron_category", "frac_canonical", "n_noncanonical",
]
FF_COLS = [
    "locus_name", "zt_label", "zn_index", "n_exons", "read_support", "tes", "sample_counts",
]


def main():
    args = parse_args()
    novel = load_novel_from_gtf(args.gtf)
    sample_counts = load_sample_counts(args.classification)
    mod_sites, mod_codes = load_mod_sites_by_gene(args.zn_long)
    splice = load_splice_gene_summary(args.splice_genes)

    by_locus = defaultdict(list)
    for zt, d in novel.items():
        by_locus[d["gene_name"]].append((zt, d))

    os.makedirs(os.path.dirname(args.out_loci) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_fragmentforms) or ".", exist_ok=True)

    with open(args.out_loci, "w") as lo, open(args.out_fragmentforms, "w") as fo:
        lo.write("\t".join(LOCI_COLS) + "\n")
        fo.write("\t".join(FF_COLS) + "\n")

        for locus in sorted(by_locus):
            forms = sorted(by_locus[locus], key=lambda kv: (-kv[1]["read_support"], kv[0]))
            chrom = forms[0][1]["chrom"]
            strand = forms[0][1]["strand"]
            start1 = min(d["exons"][0][0] for _, d in forms)
            end1 = max(d["exons"][-1][1] for _, d in forms)
            total_rs = sum(d["read_support"] for _, d in forms)
            n_ex_max = max(len(d["exons"]) for _, d in forms)
            sites = mod_sites.get(locus, set())
            codes = sorted(mod_codes.get(locus, set()))
            cat, frac_canon, n_noncanon = splice.get(locus, ("", "", ""))

            lo.write("\t".join(str(x) for x in [
                locus, chrom, strand, start1, end1, end1 - start1 + 1,
                len(forms), total_rs, n_ex_max,
                len(sites), ",".join(codes), cat, frac_canon, n_noncanon,
            ]) + "\n")

            for zt, d in forms:
                fo.write("\t".join(str(x) for x in [
                    locus, zt, d.get("zn", ""), len(d["exons"]), d["read_support"], d["tes"],
                    sample_counts.get(zt, ""),
                ]) + "\n")

    n_loci = len(by_locus)
    n_forms = len(novel)
    print(f"[ok] wrote {args.out_loci}: {n_loci} novel locus/loci ({n_forms} fragmentform(s))")
    print(f"[ok] wrote {args.out_fragmentforms}")
    if n_loci == 0 and args.verbose:
        print("[info] no NOVEL_LOCUS fragmentforms in this run "
              "(every assembled fragmentform matched a reference gene).", file=sys.stderr)

    # Uniqueness invariant: every novel locus name must be distinct (the naming fix in
    # assemble_transcripts.py guarantees this). Fail loudly if two loci ever collide.
    names = list(by_locus.keys())
    if len(names) != len(set(names)):
        sys.exit("[error] novel locus names are not unique -- naming invariant violated")


if __name__ == "__main__":
    main()
