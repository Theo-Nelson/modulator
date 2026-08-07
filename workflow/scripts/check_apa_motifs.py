#!/usr/bin/env python3
"""Polyadenylation-signal (PAS) motif check for every detected APA site.

Each assembled fragmentform's TES is a cleavage/polyadenylation site. For each one this pulls the
genomic sequence around it *in transcript orientation* and asks:

  * is there a PAS hexamer upstream (AATAAA canonical, or one of the 11 known variants)?
  * is the downstream element U-/GU-rich, as a real cleavage site should be?
  * is the genomic sequence immediately downstream A-RICH?

The last one matters: an oligo-dT primer that anneals to an internal genomic A-stretch produces a
false 3' end (internal priming) -- the classic direct-RNA / 3'-end artifact. A site with NO PAS and
an A-rich downstream genome is very likely that artifact rather than a real APA site, so the motif
check doubles as an artifact filter. Cross-check it against the measured dorado poly(A) tail
(results/polya): a genuine site carries a real tail, an internally-primed one should not.

Sequence handling mirrors classify_splice_junctions.py (pysam FastaFile + strand-aware revcomp).
"""
import argparse
import os
import sys

import pandas as pd
import pysam

# Canonical first, then the standard variants (Beaudoing et al.), roughly by prevalence.
PAS_CANONICAL = "AATAAA"
PAS_VARIANTS = ["ATTAAA", "TATAAA", "AGTAAA", "AAGAAA", "AATATA", "AATACA",
                "CATAAA", "GATAAA", "AATGAA", "ACTAAA", "AATAGA"]

OUT_COLS = ["zt_label", "gene_name", "chrom", "strand", "tes", "fragmentform_class", "read_support",
            "apa_motif_class", "pas_motif", "pas_distance_nt", "downstream_a_frac",
            "downstream_u_frac", "downstream_gu_frac", "upstream_seq", "downstream_seq"]

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def parse_args():
    ap = argparse.ArgumentParser(description="PAS motif check for every APA site (fragmentform TES).")
    ap.add_argument("--classification-summary", required=True, help="*_classification_summary.tsv (has zt_label/chrom/strand/iso_tes)")
    ap.add_argument("--reference-fa", required=True, help="Reference FASTA (indexed)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--upstream", type=int, default=60, help="nt upstream of the cleavage site to scan for a PAS")
    ap.add_argument("--downstream", type=int, default=40, help="nt downstream to score U/GU-richness")
    ap.add_argument("--pas-max-distance", type=int, default=40, help="max nt from hexamer end to cleavage site to count as that site's PAS")
    ap.add_argument("--internal-priming-a-frac", type=float, default=0.65, help="downstream A fraction above which a PAS-less site is called internal priming")
    ap.add_argument("--internal-priming-window", type=int, default=20, help="nt downstream used for the A-richness test")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def _windows(fa, chrom, tes_1based, strand, up, down):
    """Return (upstream_seq, downstream_seq) in TRANSCRIPT orientation.

    upstream_seq is 5'->3' ending immediately before the cleavage site; downstream_seq starts
    immediately after it. Returns (None, None) if the window runs off the contig.
    """
    pos0 = int(tes_1based) - 1
    try:
        clen = fa.get_reference_length(chrom)
    except Exception:
        return None, None
    if strand == "+":
        u0, u1 = pos0 - up, pos0
        d0, d1 = pos0 + 1, pos0 + 1 + down
        if u0 < 0 or d1 > clen:
            return None, None
        return fa.fetch(chrom, u0, u1).upper(), fa.fetch(chrom, d0, d1).upper()
    # minus strand: transcript runs right-to-left, so genomic-right is upstream
    u0, u1 = pos0 + 1, pos0 + 1 + up
    d0, d1 = pos0 - down, pos0
    if d0 < 0 or u1 > clen:
        return None, None
    return revcomp(fa.fetch(chrom, u0, u1).upper()), revcomp(fa.fetch(chrom, d0, d1).upper())


def _find_pas(upstream_seq, max_dist):
    """Closest PAS hexamer to the cleavage site. Canonical wins ties. -> (motif, distance) or (None, None).

    upstream_seq ends at the cleavage site, so distance = nt from the hexamer's 3' end to the site.
    """
    best = None
    for motif in [PAS_CANONICAL] + PAS_VARIANTS:
        start = 0
        while True:
            i = upstream_seq.find(motif, start)
            if i < 0:
                break
            dist = len(upstream_seq) - (i + len(motif))
            if dist <= max_dist:
                is_canon = motif == PAS_CANONICAL
                # prefer canonical, then the hexamer closest to the cleavage site
                key = (0 if is_canon else 1, dist)
                if best is None or key < best[0]:
                    best = (key, motif, dist)
            start = i + 1
    return (best[1], best[2]) if best else (None, None)


def _frac(seq, bases):
    return round(sum(seq.count(b) for b in bases) / len(seq), 4) if seq else 0.0


def main():
    args = parse_args()
    df = pd.read_csv(args.classification_summary, sep="\t", low_memory=False)
    df.columns = [str(c).lstrip("#") for c in df.columns]
    need = {"zt_label", "chrom", "strand", "iso_tes"}
    if df.empty or not need.issubset(df.columns):
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        print(f"[apa_motifs] no usable classification summary ({need - set(df.columns)}); wrote header-only",
              file=sys.stderr, flush=True)
        return

    fa = pysam.FastaFile(args.reference_fa)
    gene_col = "gtf_gene_name" if "gtf_gene_name" in df.columns else ("gene_name" if "gene_name" in df.columns else None)
    rows = []
    for r in df.itertuples(index=False):
        chrom, strand, tes = str(r.chrom), str(r.strand), r.iso_tes
        if pd.isna(tes) or strand not in ("+", "-"):
            continue
        up_seq, down_seq = _windows(fa, chrom, tes, strand, args.upstream, args.downstream)
        if up_seq is None:
            continue
        motif, dist = _find_pas(up_seq, args.pas_max_distance)
        a_frac = _frac(down_seq[:args.internal_priming_window], "A")
        if motif == PAS_CANONICAL:
            cls = "PAS_CANONICAL"
        elif motif:
            cls = "PAS_VARIANT"
        elif a_frac >= args.internal_priming_a_frac:
            cls = "PAS_NONE_INTERNAL_PRIMING"
        else:
            cls = "PAS_NONE"
        rows.append({
            "zt_label": r.zt_label, "gene_name": getattr(r, gene_col, "") if gene_col else "",
            "chrom": chrom, "strand": strand, "tes": int(tes),
            "fragmentform_class": getattr(r, "classification", ""),
            "read_support": getattr(r, "read_support", ""),
            "apa_motif_class": cls, "pas_motif": motif or "", "pas_distance_nt": dist if dist is not None else "",
            "downstream_a_frac": a_frac,
            "downstream_u_frac": _frac(down_seq, "T"), "downstream_gu_frac": _frac(down_seq, "GT"),
            "upstream_seq": up_seq, "downstream_seq": down_seq,
        })
    fa.close()

    out = pd.DataFrame(rows)
    if out.empty:
        out = pd.DataFrame(columns=OUT_COLS)
    else:
        out = out[OUT_COLS].sort_values(["apa_motif_class", "gene_name", "tes"]).reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    if args.verbose and not out.empty:
        vc = out["apa_motif_class"].value_counts()
        tot = int(vc.sum())
        print(f"[apa_motifs] {tot} APA sites -> " +
              ", ".join(f"{k}={v} ({100*v/tot:.1f}%)" for k, v in vc.items()), flush=True)


if __name__ == "__main__":
    main()
