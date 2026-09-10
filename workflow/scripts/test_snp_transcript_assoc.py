#!/usr/bin/env python3

import argparse
import json
import os

import numpy as np
import pandas as pd

from genotype_utils import (add_heterogeneity_flag, benjamini_hochberg, informative_strata, max_abs_distribution_shift,
                            open_chrom_table, run_contingency_test, stratified_max_distribution_shift,
                            stratified_primary, stratum_heterogeneity)

# Only these columns are used (grouping: sample/snp_id/allele_class/ZT; per-SNP metadata: the rest).
# The molecule_snps table is ~1.7 GB / 7.5M rows on Huh7 mock, and reading all 21 object-dtype columns
# is what pushes this to ~5.6 GB; loading just these -- with the repeated string columns as
# categoricals -- is far lighter and produces identical output. `sample` is loaded so the test can be
# stratified by sample (see main): pooling reads across replicates lets per-sample allele-rate +
# transcript-composition imbalance manufacture a Simpson's-paradox SNP->transcript association.
WANTED_COLS = [
    "sample", "snp_id", "allele_class", "ZT", "chrom", "pos1", "ref", "alt",
    "gene_names", "gene_ids", "metagene_indices",
]
# Categoricals for the repeated low-cardinality columns. allele_class and ZT stay object:
# read_csv builds the full string column transiently regardless of dtype, so categorical ZT gives
# no read-time saving here and only adds fillna/observed traps. The real cut comes from usecols
# (10 of 21 columns) + these categoricals.
CATEGORICAL = {
    "snp_id": "category", "chrom": "category", "ref": "category", "alt": "category",
    "gene_names": "category", "gene_ids": "category", "metagene_indices": "category",
}


def parse_args():
    ap = argparse.ArgumentParser(description="Test SNP allele to transcript assignment associations.")
    ap.add_argument("--molecule-snps", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-allele-reads", type=int, default=4)
    ap.add_argument("--min-transcript-reads", type=int, default=4)
    ap.add_argument("--test", choices=["auto", "fisher", "chi2"], default="auto")
    ap.add_argument("--pseudocount", type=float, default=0.5)
    return ap.parse_args()


def _rows_for_chrom(df, args, rows):
    """Per-SNP allele x transcript tables built with numpy on integer codes (the previous version ran
    three pandas groupby/pivot_table calls per SNP -- ~30 ms each, 54 min on an 8-sample mosquito
    chromosome set). Same numbers, same SNP order (first appearance), same sorted transcript columns
    and sorted sample strata."""
    if df.empty:
        return
    min_tx = int(args.min_transcript_reads)
    min_al = int(args.min_allele_reads)
    snp_codes, snp_ids = pd.factorize(df["snp_id"].astype(str))      # first-appearance order
    order = np.argsort(snp_codes, kind="stable")
    codes_sorted = snp_codes[order]
    bounds = np.flatnonzero(np.r_[True, codes_sorted[1:] != codes_sorted[:-1], True])
    zt_all = np.asarray(df["ZT"].astype(str).tolist(), dtype=object)
    al_all = (df["allele_class"].astype(str).to_numpy() == "alt").astype(np.int64)
    samp_all = np.asarray(df["sample"].astype(str).tolist(), dtype=object) if "sample" in df.columns else None
    raw = {c: df[c].to_numpy() for c in ("chrom", "pos1", "ref", "alt", "gene_names", "gene_ids", "metagene_indices")
           if c in df.columns}

    for b in range(len(bounds) - 1):
        idx = order[bounds[b]:bounds[b + 1]]
        first_i = idx[0]                       # first row of this SNP in file order (idx is sorted, stable)
        zt_b = zt_all[idx]
        al_b = al_all[idx]
        zt_u, zt_inv = np.unique(zt_b, return_inverse=True)    # sorted -> keep_tx sorted, as before
        tot = np.bincount(zt_inv, minlength=zt_u.size)
        keepm = tot >= min_tx
        if int(keepm.sum()) < 2:
            continue
        keep_tx = [str(z) for z in zt_u[keepm]]
        col = np.full(zt_u.size, -1, dtype=np.int64)
        col[keepm] = np.arange(int(keepm.sum()))
        c = col[zt_inv]
        m = c >= 0
        K = len(keep_tx)
        table = np.zeros((2, K), dtype=np.int64)
        np.add.at(table, (al_b[m], c[m]), 1)
        n_ref = int(table[0].sum()); n_alt = int(table[1].sum())
        if n_ref < min_al or n_alt < min_al:
            continue
        tt = table.astype(float)
        # POOLED test (kept as *_pooled): pools reads across samples -> confounded by replicate.
        pooled_name, pooled_stat_name, pooled_stat, pooled_p = run_contingency_test(
            tt, test=args.test, pseudocount=args.pseudocount)
        # SAMPLE-STRATIFIED CMH (primary): one 2 x len(keep_tx) allele x transcript table per sample
        # (sorted sample order; only samples with reads on a kept transcript, as before).
        strata = []
        if samp_all is not None:
            s_b = samp_all[idx][m]; al_m = al_b[m]; c_m = c[m]
            s_u, s_inv = np.unique(s_b, return_inverse=True)
            for k in range(s_u.size):
                sel = s_inv == k
                T = np.zeros((2, K), dtype=np.int64)
                np.add.at(T, (al_m[sel], c_m[sel]), 1)
                strata.append(T.astype(float))
        inf = informative_strata(strata)
        test_name, stat_name, stat_value, p_value, n_strata, _mode = stratified_primary(
            inf, lambda T: run_contingency_test(T, test=args.test, pseudocount=args.pseudocount))
        eff_strat = stratified_max_distribution_shift(inf) if _mode != "none" else float("nan")
        _hstat, het_p, _hdf, _ = stratum_heterogeneity(inf)
        strata_heterogeneous = bool(np.isfinite(het_p) and het_p < 0.05)
        per_tx = [{"ZT": tx, "ref_reads": int(table[0, jx]), "alt_reads": int(table[1, jx])}
                  for jx, tx in enumerate(keep_tx)]
        rows.append({
            "snp_id": snp_ids[b],
            "chrom": raw["chrom"][first_i] if "chrom" in raw else "",
            "pos1": int(raw["pos1"][first_i]) if "pos1" in raw else 0,
            "ref": raw["ref"][first_i] if "ref" in raw else "",
            "alt": raw["alt"][first_i] if "alt" in raw else "",
            "gene_names": raw["gene_names"][first_i] if "gene_names" in raw else "",
            "gene_ids": raw["gene_ids"][first_i] if "gene_ids" in raw else "",
            "metagene_indices": raw["metagene_indices"][first_i] if "metagene_indices" in raw else "",
            "n_reads": int(tt.sum()),
            "n_ref_reads": n_ref,
            "n_alt_reads": n_alt,
            "n_transcripts_tested": int(tt.shape[1]),
            "n_strata_informative": int(n_strata),
            "strata_heterogeneous": bool(strata_heterogeneous),
            "strata_heterogeneity_p": round(het_p, 6) if np.isfinite(het_p) else float("nan"),
            "test_name": test_name,
            "stat_name": stat_name,
            "stat_value": stat_value,
            "p_value": p_value,
            "effect_max_abs_tx_frac_diff": eff_strat,
            "test_name_pooled": pooled_name,
            "stat_value_pooled": pooled_stat,
            "p_value_pooled": pooled_p,
            "effect_max_abs_tx_frac_diff_pooled": max_abs_distribution_shift(tt),
            "per_transcript_json": json.dumps(per_tx, separators=(",", ":")),
        })


def main():
    args = parse_args()
    # One chromosome at a time (a byte range of the chrom-sorted table; 71 GiB when loaded whole on an
    # 8-sample mosquito run). The table is written in sorted-chrom order with rows in position order, so
    # concatenating per-chrom rows in that order reproduces the old first-appearance group order exactly.
    tbl = open_chrom_table(args.molecule_snps)
    usecols = [c for c in WANTED_COLS if c in tbl.header_cols]
    dtype = {c: t for c, t in CATEGORICAL.items() if c in usecols}
    rows = []
    try:
        for chrom in tbl.chroms:
            df = tbl.read(chrom, usecols=usecols, dtype=dtype)
            if df.empty:
                continue
            keep = df["allele_class"].isin(["ref", "alt"]) & df["ZT"].fillna("").astype(str).ne("")
            df = df.loc[keep].copy()
            if df.empty:
                continue
            _rows_for_chrom(df, args, rows)
            del df
    finally:
        if hasattr(tbl, "close"):
            tbl.close()

    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
        out = add_heterogeneity_flag(out)          # BH-adjust the heterogeneity flag like every other p
        out = out.sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False]).reset_index(drop=True)
    else:
        out = pd.DataFrame(columns=[
            "snp_id", "chrom", "pos1", "ref", "alt", "gene_names", "gene_ids", "metagene_indices",
            "n_reads", "n_ref_reads", "n_alt_reads", "n_transcripts_tested", "n_strata_informative",
            "strata_heterogeneous", "strata_heterogeneity_p", "strata_heterogeneity_p_adj",
            "test_name", "stat_name", "stat_value", "p_value", "effect_max_abs_tx_frac_diff",
            "test_name_pooled", "stat_value_pooled", "p_value_pooled", "effect_max_abs_tx_frac_diff_pooled",
            "per_transcript_json", "p_adj_bh"
        ])

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
