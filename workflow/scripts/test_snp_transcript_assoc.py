#!/usr/bin/env python3

import argparse
import json
import os

import pandas as pd

from genotype_utils import (benjamini_hochberg, cmh_stratified_test, max_abs_distribution_shift,
                            run_contingency_test, stratified_max_distribution_shift, tsv_header)

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


def main():
    args = parse_args()
    header = tsv_header(args.molecule_snps)
    usecols = [c for c in WANTED_COLS if c in header]
    dtype = {c: t for c, t in CATEGORICAL.items() if c in usecols}
    df = pd.read_csv(args.molecule_snps, sep="\t", usecols=usecols, dtype=dtype, low_memory=False)
    keep = df["allele_class"].isin(["ref", "alt"]) & df["ZT"].fillna("").astype(str).ne("")
    df = df.loc[keep].copy()

    rows = []
    for snp_id, sub in df.groupby("snp_id", sort=False, observed=True):
        grp = sub.groupby(["allele_class", "ZT"], as_index=False, observed=True).size()
        tx_totals = grp.groupby("ZT", observed=True)["size"].sum()
        keep_tx = sorted(tx_totals[tx_totals >= int(args.min_transcript_reads)].index)
        if len(keep_tx) < 2:
            continue
        grp = grp[grp["ZT"].isin(keep_tx)].copy()
        table = (
            grp.pivot_table(index="allele_class", columns="ZT", values="size", fill_value=0, aggfunc="sum")
               .reindex(index=["ref", "alt"], columns=keep_tx, fill_value=0)
        )
        allele_totals = table.sum(axis=1)
        if allele_totals.get("ref", 0) < int(args.min_allele_reads):
            continue
        if allele_totals.get("alt", 0) < int(args.min_allele_reads):
            continue
        tt = table.to_numpy(dtype=float)
        # POOLED test (kept as *_pooled): pools reads across samples -> confounded by replicate.
        pooled_name, pooled_stat_name, pooled_stat, pooled_p = run_contingency_test(
            tt, test=args.test, pseudocount=args.pseudocount)
        # SAMPLE-STRATIFIED CMH (primary): one 2 x len(keep_tx) allele x transcript table per sample,
        # combined by the generalized (Landis-Koch) CMH general-association statistic. Transcript column
        # order is FIXED (keep_tx) across strata; a sample lacking both alleles or >=2 transcripts is
        # dropped as uninformative. With only ONE sample the CMH statistic is an uncorrected asymptotic
        # chi2 (anti-conservative at small counts), so a lone sample keeps the exact pooled test instead.
        strata = []
        if "sample" in sub.columns and sub["sample"].nunique() > 1:
            sgrp = sub.groupby(["sample", "allele_class", "ZT"], as_index=False, observed=True).size()
            sgrp = sgrp[sgrp["ZT"].isin(keep_tx)]
            for _samp, ss in sgrp.groupby("sample", observed=True):
                st = (ss.pivot_table(index="allele_class", columns="ZT", values="size",
                                     fill_value=0, aggfunc="sum")
                        .reindex(index=["ref", "alt"], columns=keep_tx, fill_value=0))
                strata.append(st.to_numpy(dtype=float))
        if strata:
            test_name, stat_name, stat_value, p_value, n_strata = cmh_stratified_test(strata)
        else:
            test_name, stat_name, stat_value, p_value, n_strata = pooled_name, pooled_stat_name, pooled_stat, pooled_p, 0
        # sample-stratified effect: per transcript, coverage-weighted mean over samples of
        # (ref_frac - alt_frac); report the max |.| . Falls back to the pooled shift if no strata.
        eff_strat = stratified_max_distribution_shift(strata) if strata else max_abs_distribution_shift(tt)
        per_tx = []
        for tx in table.columns:
            per_tx.append({
                "ZT": tx,
                "ref_reads": int(table.loc["ref", tx]),
                "alt_reads": int(table.loc["alt", tx]),
            })
        first = sub.iloc[0]
        rows.append({
            "snp_id": snp_id,
            "chrom": first.get("chrom", ""),
            "pos1": int(first.get("pos1", 0)),
            "ref": first.get("ref", ""),
            "alt": first.get("alt", ""),
            "gene_names": first.get("gene_names", ""),
            "gene_ids": first.get("gene_ids", ""),
            "metagene_indices": first.get("metagene_indices", ""),
            "n_reads": int(tt.sum()),
            "n_ref_reads": int(allele_totals.get("ref", 0)),
            "n_alt_reads": int(allele_totals.get("alt", 0)),
            "n_transcripts_tested": int(tt.shape[1]),
            "n_strata_informative": int(n_strata),
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

    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
        out = out.sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False]).reset_index(drop=True)
    else:
        out = pd.DataFrame(columns=[
            "snp_id", "chrom", "pos1", "ref", "alt", "gene_names", "gene_ids", "metagene_indices",
            "n_reads", "n_ref_reads", "n_alt_reads", "n_transcripts_tested", "n_strata_informative",
            "test_name", "stat_name", "stat_value", "p_value", "effect_max_abs_tx_frac_diff",
            "test_name_pooled", "stat_value_pooled", "p_value_pooled", "effect_max_abs_tx_frac_diff_pooled",
            "per_transcript_json", "p_adj_bh"
        ])

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
