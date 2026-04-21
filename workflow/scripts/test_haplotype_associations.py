#!/usr/bin/env python3

import argparse
import json
import os

import pandas as pd

from genotype_utils import benjamini_hochberg, max_abs_distribution_shift, run_contingency_test


def parse_args():
    ap = argparse.ArgumentParser(description="Test haplotype to transcript and haplotype to mod-site associations.")
    ap.add_argument("--molecule-haplotypes", required=True)
    ap.add_argument("--molecule-mods", required=True)
    ap.add_argument("--out-haplotype-transcript", required=True)
    ap.add_argument("--out-haplotype-mod", required=True)
    ap.add_argument("--min-haplotype-reads", type=int, default=4)
    ap.add_argument("--min-transcript-reads", type=int, default=4)
    ap.add_argument("--min-total-reads", type=int, default=8)
    ap.add_argument("--test", choices=["auto", "fisher", "chi2"], default="auto")
    ap.add_argument("--pseudocount", type=float, default=0.5)
    return ap.parse_args()


def mod_context(row):
    mg = str(row.get("metagene_index", "")).strip()
    if mg and mg.lower() not in {"nan", "none", "null"}:
        return f"MG:{mg}"
    gene = str(row.get("gene_name", "")).strip()
    if gene and gene.lower() not in {"nan", "none", "null"}:
        return f"GENE:{gene}"
    return f"CHR:{row['chrom']}"


def main():
    args = parse_args()
    hap = pd.read_csv(args.molecule_haplotypes, sep="\t", low_memory=False)
    mod = pd.read_csv(args.molecule_mods, sep="\t", low_memory=False)

    if hap.empty:
        pd.DataFrame(columns=[
            "block_id", "context_key", "chrom", "n_reads", "n_haplotypes_tested", "n_transcripts_tested",
            "test_name", "stat_name", "stat_value", "p_value", "effect_max_abs_tx_frac_diff", "per_table_json", "p_adj_bh"
        ]).to_csv(args.out_haplotype_transcript, sep="\t", index=False)
        pd.DataFrame(columns=[
            "block_id", "mod_site_id", "context_key", "chrom", "target_mod_code", "n_reads",
            "n_haplotypes_tested", "test_name", "stat_name", "stat_value", "p_value",
            "effect_max_abs_mod_rate_diff", "per_table_json", "p_adj_bh"
        ]).to_csv(args.out_haplotype_mod, sep="\t", index=False)
        return

    hap = hap[hap["haplotype"].fillna("").astype(str).ne("")].copy()

    tx_rows = []
    tx_df = hap[hap["ZT"].fillna("").astype(str).ne("")].copy()
    for block_id, sub in tx_df.groupby("block_id", sort=False):
        hap_totals = sub.groupby("haplotype").size()
        keep_haps = sorted(hap_totals[hap_totals >= int(args.min_haplotype_reads)].index)
        if len(keep_haps) < 2:
            continue
        sub = sub[sub["haplotype"].isin(keep_haps)].copy()
        tx_totals = sub.groupby("ZT").size()
        keep_tx = sorted(tx_totals[tx_totals >= int(args.min_transcript_reads)].index)
        if len(keep_tx) < 2:
            continue
        sub = sub[sub["ZT"].isin(keep_tx)].copy()
        table = sub.groupby(["haplotype", "ZT"]).size().unstack(fill_value=0)
        tt = table.to_numpy(dtype=float)
        if int(tt.sum()) < int(args.min_total_reads):
            continue
        test_name, stat_name, stat_value, p_value = run_contingency_test(tt, test=args.test, pseudocount=args.pseudocount)
        tx_rows.append({
            "block_id": block_id,
            "context_key": sub.iloc[0].get("context_key", ""),
            "chrom": sub.iloc[0].get("chrom", ""),
            "n_reads": int(tt.sum()),
            "n_haplotypes_tested": int(tt.shape[0]),
            "n_transcripts_tested": int(tt.shape[1]),
            "test_name": test_name,
            "stat_name": stat_name,
            "stat_value": stat_value,
            "p_value": p_value,
            "effect_max_abs_tx_frac_diff": max_abs_distribution_shift(tt),
            "per_table_json": json.dumps({
                "haplotypes": list(table.index),
                "transcripts": list(table.columns),
                "counts": table.values.tolist(),
            }, separators=(",", ":")),
        })

    tx_out = pd.DataFrame(tx_rows)
    if not tx_out.empty:
        tx_out["p_adj_bh"] = benjamini_hochberg(tx_out["p_value"].values)
        tx_out = tx_out.sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False]).reset_index(drop=True)

    if "usable" in mod.columns:
        mod = mod[mod["usable"].fillna(False)].copy()
    else:
        mod = mod[(~mod["fail"].fillna(True)) & mod["within_alignment"].fillna(False)].copy()
    mod = mod[mod["state_detail"].isin(["modified", "canonical", "other_mod"])].copy()
    mod["target_state"] = mod["state_detail"].eq("modified").astype(int)
    mod["context_key"] = mod.apply(mod_context, axis=1)

    mod_rows = []
    if not mod.empty:
        merged = hap.merge(mod, on=["sample", "qname"], how="inner", suffixes=("_hap", "_mod"))
        merged = merged[(merged["context_key_hap"] == merged["context_key_mod"]) | merged["chrom_hap"].eq(merged["chrom_mod"])]
        for (block_id, mod_site_id), sub in merged.groupby(["block_id", "mod_site_id"], sort=False):
            hap_totals = sub.groupby("haplotype").size()
            keep_haps = sorted(hap_totals[hap_totals >= int(args.min_haplotype_reads)].index)
            if len(keep_haps) < 2:
                continue
            sub = sub[sub["haplotype"].isin(keep_haps)].copy()
            table = sub.groupby(["haplotype", "target_state"]).size().unstack(fill_value=0)
            table = table.reindex(columns=[1, 0], fill_value=0)
            tt = table.to_numpy(dtype=float)
            if int(tt.sum()) < int(args.min_total_reads):
                continue
            if table[1].sum() == 0:
                continue
            test_name, stat_name, stat_value, p_value = run_contingency_test(tt, test=args.test, pseudocount=args.pseudocount)
            mod_rows.append({
                "block_id": block_id,
                "mod_site_id": mod_site_id,
                "context_key": sub.iloc[0].get("context_key_hap", ""),
                "chrom": sub.iloc[0].get("chrom_hap", sub.iloc[0].get("chrom_mod", "")),
                "target_mod_code": sub.iloc[0].get("target_mod_code", ""),
                "n_reads": int(tt.sum()),
                "n_haplotypes_tested": int(tt.shape[0]),
                "test_name": test_name,
                "stat_name": stat_name,
                "stat_value": stat_value,
                "p_value": p_value,
                "effect_max_abs_mod_rate_diff": max_abs_distribution_shift(tt),
                "per_table_json": json.dumps({
                    "haplotypes": list(table.index),
                    "states": ["modified", "not_target"],
                    "counts": table.values.tolist(),
                }, separators=(",", ":")),
            })

    mod_out = pd.DataFrame(mod_rows)
    if not mod_out.empty:
        mod_out["p_adj_bh"] = benjamini_hochberg(mod_out["p_value"].values)
        mod_out = mod_out.sort_values(["p_adj_bh", "effect_max_abs_mod_rate_diff"], ascending=[True, False]).reset_index(drop=True)
    else:
        mod_out = pd.DataFrame(columns=[
            "block_id", "mod_site_id", "context_key", "chrom", "target_mod_code", "n_reads",
            "n_haplotypes_tested", "test_name", "stat_name", "stat_value", "p_value",
            "effect_max_abs_mod_rate_diff", "per_table_json", "p_adj_bh"
        ])

    if tx_out.empty:
        tx_out = pd.DataFrame(columns=[
            "block_id", "context_key", "chrom", "n_reads", "n_haplotypes_tested", "n_transcripts_tested",
            "test_name", "stat_name", "stat_value", "p_value", "effect_max_abs_tx_frac_diff", "per_table_json", "p_adj_bh"
        ])

    os.makedirs(os.path.dirname(args.out_haplotype_transcript) or ".", exist_ok=True)
    tx_out.to_csv(args.out_haplotype_transcript, sep="\t", index=False)
    mod_out.to_csv(args.out_haplotype_mod, sep="\t", index=False)


if __name__ == "__main__":
    main()
