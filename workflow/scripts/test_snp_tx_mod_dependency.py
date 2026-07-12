#!/usr/bin/env python3

import argparse
import json
import os

import pandas as pd

from genotype_utils import (
    benjamini_hochberg,
    binary_rate_delta,
    cmh_test_2x2xk,
    context_key_from_row,
    context_key_from_snp_row,
    load_molecule_mods_for_pairing,
    read_keys_of,
    stream_filter_by_read_keys,
    tsv_header,
)

# Columns the merge / row-builder touches. chrom, start0, end0 and ZT are kept on BOTH sides so
# pandas still emits the "_snp"/"_mod" suffixes the code reads (chrom_snp, start0_mod, ZT_mod).
SNP_USECOLS = [
    "sample", "qname", "snp_id", "chrom", "pos1", "start0", "end0",
    "allele_class", "gene_names", "metagene_indices", "ZT",
]


def parse_args():
    ap = argparse.ArgumentParser(description="Test whether SNP-mod associations persist after conditioning on transcript assignment.")
    ap.add_argument("--molecule-snps", required=True)
    ap.add_argument("--molecule-mods", required=True)
    ap.add_argument("--snp-transcript-assoc", required=True)
    ap.add_argument("--snp-mod-assoc", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-stratum-reads", type=int, default=4)
    return ap.parse_args()

def main():
    args = parse_args()
    snp_tx = pd.read_csv(args.snp_transcript_assoc, sep="\t", low_memory=False) if os.path.exists(args.snp_transcript_assoc) and os.path.getsize(args.snp_transcript_assoc) else pd.DataFrame()
    snp_mod = pd.read_csv(args.snp_mod_assoc, sep="\t", low_memory=False) if os.path.exists(args.snp_mod_assoc) and os.path.getsize(args.snp_mod_assoc) else pd.DataFrame()

    if snp_mod.empty:
        out = pd.DataFrame(columns=[
            "snp_id", "mod_site_id", "n_reads", "n_transcripts_tested", "cmh_stat", "cmh_p_value",
            "cmh_common_odds_ratio", "overall_effect_abs_delta_mod_frac", "weighted_within_tx_effect",
            "classification", "cmh_p_adj_bh"
        ])
        os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
        out.to_csv(args.out_tsv, sep="\t", index=False)
        return

    candidate_pairs = set(zip(snp_mod["snp_id"], snp_mod["mod_site_id"]))
    # Only (snp_id, mod_site_id) pairs from snp_mod_assoc are ever tested, and the merge is an
    # inner join on (sample, qname). So restrict BOTH sides to the ids that appear in a candidate
    # pair, then stream the large SNP table keeping only reads that carry a usable mod call.
    # Every dropped row is one the original code would have discarded anyway.
    cand_snp_ids = {s for s, _ in candidate_pairs}
    cand_mod_ids = {m for _, m in candidate_pairs}

    mod_df = load_molecule_mods_for_pairing(args.molecule_mods, extra_cols=["ZT"])
    if not mod_df.empty:
        mod_df = mod_df[mod_df["mod_site_id"].isin(cand_mod_ids)].copy()
    mod_keys = read_keys_of(mod_df)

    snp_usecols = [c for c in SNP_USECOLS if c in tsv_header(args.molecule_snps)]
    snp_df = stream_filter_by_read_keys(
        args.molecule_snps, snp_usecols, mod_keys,
        row_filter=lambda ch: ch["allele_class"].isin(["ref", "alt"]) & ch["snp_id"].isin(cand_snp_ids),
    )
    if snp_df.empty or mod_df.empty:
        out = pd.DataFrame(columns=[
            "snp_id", "mod_site_id", "n_reads", "n_transcripts_tested", "cmh_stat", "cmh_p_value",
            "cmh_common_odds_ratio", "overall_effect_abs_delta_mod_frac", "weighted_within_tx_effect",
            "classification", "cmh_p_adj_bh"
        ])
        os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
        out.to_csv(args.out_tsv, sep="\t", index=False)
        return

    snp_df["context_key"] = snp_df.apply(context_key_from_snp_row, axis=1)
    mod_df["context_key"] = mod_df.apply(context_key_from_row, axis=1)

    merged = snp_df.merge(mod_df, on=["sample", "qname"], how="inner", suffixes=("_snp", "_mod"))
    merged = merged[
        (merged["chrom_snp"] == merged["chrom_mod"]) &
        (merged["context_key_snp"] == merged["context_key_mod"])
    ].copy()
    merged["transcript_label"] = merged["ZT_mod"].fillna(merged["ZT_snp"]).astype(str)
    merged = merged[merged["transcript_label"].ne("")]
    merged = merged[merged[["snp_id", "mod_site_id"]].apply(tuple, axis=1).isin(candidate_pairs)]

    snp_tx_map = {}
    if not snp_tx.empty and "snp_id" in snp_tx.columns:
        snp_tx_map = snp_tx.set_index("snp_id")[["p_value", "p_adj_bh", "effect_max_abs_tx_frac_diff"]].to_dict("index")
    snp_mod_map = snp_mod.set_index(["snp_id", "mod_site_id"])[["p_value", "p_adj_bh", "effect_abs_delta_mod_frac"]].to_dict("index")

    rows = []
    # sort=True on both groupbys: a deterministic stratum order makes the CMH
    # float summation (cmh_test_2x2xk sums over `strata`) reproducible regardless of
    # upstream (BAM x chrom) shard order; floating-point addition is not associative.
    for (snp_id, mod_site_id), sub in merged.groupby(["snp_id", "mod_site_id"], sort=True):
        strata = []
        stratum_details = []
        for tx, tx_sub in sub.groupby("transcript_label", sort=True):
            grp = tx_sub.groupby(["allele_class", "target_state"]).size()
            a = int(grp.get(("ref", 1), 0))
            b = int(grp.get(("ref", 0), 0))
            c = int(grp.get(("alt", 1), 0))
            d = int(grp.get(("alt", 0), 0))
            total = a + b + c + d
            if total < int(args.min_stratum_reads):
                continue
            if (a + b) == 0 or (c + d) == 0:
                continue
            table = [[a, b], [c, d]]
            strata.append(table)
            stratum_details.append({
                "ZT": tx,
                "ref_modified": a,
                "ref_not_target": b,
                "alt_modified": c,
                "alt_not_target": d,
                "effect": binary_rate_delta(table),
                "n_reads": total,
            })

        if not strata:
            continue

        overall_grp = sub.groupby(["allele_class", "target_state"]).size()
        overall = [
            [int(overall_grp.get(("ref", 1), 0)), int(overall_grp.get(("ref", 0), 0))],
            [int(overall_grp.get(("alt", 1), 0)), int(overall_grp.get(("alt", 0), 0))],
        ]
        cmh_stat, cmh_p_value, cmh_or = cmh_test_2x2xk(strata)
        weighted_effect = 0.0
        total_weight = 0
        for item in stratum_details:
            weighted_effect += item["effect"] * item["n_reads"]
            total_weight += item["n_reads"]
        weighted_effect = round(float(weighted_effect / total_weight), 6) if total_weight > 0 else 0.0

        first = sub.iloc[0]
        snp_tx_meta = snp_tx_map.get(snp_id, {})
        snp_mod_meta = snp_mod_map.get((snp_id, mod_site_id), {})
        rows.append({
            "snp_id": snp_id,
            "mod_site_id": mod_site_id,
            "chrom": first.get("chrom_snp", ""),
            "pos1": int(first.get("pos1", 0)),
            "mod_start0": int(first.get("start0_mod", 0)),
            "target_mod_code": first.get("target_mod_code", ""),
            "n_reads": int(sum(sum(sum(row) for row in table) for table in strata)),
            "n_transcripts_tested": int(len(strata)),
            "cmh_stat": cmh_stat,
            "cmh_p_value": cmh_p_value,
            "cmh_common_odds_ratio": cmh_or,
            "overall_effect_abs_delta_mod_frac": binary_rate_delta(overall),
            "weighted_within_tx_effect": weighted_effect,
            "snp_tx_p_value": snp_tx_meta.get("p_value"),
            "snp_tx_p_adj_bh": snp_tx_meta.get("p_adj_bh"),
            "snp_tx_effect": snp_tx_meta.get("effect_max_abs_tx_frac_diff"),
            "snp_mod_p_value": snp_mod_meta.get("p_value"),
            "snp_mod_p_adj_bh": snp_mod_meta.get("p_adj_bh"),
            "snp_mod_effect": snp_mod_meta.get("effect_abs_delta_mod_frac"),
            "strata_json": json.dumps(sorted(stratum_details, key=lambda d: str(d["ZT"])), separators=(",", ":")),
        })

    out = pd.DataFrame(rows)
    if not out.empty:
        out["cmh_p_adj_bh"] = benjamini_hochberg(out["cmh_p_value"].values)
        classifications = []
        for row in out.itertuples(index=False):
            if row.n_transcripts_tested < 2:
                classifications.append("insufficient_transcript_strata")
            elif row.cmh_p_adj_bh <= 0.05:
                classifications.append("snp_mod_persists_after_transcript_control")
            elif pd.notna(row.snp_mod_p_adj_bh) and row.snp_mod_p_adj_bh <= 0.05 and pd.notna(row.snp_tx_p_adj_bh) and row.snp_tx_p_adj_bh <= 0.05:
                classifications.append("likely_transcript_mediated")
            elif pd.notna(row.snp_mod_p_adj_bh) and row.snp_mod_p_adj_bh <= 0.05:
                classifications.append("snp_mod_not_stable_across_transcripts")
            else:
                classifications.append("unclear")
        out["classification"] = classifications
        out = out.sort_values(
            ["cmh_p_adj_bh", "weighted_within_tx_effect", "snp_id", "mod_site_id"],
            ascending=[True, False, True, True],
        ).reset_index(drop=True)
    else:
        out = pd.DataFrame(columns=[
            "snp_id", "mod_site_id", "chrom", "pos1", "mod_start0", "target_mod_code",
            "n_reads", "n_transcripts_tested", "cmh_stat", "cmh_p_value", "cmh_common_odds_ratio",
            "overall_effect_abs_delta_mod_frac", "weighted_within_tx_effect",
            "snp_tx_p_value", "snp_tx_p_adj_bh", "snp_tx_effect",
            "snp_mod_p_value", "snp_mod_p_adj_bh", "snp_mod_effect",
            "strata_json", "cmh_p_adj_bh", "classification"
        ])

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
