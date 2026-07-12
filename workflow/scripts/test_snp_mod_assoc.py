#!/usr/bin/env python3

import argparse
import json
import os

import pandas as pd

from genotype_utils import (
    benjamini_hochberg,
    binary_rate_delta,
    context_key_from_row,
    context_key_from_snp_row,
    load_molecule_mods_for_pairing,
    read_keys_of,
    run_contingency_test,
    stream_filter_by_read_keys,
    tsv_header,
)


def parse_args():
    ap = argparse.ArgumentParser(description="Test SNP allele to mod-site associations on the same molecules.")
    ap.add_argument("--molecule-snps", required=True)
    ap.add_argument("--molecule-mods", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-allele-reads", type=int, default=4)
    ap.add_argument("--min-total-reads", type=int, default=8)
    ap.add_argument("--test", choices=["auto", "fisher", "chi2"], default="auto")
    ap.add_argument("--pseudocount", type=float, default=0.5)
    return ap.parse_args()

SNP_USECOLS = [
    # Every column the merge/row-builder touches. chrom/start0/end0 are kept even though the SNP
    # side does not use them directly, because they must still COLLIDE with the mod side so pandas
    # emits the "_snp"/"_mod" suffixes the row builder reads (e.g. start0_mod, chrom_snp).
    "sample", "qname", "snp_id", "chrom", "pos1", "start0", "end0",
    "allele_class", "gene_names", "metagene_indices",
]


def main():
    args = parse_args()

    # Load the SMALL mod table first, filter it, and take its read keys. The merge below is an
    # inner join on (sample, qname), so a SNP row whose read has no usable mod call can never
    # survive it -- streaming molecule_snps and keeping only matching reads is therefore exactly
    # lossless, and avoids materializing the whole (multi-GB, 7.5M-row) SNP table.
    mod_df = load_molecule_mods_for_pairing(args.molecule_mods)
    mod_keys = read_keys_of(mod_df)

    snp_usecols = [c for c in SNP_USECOLS if c in tsv_header(args.molecule_snps)]
    snp_df = stream_filter_by_read_keys(
        args.molecule_snps, snp_usecols, mod_keys,
        row_filter=lambda ch: ch["allele_class"].isin(["ref", "alt"]),
    )

    if snp_df.empty or mod_df.empty:
        out = pd.DataFrame(columns=[
            "snp_id", "mod_site_id", "chrom", "pos1", "mod_start0", "mod_end0", "target_mod_code",
            "n_reads", "n_ref_reads", "n_alt_reads", "n_modified", "n_not_target",
            "test_name", "stat_name", "stat_value", "p_value", "effect_abs_delta_mod_frac",
            "per_state_json", "p_adj_bh"
        ])
        os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
        out.to_csv(args.out_tsv, sep="\t", index=False)
        return

    snp_df["context_key"] = snp_df.apply(context_key_from_snp_row, axis=1)
    mod_df["context_key"] = mod_df.apply(context_key_from_row, axis=1)

    merged = snp_df.merge(
        mod_df,
        on=["sample", "qname"],
        how="inner",
        suffixes=("_snp", "_mod"),
    )
    merged = merged[(merged["chrom_snp"] == merged["chrom_mod"]) & (merged["context_key_snp"] == merged["context_key_mod"])].copy()
    if merged.empty:
        out = pd.DataFrame(columns=[
            "snp_id", "mod_site_id", "chrom", "pos1", "mod_start0", "mod_end0", "target_mod_code",
            "n_reads", "n_ref_reads", "n_alt_reads", "n_modified", "n_not_target",
            "test_name", "stat_name", "stat_value", "p_value", "effect_abs_delta_mod_frac",
            "per_state_json", "p_adj_bh"
        ])
        os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
        out.to_csv(args.out_tsv, sep="\t", index=False)
        return

    rows = []
    for (snp_id, mod_site_id), sub in merged.groupby(["snp_id", "mod_site_id"], sort=False):
        grp = sub.groupby(["allele_class", "target_state"], as_index=False).size()
        table = grp.pivot_table(index="allele_class", columns="target_state", values="size", fill_value=0, aggfunc="sum")
        table = table.reindex(index=["ref", "alt"], fill_value=0).reindex(columns=[1, 0], fill_value=0)
        allele_totals = table.sum(axis=1)
        if allele_totals.get("ref", 0) < int(args.min_allele_reads):
            continue
        if allele_totals.get("alt", 0) < int(args.min_allele_reads):
            continue
        if int(table.to_numpy().sum()) < int(args.min_total_reads):
            continue
        tt = table.to_numpy(dtype=float)
        test_name, stat_name, stat_value, p_value = run_contingency_test(tt, test=args.test, pseudocount=args.pseudocount)
        first = sub.iloc[0]
        rows.append({
            "snp_id": snp_id,
            "mod_site_id": mod_site_id,
            "chrom": first.get("chrom_snp", ""),
            "pos1": int(first.get("pos1", 0)),
            "mod_start0": int(first.get("start0_mod", 0)),
            "mod_end0": int(first.get("end0_mod", 0)),
            "target_mod_code": first.get("target_mod_code", ""),
            "gene_names": first.get("gene_names", first.get("gene_name_mod", "")),
            "metagene_indices": first.get("metagene_indices", first.get("metagene_index_mod", "")),
            "n_reads": int(tt.sum()),
            "n_ref_reads": int(allele_totals.get("ref", 0)),
            "n_alt_reads": int(allele_totals.get("alt", 0)),
            "n_modified": int(table[1].sum()),
            "n_not_target": int(table[0].sum()),
            "test_name": test_name,
            "stat_name": stat_name,
            "stat_value": stat_value,
            "p_value": p_value,
            "effect_abs_delta_mod_frac": binary_rate_delta(tt),
            "per_state_json": json.dumps({
                "ref_modified": int(table.loc["ref", 1]),
                "ref_not_target": int(table.loc["ref", 0]),
                "alt_modified": int(table.loc["alt", 1]),
                "alt_not_target": int(table.loc["alt", 0]),
            }, separators=(",", ":")),
        })

    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
        out = out.sort_values(["p_adj_bh", "effect_abs_delta_mod_frac"], ascending=[True, False]).reset_index(drop=True)
    else:
        out = pd.DataFrame(columns=[
            "snp_id", "mod_site_id", "chrom", "pos1", "mod_start0", "mod_end0", "target_mod_code",
            "gene_names", "metagene_indices", "n_reads", "n_ref_reads", "n_alt_reads", "n_modified", "n_not_target",
            "test_name", "stat_name", "stat_value", "p_value", "effect_abs_delta_mod_frac",
            "per_state_json", "p_adj_bh"
        ])

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
