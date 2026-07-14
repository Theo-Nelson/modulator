#!/usr/bin/env python3
"""Roaring-bitmap engine for SNP x mod dependency stratified by transcript (CMH).

snp_mod + a transcript axis: per metagene, each SNP allele, mod state, and transcript (ZT) is a
roaring bitmap over the metagene's reads. A per-transcript 2x2 cell is a 3-way intersection
|ref & mod & tx|; CMH runs across the transcript strata. Boolean bitmaps dedupe the overlapping-gene
duplicates the merge double-counts (so this is the corrected result). Per chromosome -> bounded RAM.
Derived fields + CMH copied from the production script.
"""
import argparse
import json
import os
import shutil
import tempfile
import numpy as np
import pandas as pd
from pyroaring import BitMap

from genotype_utils import (benjamini_hochberg, binary_rate_delta, cmh_test_2x2xk,
                            context_key_from_row, context_key_from_snp_row,
                            load_molecule_mods_for_pairing, read_keys_of, shard_tsv_by_chrom,
                            stream_filter_by_read_keys, tsv_header)

SNP_USECOLS = ["sample", "qname", "snp_id", "chrom", "pos1", "start0", "end0",
               "allele_class", "gene_names", "metagene_indices", "ZT"]
EMPTY_COLS_SHORT = [
    "snp_id", "mod_site_id", "n_reads", "n_transcripts_tested", "cmh_stat", "cmh_p_value",
    "cmh_common_odds_ratio", "overall_effect_abs_delta_mod_frac", "weighted_within_tx_effect",
    "classification", "cmh_p_adj_bh",
]
EMPTY_COLS_FULL = [
    "snp_id", "mod_site_id", "chrom", "pos1", "mod_start0", "target_mod_code",
    "n_reads", "n_transcripts_tested", "cmh_stat", "cmh_p_value", "cmh_common_odds_ratio",
    "overall_effect_abs_delta_mod_frac", "weighted_within_tx_effect",
    "snp_tx_p_value", "snp_tx_p_adj_bh", "snp_tx_effect",
    "snp_mod_p_value", "snp_mod_p_adj_bh", "snp_mod_effect",
    "strata_json", "cmh_p_adj_bh", "classification",
]


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--molecule-snps", required=True)
    ap.add_argument("--molecule-mods", required=True)
    ap.add_argument("--snp-transcript-assoc", required=True)
    ap.add_argument("--snp-mod-assoc", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-stratum-reads", type=int, default=4)
    return ap.parse_args()


def _write(out, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    out.to_csv(path, sep="\t", index=False)


def _rk(df):
    return (df["sample"].astype(str) + "\x00" + df["qname"].astype(str))


def _strata_for_one_chrom(mod_path, snp_path, args, candidate_pairs, cand_snp_ids, cand_mod_ids,
                          snp_tx_map, snp_mod_map):
    mod_df = load_molecule_mods_for_pairing(mod_path, extra_cols=["ZT"])
    if not mod_df.empty:
        mod_df = mod_df[mod_df["mod_site_id"].isin(cand_mod_ids)].copy()
    mod_keys = read_keys_of(mod_df)
    snp_uc = [c for c in SNP_USECOLS if c in tsv_header(snp_path)]
    snp_df = stream_filter_by_read_keys(
        snp_path, snp_uc, mod_keys,
        row_filter=lambda ch: ch["allele_class"].isin(["ref", "alt"]) & ch["snp_id"].isin(cand_snp_ids),
    )
    if snp_df.empty or mod_df.empty:
        return [], False

    mod_df["context_key"] = mod_df.apply(context_key_from_row, axis=1)
    snp_df["context_key"] = snp_df.apply(context_key_from_snp_row, axis=1)
    snp_by_ctx = {k: v for k, v in snp_df.groupby("context_key", sort=False)}

    rows = []
    for ck, mod_meta in mod_df.groupby("context_key", sort=False):
        snp_meta = snp_by_ctx.get(ck)
        if snp_meta is None:
            continue
        mk = _rk(mod_meta)
        ridx = {k: i for i, k in enumerate(pd.unique(mk))}
        mm = mod_meta.assign(_ri=mk.map(ridx).to_numpy())

        site_mod = {}
        for site, g in mm.groupby("mod_site_id", sort=False):
            ri = g["_ri"].to_numpy(); ts = g["target_state"].to_numpy()
            f = g.iloc[0]
            site_mod[site] = (BitMap(ri[ts == 1].astype(np.uint32)), BitMap(ri[ts == 0].astype(np.uint32)),
                              int(f["start0"]), f.get("target_mod_code", ""))
        # transcript (ZT) bitmaps: one ZT per read
        rz = mm.drop_duplicates("_ri")
        tx_bm = {}
        for zt, g in rz.groupby(rz["ZT"].astype(str), sort=False):
            if zt == "" or zt == "nan":
                continue
            tx_bm[zt] = BitMap(g["_ri"].to_numpy().astype(np.uint32))
        tx_order = sorted(tx_bm)

        sk = _rk(snp_meta)
        sm = snp_meta.assign(_ri=sk.map(ridx).to_numpy())
        sm = sm[sm["_ri"].notna()]
        if sm.empty:
            continue
        sm = sm.assign(_ri=sm["_ri"].astype(int))
        snp_bm = {}
        for s, g in sm.groupby("snp_id", sort=False):
            ri = g["_ri"].to_numpy(); ac = g["allele_class"].to_numpy()
            f = g.iloc[0]
            snp_bm[s] = (BitMap(ri[ac == "ref"].astype(np.uint32)), BitMap(ri[ac == "alt"].astype(np.uint32)),
                         f.get("chrom", ""), int(f.get("pos1", 0)))

        for s, (ref_bm, alt_bm, chrom, pos1) in snp_bm.items():
            for m, (mod_bm, unmod_bm, s0, code) in site_mod.items():
                if (s, m) not in candidate_pairs:
                    continue
                rm = ref_bm & mod_bm; ru = ref_bm & unmod_bm
                am = alt_bm & mod_bm; au = alt_bm & unmod_bm
                strata = []
                stratum_details = []
                for tx in tx_order:
                    tb = tx_bm[tx]
                    a = rm.intersection_cardinality(tb); b = ru.intersection_cardinality(tb)
                    c = am.intersection_cardinality(tb); d = au.intersection_cardinality(tb)
                    total = a + b + c + d
                    if total < int(args.min_stratum_reads):
                        continue
                    if (a + b) == 0 or (c + d) == 0:
                        continue
                    table = [[a, b], [c, d]]
                    strata.append(table)
                    stratum_details.append({
                        "ZT": tx, "ref_modified": a, "ref_not_target": b,
                        "alt_modified": c, "alt_not_target": d,
                        "effect": binary_rate_delta(table), "n_reads": total,
                    })
                if not strata:
                    continue
                overall = [[len(rm), len(ru)], [len(am), len(au)]]
                cmh_stat, cmh_p_value, cmh_or = cmh_test_2x2xk(strata)
                weighted_effect = 0.0
                total_weight = 0
                for item in stratum_details:
                    weighted_effect += item["effect"] * item["n_reads"]
                    total_weight += item["n_reads"]
                weighted_effect = round(float(weighted_effect / total_weight), 6) if total_weight > 0 else 0.0
                snp_tx_meta = snp_tx_map.get(s, {})
                snp_mod_meta = snp_mod_map.get((s, m), {})
                rows.append({
                    "snp_id": s, "mod_site_id": m, "chrom": chrom, "pos1": pos1,
                    "mod_start0": s0, "target_mod_code": code,
                    "n_reads": int(sum(sum(sum(r) for r in t) for t in strata)),
                    "n_transcripts_tested": int(len(strata)),
                    "cmh_stat": cmh_stat, "cmh_p_value": cmh_p_value, "cmh_common_odds_ratio": cmh_or,
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
    return rows, True


def main():
    args = parse_args()
    snp_tx = pd.read_csv(args.snp_transcript_assoc, sep="\t", low_memory=False) if os.path.exists(args.snp_transcript_assoc) and os.path.getsize(args.snp_transcript_assoc) else pd.DataFrame()
    snp_mod = pd.read_csv(args.snp_mod_assoc, sep="\t", low_memory=False) if os.path.exists(args.snp_mod_assoc) and os.path.getsize(args.snp_mod_assoc) else pd.DataFrame()
    if snp_mod.empty:
        _write(pd.DataFrame(columns=EMPTY_COLS_SHORT), args.out_tsv)
        return
    candidate_pairs = set(zip(snp_mod["snp_id"], snp_mod["mod_site_id"]))
    cand_snp_ids = {s for s, _ in candidate_pairs}
    cand_mod_ids = {m for _, m in candidate_pairs}
    snp_tx_map = {}
    if not snp_tx.empty and "snp_id" in snp_tx.columns:
        snp_tx_map = snp_tx.set_index("snp_id")[["p_value", "p_adj_bh", "effect_max_abs_tx_frac_diff"]].to_dict("index")
    snp_mod_map = snp_mod.set_index(["snp_id", "mod_site_id"])[["p_value", "p_adj_bh", "effect_abs_delta_mod_frac"]].to_dict("index")

    if not (os.path.exists(args.molecule_mods) and os.path.getsize(args.molecule_mods)
            and os.path.exists(args.molecule_snps) and os.path.getsize(args.molecule_snps)):
        _write(pd.DataFrame(columns=EMPTY_COLS_SHORT), args.out_tsv)
        return

    tmp = tempfile.mkdtemp(prefix=".snptx_bs_", dir=os.path.dirname(args.out_tsv) or ".")
    rows = []
    had_any = False
    try:
        mod_shards = shard_tsv_by_chrom(args.molecule_mods, os.path.join(tmp, "mod"))
        snp_shards = shard_tsv_by_chrom(args.molecule_snps, os.path.join(tmp, "snp"))
        for chrom in sorted(set(mod_shards) & set(snp_shards)):
            r, had = _strata_for_one_chrom(mod_shards[chrom], snp_shards[chrom], args,
                                           candidate_pairs, cand_snp_ids, cand_mod_ids,
                                           snp_tx_map, snp_mod_map)
            rows.extend(r); had_any = had_any or had
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if not had_any:
        _write(pd.DataFrame(columns=EMPTY_COLS_SHORT), args.out_tsv)
        return

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
        out = out.sort_values(["cmh_p_adj_bh", "weighted_within_tx_effect", "snp_id", "mod_site_id"],
                              ascending=[True, False, True, True]).reset_index(drop=True)
    else:
        out = pd.DataFrame(columns=EMPTY_COLS_FULL)
    _write(out, args.out_tsv)


if __name__ == "__main__":
    main()
