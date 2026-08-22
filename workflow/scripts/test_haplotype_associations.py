#!/usr/bin/env python3
"""Haplotype associations -- roaring hap x mod, pandas hap x transcript.

hap x transcript touches only the (small) haplotype table, so it stays pandas (byte-identical, no
roaring gain). hap x mod merges the big mod table -> per metagene, each (block, haplotype) and each
mod site is a roaring bitmap; a (block, mod_site) contingency cell is |block_hap & mod_state|.
Boolean bitmaps dedupe the overlapping-gene duplicates the merge double-counts (corrected result).
Per chromosome -> bounded RAM.
"""
import argparse
import json
import os
import shutil
import tempfile
import numpy as np
import pandas as pd
from pyroaring import BitMap

from genotype_utils import (benjamini_hochberg, context_key_from_row, drop_unassigned_reads,
                            max_abs_distribution_shift, run_contingency_test, shard_tsv_by_chrom, tsv_header)

TX_COLS = [
    "block_id", "context_key", "chrom", "n_reads", "n_haplotypes_tested", "n_transcripts_tested",
    "test_name", "stat_name", "stat_value", "p_value", "effect_max_abs_tx_frac_diff", "per_table_json", "p_adj_bh",
]
MOD_COLS = [
    "block_id", "mod_site_id", "context_key", "chrom", "target_mod_code", "n_reads",
    "n_haplotypes_tested", "test_name", "stat_name", "stat_value", "p_value",
    "effect_max_abs_mod_rate_diff", "per_table_json", "p_adj_bh",
]
_MOD_WANT = ["sample", "qname", "mod_site_id", "chrom", "state_detail", "target_mod_code",
             "gene_name", "metagene_index", "usable", "fail", "within_alignment"]


def parse_args():
    ap = argparse.ArgumentParser()
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


def _rk(df):
    return (df["sample"].astype(str) + "\x00" + df["qname"].astype(str))


def _hap_for_one_chrom(hap_path, mod_path, args):
    hap = pd.read_csv(hap_path, sep="\t", low_memory=False)
    _hl = hap["haplotype"].fillna("").astype(str)
    # drop empty AND the pooled "OTHER" bucket: build_haplotype_blocks merges all sub-threshold
    # haplotypes into one "OTHER" label, which is not a real allele string and must not be tested as a
    # haplotype (it would enter the contingency table as an uninterpretable mixed row).
    hap = hap[_hl.ne("") & _hl.ne("OTHER")].copy()

    # ---- hap x transcript: pandas (unchanged, small) ----
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
            "block_id": block_id, "context_key": sub.iloc[0].get("context_key", ""),
            "chrom": sub.iloc[0].get("chrom", ""), "n_reads": int(tt.sum()),
            "n_haplotypes_tested": int(tt.shape[0]), "n_transcripts_tested": int(tt.shape[1]),
            "test_name": test_name, "stat_name": stat_name, "stat_value": stat_value, "p_value": p_value,
            "effect_max_abs_tx_frac_diff": max_abs_distribution_shift(tt),
            "per_table_json": json.dumps({"haplotypes": list(table.index), "transcripts": list(table.columns),
                                          "counts": table.values.tolist()}, separators=(",", ":")),
        })

    # ---- hap x mod: roaring ----
    mod_rows = []
    if mod_path is not None:
        _hdr = tsv_header(mod_path)
        mod = pd.read_csv(mod_path, sep="\t", low_memory=False, usecols=[c for c in _MOD_WANT if c in _hdr])
        if "usable" in mod.columns:
            mod = mod[mod["usable"].fillna(False)].copy()
        else:
            mod = mod[(~mod["fail"].fillna(True)) & mod["within_alignment"].fillna(False)].copy()
        mod = mod[mod["state_detail"].isin(["modified", "canonical", "other_mod"])].copy()
        mod = drop_unassigned_reads(mod)
        if not mod.empty:
            mod["target_state"] = mod["state_detail"].eq("modified").astype(int)
            mod["context_key"] = mod.apply(context_key_from_row, axis=1)
            hap_by_ctx = {k: v for k, v in hap.groupby("context_key", sort=False)} if "context_key" in hap.columns else {}
            for ck, mm in mod.groupby("context_key", sort=False):
                hm = hap_by_ctx.get(ck)
                if hm is None:
                    continue
                mk = _rk(mm)
                ridx = {k: i for i, k in enumerate(pd.unique(mk))}
                mm = mm.assign(_ri=mk.map(ridx).to_numpy())
                site_bits = {}
                for site, g in mm.groupby("mod_site_id", sort=False):
                    ri = g["_ri"].to_numpy(); ts = g["target_state"].to_numpy()
                    site_bits[site] = (BitMap(ri[ts == 1].astype(np.uint32)), BitMap(ri[ts == 0].astype(np.uint32)),
                                       g.iloc[0].get("target_mod_code", ""))
                hk = _rk(hm)
                hm = hm.assign(_ri=hk.map(ridx).to_numpy())
                hm = hm[hm["_ri"].notna()]
                if hm.empty:
                    continue
                hm = hm.assign(_ri=hm["_ri"].astype(int))
                block_hap = {}   # block -> {haplotype -> BitMap}
                block_meta = {}  # block -> (context_key_hap, chrom_hap)
                for (block, h), g in hm.groupby(["block_id", "haplotype"], sort=False):
                    block_hap.setdefault(block, {})[h] = BitMap(g["_ri"].to_numpy().astype(np.uint32))
                    if block not in block_meta:
                        f = g.iloc[0]
                        block_meta[block] = (f.get("context_key", ""), f.get("chrom", ""))
                for block, hbm in block_hap.items():
                    ck_hap, chrom_hap = block_meta[block]
                    for site, (mod_bm, unmod_bm, code) in site_bits.items():
                        cover = mod_bm | unmod_bm
                        keep = sorted(h for h, bm in hbm.items() if bm.intersection_cardinality(cover) >= int(args.min_haplotype_reads))
                        if len(keep) < 2:
                            continue
                        counts = [[hbm[h].intersection_cardinality(mod_bm), hbm[h].intersection_cardinality(unmod_bm)] for h in keep]
                        tt = np.array(counts, dtype=float)   # rows=haplotypes(sorted), cols=[modified, not]
                        if int(tt.sum()) < int(args.min_total_reads):
                            continue
                        if tt[:, 0].sum() == 0:
                            continue
                        test_name, stat_name, stat_value, p_value = run_contingency_test(tt, test=args.test, pseudocount=args.pseudocount)
                        mod_rows.append({
                            "block_id": block, "mod_site_id": site, "context_key": ck_hap, "chrom": chrom_hap,
                            "target_mod_code": code, "n_reads": int(tt.sum()), "n_haplotypes_tested": int(tt.shape[0]),
                            "test_name": test_name, "stat_name": stat_name, "stat_value": stat_value, "p_value": p_value,
                            "effect_max_abs_mod_rate_diff": max_abs_distribution_shift(tt),
                            "per_table_json": json.dumps({"haplotypes": keep, "states": ["modified", "not_target"],
                                                          "counts": tt.astype(int).tolist()}, separators=(",", ":")),
                        })
    return tx_rows, mod_rows


def main():
    args = parse_args()
    if not (os.path.exists(args.molecule_haplotypes) and os.path.getsize(args.molecule_haplotypes)):
        pd.DataFrame(columns=TX_COLS).to_csv(args.out_haplotype_transcript, sep="\t", index=False)
        pd.DataFrame(columns=MOD_COLS).to_csv(args.out_haplotype_mod, sep="\t", index=False)
        return
    tmp = tempfile.mkdtemp(prefix=".hap_bs_", dir=os.path.dirname(args.out_haplotype_mod) or ".")
    tx_rows, mod_rows = [], []
    try:
        hap_shards = shard_tsv_by_chrom(args.molecule_haplotypes, os.path.join(tmp, "hap"))
        mod_shards = ({} if not (os.path.exists(args.molecule_mods) and os.path.getsize(args.molecule_mods))
                      else shard_tsv_by_chrom(args.molecule_mods, os.path.join(tmp, "mod")))
        for chrom in sorted(hap_shards):
            tr, mr = _hap_for_one_chrom(hap_shards[chrom], mod_shards.get(chrom), args)
            tx_rows.extend(tr); mod_rows.extend(mr)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    tx_out = pd.DataFrame(tx_rows)
    if not tx_out.empty:
        tx_out["p_adj_bh"] = benjamini_hochberg(tx_out["p_value"].values)
        tx_out = tx_out.sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False]).reset_index(drop=True)
    else:
        tx_out = pd.DataFrame(columns=TX_COLS)
    mod_out = pd.DataFrame(mod_rows)
    if not mod_out.empty:
        mod_out["p_adj_bh"] = benjamini_hochberg(mod_out["p_value"].values)
        mod_out = mod_out.sort_values(["p_adj_bh", "effect_max_abs_mod_rate_diff", "block_id", "mod_site_id"],
                                      ascending=[True, False, True, True]).reset_index(drop=True)
    else:
        mod_out = pd.DataFrame(columns=MOD_COLS)
    os.makedirs(os.path.dirname(args.out_haplotype_transcript) or ".", exist_ok=True)
    tx_out.to_csv(args.out_haplotype_transcript, sep="\t", index=False)
    mod_out.to_csv(args.out_haplotype_mod, sep="\t", index=False)


if __name__ == "__main__":
    main()
