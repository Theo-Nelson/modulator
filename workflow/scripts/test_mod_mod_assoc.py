#!/usr/bin/env python3
"""Roaring-bitmap engine for co-localized mod x mod dependency -- replaces test_mod_mod_assoc.py.

Per metagene, each mod site gets two roaring bitmaps over the metagene's reads: {reads modified} and
{reads not-modified}. A site-pair's 2x2 table is then 4 intersection_cardinality calls. The current
test already dedups per (read, site) via a recs dict, so boolean bitmaps reproduce it exactly (no
overlap double-count here). Sites are swept in genomic order with a distance window (break past
--max-distance). One chromosome at a time -> bounded memory. Derived-field math copied verbatim from
the production script so the statistics are identical by construction.
"""
import argparse
import json
import math
import os
import shutil
import sys
import tempfile
import numpy as np
import pandas as pd
from pyroaring import BitMap

from genotype_utils import (benjamini_hochberg, binary_rate_delta, cmh_stratified_test,
                            context_key_from_row, drop_unassigned_reads, mh_common_odds_ratio,
                            mh_stratified_effect, run_contingency_test, shard_tsv_by_chrom,
                            strata_effect_reversed, stratified_primary, tsv_header)

# Columns the current test_mod_mod_assoc uses (incl. strand, and gene_name/metagene_index for context_key).
_MOD_WANT = ["sample", "qname", "mod_site_id", "chrom", "start0", "strand", "target_mod_code",
             "state_detail", "gene_name", "gene_names", "metagene_index",
             "usable", "fail", "within_alignment"]

OUT_COLS = [
    "mod_site_id_a", "mod_site_id_b", "chrom", "start0_a", "start0_b", "distance_bp", "strand",
    "mod_code_a", "mod_code_b", "context_key", "gene_names",
    "n_reads", "n_a_modified", "n_a_unmodified", "n_b_modified",
    "n_both_modified", "n_a_only", "n_b_only", "n_neither",
    "odds_ratio", "log2_odds_ratio", "direction",
    "n_strata_informative", "strata_effect_reversed", "test_name", "stat_name", "stat_value", "p_value",
    "effect_abs_delta_mod_frac",
    "test_name_pooled", "stat_value_pooled", "p_value_pooled", "effect_abs_delta_mod_frac_pooled",
    "odds_ratio_pooled", "log2_odds_ratio_pooled", "direction_pooled",
    "exp_both_modified_pooled", "exp_neither_pooled",
    "concordant_frac_obs_pooled", "concordant_frac_exp_pooled", "concordance_log2_ratio_pooled",
    "jaccard_both_pooled", "per_state_json", "p_adj_bh",
]


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--molecule-mods", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--max-distance", type=int, default=1000)
    ap.add_argument("--min-pair-reads", type=int, default=8)
    ap.add_argument("--min-state-reads", type=int, default=4)
    ap.add_argument("--max-sites-per-read", type=int, default=200)
    ap.add_argument("--test", choices=["auto", "fisher", "chi2"], default="auto")
    ap.add_argument("--pseudocount", type=float, default=0.5)
    return ap.parse_args()


def _rk(df):
    return (df["sample"].astype(str) + "\x00" + df["qname"].astype(str))


def _pair_row(sid_a, sid_b, counts, site_meta, args, strata_counts=None):
    n_both, n_a_only, n_b_only, n_neither = counts
    n_reads = n_both + n_a_only + n_b_only + n_neither
    if n_reads < int(args.min_pair_reads):
        return None
    n_a_mod = n_both + n_a_only
    n_a_unmod = n_b_only + n_neither
    if n_a_mod < int(args.min_state_reads) or n_a_unmod < int(args.min_state_reads):
        return None
    tt = [[float(n_both), float(n_a_only)], [float(n_b_only), float(n_neither)]]
    chrom_a, s0_a, strand_a, code_a, ctx, genes = site_meta[sid_a]
    chrom_b, s0_b, _sb, code_b, _cb, _g = site_meta[sid_b]
    n_b_mod = n_both + n_b_only
    n_b_unmod = n_a_only + n_neither
    union_mod = n_both + n_a_only + n_b_only

    # ---- POOLED descriptive statistics (single 2x2 over all reads) -- kept ONLY as *_pooled ----
    # These are confounded by sample when reads are pooled (both sites' rates tracking sample fakes a
    # dependency), so every reader-facing statistic below has a stratified primary computed further down.
    pooled_name, pooled_sn, pooled_stat, pooled_p = run_contingency_test(
        tt, test=args.test, pseudocount=args.pseudocount)
    pooled_eff = binary_rate_delta(tt)
    exp_both = n_a_mod * n_b_mod / n_reads
    exp_neither = n_a_unmod * n_b_unmod / n_reads
    or_pooled = ((n_both + 0.5) * (n_neither + 0.5)) / ((n_a_only + 0.5) * (n_b_only + 0.5))
    log2_or_pooled = math.log2(or_pooled)
    conc_obs = (n_both + n_neither) / n_reads
    conc_exp = (n_a_mod * n_b_mod + n_a_unmod * n_b_unmod) / (n_reads * n_reads)
    conc_log2 = math.log2(conc_obs / conc_exp) if conc_obs > 0 and conc_exp > 0 else 0.0
    _dir = lambda l2: ("CONCORDANT" if l2 > 0.32 else "MUTUALLY_EXCLUSIVE" if l2 < -0.32 else "INDEPENDENT")
    direction_pooled = _dir(log2_or_pooled)

    # ---- PRIMARY statistics ----
    # strata_counts is passed only when >1 sample exists (a cheap short-circuit in _modmod_for_chrom),
    # but the CMH is adopted as primary ONLY when >=2 strata are INFORMATIVE (stratified_primary): a
    # single informative stratum -- one sample, OR many samples where only one survives filtering -- is
    # an uncorrected asymptotic chi2, anti-conservative at these cell counts ([[6,1],[1,6]] -> 0.010 vs
    # Fisher's exact 0.029), so it falls back to the exact pooled test.
    strata = [[[float(cb), float(ca)], [float(cbb), float(cn)]]
              for (cb, ca, cbb, cn) in strata_counts] if strata_counts else []
    cmh = cmh_stratified_test(strata) if strata else None
    test_name, stat_name, stat_value, p_value, n_strata, _used = stratified_primary(
        cmh, (pooled_name, pooled_sn, pooled_stat, pooled_p))
    if _used:
        eff = mh_stratified_effect(strata)                      # sample-adjusted rate difference
        or_mh = mh_common_odds_ratio(strata)                    # sample-adjusted (Mantel-Haenszel) OR
        if math.isfinite(or_mh) and or_mh > 0:
            odds_ratio, log2_or = or_mh, math.log2(or_mh)
        else:                                                   # no discordant mass -> fall back to pooled OR
            odds_ratio, log2_or = or_pooled, log2_or_pooled
        direction = _dir(log2_or)
    else:
        eff = pooled_eff
        odds_ratio, log2_or, direction = or_pooled, log2_or_pooled, direction_pooled
    # heterogeneity: does the co-modification effect REVERSE sign across samples? The CMH averages such
    # strata to ~0/p~1, so flag it (mirror of the pooling confound; only meaningful once stratified).
    strata_reversed = strata_effect_reversed(strata) if _used else False
    return {
        "odds_ratio": round(odds_ratio, 4), "log2_odds_ratio": round(log2_or, 4), "direction": direction,
        "exp_both_modified_pooled": round(exp_both, 2), "exp_neither_pooled": round(exp_neither, 2),
        "odds_ratio_pooled": round(or_pooled, 4), "log2_odds_ratio_pooled": round(log2_or_pooled, 4),
        "direction_pooled": direction_pooled,
        "concordant_frac_obs_pooled": round(conc_obs, 4), "concordant_frac_exp_pooled": round(conc_exp, 4),
        "concordance_log2_ratio_pooled": round(conc_log2, 4),
        "mod_site_id_a": sid_a, "mod_site_id_b": sid_b,
        "chrom": chrom_a, "start0_a": s0_a, "start0_b": s0_b,
        "distance_bp": abs(s0_b - s0_a), "strand": strand_a,
        "mod_code_a": code_a, "mod_code_b": code_b,
        "context_key": ctx, "gene_names": genes,
        "n_reads": n_reads,
        "n_a_modified": n_a_mod, "n_a_unmodified": n_a_unmod, "n_b_modified": n_b_mod,
        "n_both_modified": n_both, "n_a_only": n_a_only, "n_b_only": n_b_only, "n_neither": n_neither,
        "n_strata_informative": int(n_strata), "strata_effect_reversed": bool(strata_reversed),
        "test_name": test_name, "stat_name": stat_name, "stat_value": stat_value, "p_value": p_value,
        "effect_abs_delta_mod_frac": eff,
        "test_name_pooled": pooled_name, "stat_value_pooled": pooled_stat, "p_value_pooled": pooled_p,
        "effect_abs_delta_mod_frac_pooled": pooled_eff,
        "jaccard_both_pooled": round(n_both / union_mod, 6) if union_mod else 0.0,
        "per_state_json": json.dumps({
            "both_modified": n_both, "a_only": n_a_only,
            "b_only": n_b_only, "neither": n_neither,
        }, separators=(",", ":")),
    }


def _modmod_for_chrom(mod_path, args):
    hdr = tsv_header(mod_path)
    mod_df = pd.read_csv(mod_path, sep="\t", low_memory=False, usecols=[c for c in _MOD_WANT if c in hdr])
    if mod_df.empty:
        return [], 0
    # same usability + state filter as the production test
    if "usable" in mod_df.columns:
        mod_df = mod_df[mod_df["usable"].fillna(False)].copy()
    else:
        mod_df = mod_df[(~mod_df["fail"].fillna(True)) & mod_df["within_alignment"].fillna(False)].copy()
    mod_df = mod_df[mod_df["state_detail"].isin(["modified", "canonical", "other_mod"])].copy()
    mod_df = drop_unassigned_reads(mod_df)
    if mod_df.empty:
        return [], 0
    mod_df["target_state"] = mod_df["state_detail"].eq("modified").astype(int)
    mod_df["context_key"] = mod_df.apply(context_key_from_row, axis=1)
    gene_col = "gene_names" if "gene_names" in mod_df.columns else ("gene_name" if "gene_name" in mod_df.columns else None)

    rows = []
    n_skipped_dense = 0
    for ck, mm in mod_df.groupby("context_key", sort=False):
        mk = _rk(mm)
        ridx = {k: i for i, k in enumerate(pd.unique(mk))}
        mm = mm.assign(_ri=mk.map(ridx).to_numpy())
        # one record per (read, site) -- the recs-dict dedup of the production test
        mm = mm.drop_duplicates(["_ri", "mod_site_id"], keep="first")
        # drop reads carrying too many sites (the max_sites_per_read cap: those reads form no pairs)
        if args.max_sites_per_read:
            per_read = mm.groupby("_ri")["mod_site_id"].transform("size")
            dense = mm["_ri"][per_read > args.max_sites_per_read].nunique()
            n_skipped_dense += int(dense)
            mm = mm[per_read <= args.max_sites_per_read]
        if mm.empty:
            continue

        site_bits = {}
        site_meta = {}
        for site, g in mm.groupby("mod_site_id", sort=False):
            ri = g["_ri"].to_numpy(); ts = g["target_state"].to_numpy()
            site_bits[site] = (BitMap(ri[ts == 1].astype(np.uint32)),
                               BitMap(ri[ts == 0].astype(np.uint32)))
            f = g.iloc[0]
            site_meta[site] = (f["chrom"], int(f["start0"]), f["strand"], f["target_mod_code"],
                               ck, (f.get(gene_col, "") if gene_col else ""))

        # per-sample read bitmaps over the SAME _ri space, for sample-stratified CMH. Each _ri belongs
        # to exactly one sample (it is keyed by sample\x00qname). Sorted so strata order is stable.
        sample_bits = {}
        if "sample" in mm.columns:
            for samp, sg in mm.groupby("sample", sort=True):
                sample_bits[samp] = BitMap(pd.unique(sg["_ri"]).astype(np.uint32))

        # genomic-order sweep with a distance window (chrom constant within a metagene)
        order = sorted(site_meta, key=lambda s: (site_meta[s][0], site_meta[s][1], s))
        for i, sa in enumerate(order):
            s0a = site_meta[sa][1]; ma, ua = site_bits[sa]
            for j in range(i + 1, len(order)):
                sb = order[j]
                if site_meta[sb][0] != site_meta[sa][0]:
                    continue
                if site_meta[sb][1] - s0a > args.max_distance:
                    break
                mb, ub = site_bits[sb]
                b_both = ma & mb; b_ao = ma & ub; b_bo = ua & mb; b_ne = ua & ub
                counts = (len(b_both), len(b_ao), len(b_bo), len(b_ne))
                strata_counts = None
                if len(sample_bits) > 1:   # single sample -> pooled EXACT test (asymptotic CMH is anti-conservative)
                    strata_counts = [
                        (b_both.intersection_cardinality(sbm), b_ao.intersection_cardinality(sbm),
                         b_bo.intersection_cardinality(sbm), b_ne.intersection_cardinality(sbm))
                        for sbm in sample_bits.values()]
                r = _pair_row(sa, sb, counts, site_meta, args, strata_counts=strata_counts)
                if r is not None:
                    rows.append(r)
    return rows, n_skipped_dense


def main():
    args = parse_args()
    if not (os.path.exists(args.molecule_mods) and os.path.getsize(args.molecule_mods)):
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return
    tmp = tempfile.mkdtemp(prefix=".modmod_bs_", dir=os.path.dirname(args.out_tsv) or ".")
    rows = []
    n_skipped = 0
    try:
        mod_shards = shard_tsv_by_chrom(args.molecule_mods, os.path.join(tmp, "mod"))
        for chrom in sorted(mod_shards):
            r, nsk = _modmod_for_chrom(mod_shards[chrom], args)
            rows.extend(r); n_skipped += nsk
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if n_skipped:
        # Surface the cap's effect PERSISTENTLY, not just on stdout: dense reads are long reads that
        # carry most of the co-occurrence information, so a large skip count means the mod-mod result
        # under-represents them. Written to a sidecar next to the output so it survives the run.
        msg = (f"skipped pairing on {n_skipped} read(s) that exceeded --max-sites-per-read="
               f"{args.max_sites_per_read}; these are long, site-dense reads and carry disproportionate "
               f"co-occurrence signal -- raise --max-sites-per-read to include them.")
        print(f"[warn] {msg}", file=sys.stderr, flush=True)
        try:
            with open(args.out_tsv + ".skipped_dense_reads.txt", "w") as _fh:
                _fh.write(f"n_reads_skipped\t{n_skipped}\nmax_sites_per_read\t{args.max_sites_per_read}\n{msg}\n")
        except OSError:
            pass
    else:
        # No reads skipped this run -> remove any STALE sidecar from a previous run, or a resumed run
        # would keep reporting a skip count that no longer applies to the current output.
        try:
            os.remove(args.out_tsv + ".skipped_dense_reads.txt")
        except OSError:
            pass

    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
        out = out.sort_values(["p_adj_bh", "effect_abs_delta_mod_frac"],
                              ascending=[True, False]).reset_index(drop=True)
        out = out[OUT_COLS]
    else:
        out = pd.DataFrame(columns=OUT_COLS)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
