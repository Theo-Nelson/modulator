#!/usr/bin/env python3
"""Bitset engine v2 for SNP x mod association -- corrected + memory-bounded.

Differences from the byte-identical PoC:
  * BOOLEAN membership (not integer multiplicity): a read is one molecule, so it contributes to a
    (SNP, site) cell at most once. This DEDUPS the duplicate (read, site) rows that overlapping-gene
    annotation injects into the mod table -- i.e. it FIXES the ~1% double-count bug that the
    many-to-many merge (and the byte-identical PoC) reproduce. So results differ from the current
    pipeline only on overlap-region pairs, and are the *correct* counts.
  * Processed one chromosome at a time (shard_tsv_by_chrom), so peak RAM is one chromosome's data
    and never the whole genome-wide table.

Validation target: running the CURRENT test on a deduped input == running this engine on the raw
input (both dedupe; identical logic otherwise).
"""
import argparse
import json
import os
import shutil
import tempfile
import numpy as np
from pyroaring import BitMap
import pandas as pd

from genotype_utils import (benjamini_hochberg, binary_rate_delta, cmh_stratified_test,
                            context_key_from_row, context_key_from_snp_row, context_keys_from_snp_row,
                            load_molecule_mods_for_pairing, mh_stratified_effect,
                            run_contingency_test, shard_tsv_by_chrom, tsv_header)

SNP_USECOLS = ["sample", "qname", "snp_id", "chrom", "pos1", "start0", "end0",
               "allele_class", "gene_names", "metagene_indices"]
OUT_COLS = [
    "snp_id", "mod_site_id", "chrom", "pos1", "mod_start0", "mod_end0", "target_mod_code",
    "gene_names", "metagene_indices", "n_reads", "n_ref_reads", "n_alt_reads", "n_modified",
    "n_not_target", "test_name", "stat_name", "stat_value", "p_value", "n_strata_informative",
    "effect_abs_delta_mod_frac", "test_name_pooled", "p_value_pooled", "effect_pooled",
    "per_state_json", "p_adj_bh",
]


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--molecule-snps", required=True)
    ap.add_argument("--molecule-mods", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-allele-reads", type=int, default=4)
    ap.add_argument("--min-total-reads", type=int, default=8)
    ap.add_argument("--test", choices=["auto", "fisher", "chi2"], default="auto")
    ap.add_argument("--pseudocount", type=float, default=0.5)
    return ap.parse_args()


def _rk(df):
    return (df["sample"].astype(str) + "\x00" + df["qname"].astype(str))


def _pairs_for_one_chrom(mod_path, snp_path, args):
    """Boolean-bitset SNP x mod pairing for ONE chromosome's shard. Returns list of row dicts.
    Loads only this chromosome, groups by metagene, and for each metagene builds boolean membership
    vectors over its reads -- a cell count is then a vector AND + popcount (dedup by construction)."""
    mod_df = load_molecule_mods_for_pairing(mod_path)
    if mod_df.empty:
        return []
    snp_uc = [c for c in SNP_USECOLS if c in tsv_header(snp_path)]
    snp_df = pd.read_csv(snp_path, sep="\t", usecols=snp_uc, low_memory=False)
    snp_df = snp_df[snp_df["allele_class"].isin(["ref", "alt"])]
    if snp_df.empty:
        return []
    mod_df["context_key"] = mod_df.apply(context_key_from_row, axis=1)
    # Fan out each SNP over EVERY context (metagene) it spans, so a SNP overlapping >1 gene is paired
    # with modifications in each -- collapsing to one CHR: key silently drops multi-metagene SNPs
    # (0/123 reached snp_mod_assoc on real data before this).
    snp_df = snp_df.copy()
    snp_df["context_key"] = snp_df.apply(context_keys_from_snp_row, axis=1)
    snp_df = snp_df.explode("context_key", ignore_index=True)
    snp_by_ctx = {k: v for k, v in snp_df.groupby("context_key", sort=False)}

    rows = []
    for ck, mod_meta in mod_df.groupby("context_key", sort=False):
        snp_meta = snp_by_ctx.get(ck)
        if snp_meta is None:
            continue
        mk = _rk(mod_meta)
        ridx = {k: i for i, k in enumerate(pd.unique(mk))}
        R = len(ridx)

        mm = mod_meta.assign(_ri=mk.map(ridx).to_numpy())
        # per-SAMPLE read-id bitmaps, so the SNP x mod contingency table can be stratified by sample
        # (BLOCKER: a sample-pooled Fisher/chi2 manufactures a Simpson's-paradox association from
        # per-sample allele-composition + baseline-rate imbalance -- the same confound the stoichiometry
        # and tail tests were fixed for). read-id universe = the mod reads (ridx is built from them).
        samp_bits = {}
        if "sample" in mm.columns:
            for _samp, _sg in mm.groupby("sample", sort=False):
                samp_bits[str(_samp)] = BitMap(_sg["_ri"].to_numpy().astype(np.uint32))
        site_bits = {}
        for site, g in mm.groupby("mod_site_id", sort=False):
            ri = g["_ri"].to_numpy(); ts = g["target_state"].to_numpy()
            ma = BitMap(ri[ts == 1].astype(np.uint32))          # roaring set of read-ids (dedups)
            ua = BitMap(ri[ts == 0].astype(np.uint32))
            f = g.iloc[0]
            site_bits[site] = (ma, ua, int(f["start0"]), int(f.get("end0", f["start0"] + 1)),
                               f.get("target_mod_code", ""))

        sk = _rk(snp_meta)
        sm = snp_meta.assign(_ri=sk.map(ridx).to_numpy())
        sm = sm[sm["_ri"].notna()]
        if sm.empty:
            continue
        sm = sm.assign(_ri=sm["_ri"].astype(int))
        snp_bits = {}
        for s, g in sm.groupby("snp_id", sort=False):
            ri = g["_ri"].to_numpy(); ac = g["allele_class"].to_numpy()
            ra = BitMap(ri[ac == "ref"].astype(np.uint32))
            aa = BitMap(ri[ac == "alt"].astype(np.uint32))
            f = g.iloc[0]
            snp_bits[s] = (ra, aa, f.get("chrom", ""), int(f.get("pos1", 0)),
                           f.get("gene_names", ""), f.get("metagene_indices", ""))

        for s, (ra, aa, chrom, pos1, genes, metas) in snp_bits.items():
            for m, (ma, ua, s0, e0, code) in site_bits.items():
                ref_mod = ra.intersection_cardinality(ma)      # |A ∩ B| via SIMD, no materialisation
                ref_unmod = ra.intersection_cardinality(ua)
                alt_mod = aa.intersection_cardinality(ma)
                alt_unmod = aa.intersection_cardinality(ua)
                n_ref = ref_mod + ref_unmod
                n_alt = alt_mod + alt_unmod
                if n_ref < int(args.min_allele_reads) or n_alt < int(args.min_allele_reads):
                    continue
                if (n_ref + n_alt) < int(args.min_total_reads):
                    continue
                tt = np.array([[ref_mod, ref_unmod], [alt_mod, alt_unmod]], dtype=float)
                # pooled test (kept as *_pooled companion -- the confounded statistic)
                p_test_name, p_stat_name, p_stat_value, p_pooled = run_contingency_test(
                    tt, test=args.test, pseudocount=args.pseudocount)
                # PRIMARY: sample-stratified CMH. Build one 2x2 [[ref_mod,ref_unmod],[alt_mod,alt_unmod]]
                # per sample and combine; falls back to the pooled table when there is no sample column.
                if samp_bits:
                    strata = []
                    for _sb in samp_bits.values():
                        rm = len(ra & ma & _sb); ru = len(ra & ua & _sb)
                        am = len(aa & ma & _sb); au = len(aa & ua & _sb)
                        strata.append([[rm, ru], [am, au]])
                    test_name, stat_name, stat_value, p_value, n_strata = cmh_stratified_test(strata)
                    effect = mh_stratified_effect(strata)
                else:
                    test_name, stat_name, stat_value, p_value, n_strata = (
                        p_test_name, p_stat_name, p_stat_value, p_pooled, 0)
                    effect = binary_rate_delta(tt)
                rows.append({
                    "snp_id": s, "mod_site_id": m, "chrom": chrom, "pos1": pos1,
                    "mod_start0": s0, "mod_end0": e0, "target_mod_code": code,
                    "gene_names": genes, "metagene_indices": metas,
                    "n_reads": int(tt.sum()), "n_ref_reads": n_ref, "n_alt_reads": n_alt,
                    "n_modified": ref_mod + alt_mod, "n_not_target": ref_unmod + alt_unmod,
                    "test_name": test_name, "stat_name": stat_name, "stat_value": stat_value,
                    "p_value": p_value, "n_strata_informative": int(n_strata),
                    "effect_abs_delta_mod_frac": effect,
                    "test_name_pooled": p_test_name, "p_value_pooled": p_pooled,
                    "effect_pooled": binary_rate_delta(tt),
                    "per_state_json": json.dumps({
                        "ref_modified": ref_mod, "ref_not_target": ref_unmod,
                        "alt_modified": alt_mod, "alt_not_target": alt_unmod,
                    }, separators=(",", ":")),
                })
    return rows


def main():
    args = parse_args()
    for p in (args.molecule_mods, args.molecule_snps):
        if not (os.path.exists(p) and os.path.getsize(p)):
            pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
            return
    tmp = tempfile.mkdtemp(prefix=".bitset_shards_", dir=os.path.dirname(args.out_tsv) or ".")
    rows = []
    try:
        mod_shards = shard_tsv_by_chrom(args.molecule_mods, os.path.join(tmp, "mod"))
        snp_shards = shard_tsv_by_chrom(args.molecule_snps, os.path.join(tmp, "snp"))
        for chrom in sorted(set(mod_shards) & set(snp_shards)):
            rows.extend(_pairs_for_one_chrom(mod_shards[chrom], snp_shards[chrom], args))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
        out = out.sort_values(["p_adj_bh", "effect_abs_delta_mod_frac"],
                              ascending=[True, False]).reset_index(drop=True)
    else:
        out = pd.DataFrame(columns=OUT_COLS)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
