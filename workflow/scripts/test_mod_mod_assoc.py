#!/usr/bin/env python3
"""
test_mod_mod_assoc.py -- co-localized modification dependency (mod x mod).

The modification analogue of test_snp_mod_assoc.py: instead of asking "does a SNP allele change
the modification state at a nearby site", it asks "do two nearby modification sites co-occur on
the SAME molecule more (or less) often than expected".

For every pair of candidate mod sites (A, B) that
  * are seen on at least one shared read,
  * lie in the same context (same metagene/gene, via context_key), and
  * are within --max-distance bp of each other,
we build the per-read 2x2 table

                        B modified   B not-modified
     A modified            n11            n10
     A not-modified        n01            n00

and test independence (Fisher 2x2 / chi-square), BH-correcting across all tested pairs. The effect
size is the absolute difference in P(B modified) between A-modified and A-not-modified reads --
i.e. how much knowing A's state tells you about B.

IMPORTANT (memory): we deliberately do NOT self-merge the molecule table on (sample, qname). A
read carrying m mod calls would emit m^2 rows, so a global many-to-many merge blows up
superlinearly on dense loci. Instead we stream read-by-read and accumulate only the 2x2 counts for
pairs that pass the distance + context filters, so peak memory is O(#kept pairs), not O(sum m^2).
"""

import argparse
import json
import math
import os
from collections import defaultdict
from itertools import combinations

import pandas as pd

from genotype_utils import benjamini_hochberg, binary_rate_delta, context_key_from_row, run_contingency_test


OUT_COLS = [
    "mod_site_id_a", "mod_site_id_b", "chrom", "start0_a", "start0_b", "distance_bp", "strand",
    "mod_code_a", "mod_code_b", "context_key", "gene_names",
    "n_reads", "n_a_modified", "n_a_unmodified", "n_b_modified",
    "n_both_modified", "n_a_only", "n_b_only", "n_neither",
    # observed vs expected-under-independence for the concordant cells (both-modified / neither),
    # so you can read off "are they modified together / unmodified together more than chance?"
    "exp_both_modified", "exp_neither",
    "odds_ratio", "log2_odds_ratio",
    "concordant_frac_obs", "concordant_frac_exp", "concordance_log2_ratio", "direction",
    "test_name", "stat_name", "stat_value", "p_value",
    "effect_abs_delta_mod_frac", "jaccard_both", "per_state_json", "p_adj_bh",
]


def parse_args():
    ap = argparse.ArgumentParser(description="Test co-localized modification (mod x mod) dependency on shared molecules.")
    ap.add_argument("--molecule-mods", required=True, help="molecule_mod_calls.tsv (per-read mod calls)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--max-distance", type=int, default=1000,
                    help="Max genomic distance (bp) between two mod sites to call them co-localized. Default 1000.")
    ap.add_argument("--min-pair-reads", type=int, default=8,
                    help="Min shared reads covering BOTH sites to test a pair.")
    ap.add_argument("--min-state-reads", type=int, default=4,
                    help="Min reads in each state of site A (modified / not-modified).")
    ap.add_argument("--max-sites-per-read", type=int, default=200,
                    help="Safety cap: skip pairing for reads carrying more than this many usable mod calls "
                         "(pairs grow quadratically). 0 disables the cap.")
    ap.add_argument("--test", choices=["auto", "fisher", "chi2"], default="auto")
    ap.add_argument("--pseudocount", type=float, default=0.5)
    return ap.parse_args()


def _write_empty(path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    pd.DataFrame(columns=OUT_COLS).to_csv(path, sep="\t", index=False)


def main():
    args = parse_args()

    if not os.path.exists(args.molecule_mods) or os.path.getsize(args.molecule_mods) == 0:
        _write_empty(args.out_tsv)
        return

    mod_df = pd.read_csv(args.molecule_mods, sep="\t", low_memory=False)
    if mod_df.empty:
        _write_empty(args.out_tsv)
        return

    # Same usability filter as test_snp_mod_assoc.py so the two analyses see the same molecules.
    if "usable" in mod_df.columns:
        mod_df = mod_df[mod_df["usable"].fillna(False)].copy()
    else:
        mod_df = mod_df[(~mod_df["fail"].fillna(True)) & mod_df["within_alignment"].fillna(False)].copy()
    mod_df = mod_df[mod_df["state_detail"].isin(["modified", "canonical", "other_mod"])].copy()
    if mod_df.empty:
        _write_empty(args.out_tsv)
        return

    mod_df["target_state"] = mod_df["state_detail"].eq("modified").astype(int)
    mod_df["context_key"] = mod_df.apply(context_key_from_row, axis=1)

    gene_col = "gene_names" if "gene_names" in mod_df.columns else ("gene_name" if "gene_name" in mod_df.columns else None)

    # site_id -> static metadata (first occurrence wins; a site has one chrom/pos/strand/mod_code)
    site_meta = {}
    # (site_a, site_b) -> [n_both, n_a_only, n_b_only, n_neither]
    pair_counts = defaultdict(lambda: [0, 0, 0, 0])

    cols = ["sample", "qname", "mod_site_id", "chrom", "start0", "strand",
            "target_mod_code", "target_state", "context_key"]
    if gene_col:
        cols.append(gene_col)
    sub = mod_df[cols]

    n_skipped_dense = 0
    for (_sample, _qname), grp in sub.groupby(["sample", "qname"], sort=False):
        # Deduplicate: one call per (read, site).
        recs = {}
        for t in grp.itertuples(index=False):
            sid = t.mod_site_id
            if sid not in recs:
                recs[sid] = t
                if sid not in site_meta:
                    site_meta[sid] = (t.chrom, int(t.start0), t.strand, t.target_mod_code,
                                      t.context_key, (getattr(t, gene_col) if gene_col else ""))
        if len(recs) < 2:
            continue
        if args.max_sites_per_read and len(recs) > args.max_sites_per_read:
            n_skipped_dense += 1
            continue

        items = sorted(recs.values(), key=lambda t: (t.chrom, int(t.start0), t.mod_site_id))
        for ta, tb in combinations(items, 2):
            if ta.chrom != tb.chrom:
                continue
            if ta.context_key != tb.context_key:
                continue
            d = abs(int(tb.start0) - int(ta.start0))
            if d > args.max_distance:
                continue
            a_mod = int(ta.target_state)
            b_mod = int(tb.target_state)
            cell = 0 if (a_mod and b_mod) else 1 if (a_mod and not b_mod) else 2 if (b_mod and not a_mod) else 3
            pair_counts[(ta.mod_site_id, tb.mod_site_id)][cell] += 1

    if n_skipped_dense:
        print(f"[warn] skipped pairing on {n_skipped_dense} read(s) exceeding --max-sites-per-read="
              f"{args.max_sites_per_read}", flush=True)

    rows = []
    for (sid_a, sid_b), (n_both, n_a_only, n_b_only, n_neither) in pair_counts.items():
        n_reads = n_both + n_a_only + n_b_only + n_neither
        if n_reads < int(args.min_pair_reads):
            continue
        n_a_mod = n_both + n_a_only
        n_a_unmod = n_b_only + n_neither
        if n_a_mod < int(args.min_state_reads) or n_a_unmod < int(args.min_state_reads):
            continue

        # rows = A state (modified, not), cols = B state (modified, not)
        tt = [[float(n_both), float(n_a_only)], [float(n_b_only), float(n_neither)]]
        test_name, stat_name, stat_value, p_value = run_contingency_test(
            tt, test=args.test, pseudocount=args.pseudocount)

        chrom_a, s0_a, strand_a, code_a, ctx, genes = site_meta[sid_a]
        chrom_b, s0_b, _strand_b, code_b, _ctx_b, _g = site_meta[sid_b]
        n_b_mod = n_both + n_b_only
        n_b_unmod = n_a_only + n_neither
        union_mod = n_both + n_a_only + n_b_only

        # --- Concordance calibration: observed vs expected-under-independence ---------------
        # If the two sites were modified independently, the expected count in each 2x2 cell is
        # row_total * col_total / n. Co-regulation shows up as the CONCORDANT cells -- both-modified
        # (n_both) and both-unmodified (n_neither) -- being LARGER than expected (and the discordant
        # cells smaller). We report that directly:
        #   exp_both_modified / exp_neither : the independence expectation for the two concordant cells
        #   odds_ratio = (both * neither) / (a_only * b_only), 0.5-smoothed : >1 => concordant
        #   concordant_frac_obs vs _exp     : P(same state) observed vs under independence
        #   direction                       : CONCORDANT (co-modified) / MUTUALLY_EXCLUSIVE / INDEPENDENT
        exp_both = n_a_mod * n_b_mod / n_reads
        exp_neither = n_a_unmod * n_b_unmod / n_reads
        odds_ratio = ((n_both + 0.5) * (n_neither + 0.5)) / ((n_a_only + 0.5) * (n_b_only + 0.5))
        log2_or = math.log2(odds_ratio)
        conc_obs = (n_both + n_neither) / n_reads
        conc_exp = (n_a_mod * n_b_mod + n_a_unmod * n_b_unmod) / (n_reads * n_reads)
        conc_log2 = math.log2(conc_obs / conc_exp) if conc_obs > 0 and conc_exp > 0 else 0.0
        direction = ("CONCORDANT" if log2_or > 0.32 else            # OR > ~1.25
                     "MUTUALLY_EXCLUSIVE" if log2_or < -0.32 else "INDEPENDENT")

        rows.append({
            "exp_both_modified": round(exp_both, 2), "exp_neither": round(exp_neither, 2),
            "odds_ratio": round(odds_ratio, 4), "log2_odds_ratio": round(log2_or, 4),
            "concordant_frac_obs": round(conc_obs, 4), "concordant_frac_exp": round(conc_exp, 4),
            "concordance_log2_ratio": round(conc_log2, 4), "direction": direction,
            "mod_site_id_a": sid_a, "mod_site_id_b": sid_b,
            "chrom": chrom_a, "start0_a": s0_a, "start0_b": s0_b,
            "distance_bp": abs(s0_b - s0_a), "strand": strand_a,
            "mod_code_a": code_a, "mod_code_b": code_b,
            "context_key": ctx, "gene_names": genes,
            "n_reads": n_reads,
            "n_a_modified": n_a_mod, "n_a_unmodified": n_a_unmod, "n_b_modified": n_b_mod,
            "n_both_modified": n_both, "n_a_only": n_a_only, "n_b_only": n_b_only, "n_neither": n_neither,
            "test_name": test_name, "stat_name": stat_name, "stat_value": stat_value, "p_value": p_value,
            # |P(B mod | A mod) - P(B mod | A not-mod)|
            "effect_abs_delta_mod_frac": binary_rate_delta(tt),
            # co-modification Jaccard: reads with both / reads with either
            "jaccard_both": round(n_both / union_mod, 6) if union_mod else 0.0,
            "per_state_json": json.dumps({
                "both_modified": n_both, "a_only": n_a_only,
                "b_only": n_b_only, "neither": n_neither,
            }, separators=(",", ":")),
        })

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
    print(f"[ok] wrote {args.out_tsv}: {len(out)} co-localized mod-site pair(s) tested "
          f"(max_distance={args.max_distance})")


if __name__ == "__main__":
    main()
