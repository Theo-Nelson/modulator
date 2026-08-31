#!/usr/bin/env python3
"""Value-asserting tests for the genotype_utils statistical primitives that ran in every pipeline but
had no test checking their RETURNED VALUE (only their existence). This is exactly the surface both
heterogeneity MAJORs hid in: a function can execute, return a plausible number, and be silently wrong.
Where a canonical reference exists the value is pinned to it (van Elteren single-stratum == scipy
Kruskal-Wallis H; benjamini_hochberg == scipy.false_discovery_control); otherwise to a hand-computed
expected value. Covers: benjamini_hochberg, max_abs_distribution_shift, stratified_max_distribution_shift,
binary_rate_delta, mh_common_odds_ratio, mh_stratified_effect, weighted_within_stratum_median_range,
van_elteren_kw, van_elteren_heterogeneity, add_heterogeneity_flag, stratified_primary."""
import os
import sys

import numpy as np
import pandas as pd
import scipy.stats as ss

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "workflow", "scripts"))
import genotype_utils as gu  # noqa: E402


def approx(a, b, tol=1e-6):
    return np.isfinite(a) and np.isfinite(b) and abs(float(a) - float(b)) <= tol


def main():
    checks = []

    # --- benjamini_hochberg: pinned to scipy.false_discovery_control; NaN passthrough ---
    pv = [0.001, 0.5, 0.04, 0.2]
    bh = gu.benjamini_hochberg(pv)
    ref = ss.false_discovery_control(pv, method="bh")
    checks.append(("benjamini_hochberg == scipy.false_discovery_control(bh)", np.allclose(bh, ref, atol=1e-9)))
    bh_nan = gu.benjamini_hochberg([0.001, np.nan, 0.04])
    checks.append(("benjamini_hochberg keeps NaN as NaN and ranks only finite (0.001->0.002)",
                   np.isnan(bh_nan[1]) and approx(bh_nan[0], 0.002) and approx(bh_nan[2], 0.04)))
    checks.append(("benjamini_hochberg clips to <=1 (all large p)",
                   float(np.max(gu.benjamini_hochberg([0.9, 0.95, 0.99]))) <= 1.0))

    # --- max_abs_distribution_shift: [[8,2],[2,8]] -> |0.8-0.2| = 0.6 ---
    checks.append(("max_abs_distribution_shift [[8,2],[2,8]] == 0.6", approx(gu.max_abs_distribution_shift([[8, 2], [2, 8]]), 0.6)))
    checks.append(("max_abs_distribution_shift degenerate (1 row) == 0.0", approx(gu.max_abs_distribution_shift([[5, 5]]), 0.0)))

    # --- stratified_max_distribution_shift: identical strata reduce to the single-table shift; within-0 -> 0 ---
    checks.append(("stratified_max_distribution_shift 2x[[8,2],[2,8]] == 0.6",
                   approx(gu.stratified_max_distribution_shift([[[8, 2], [2, 8]], [[8, 2], [2, 8]]]), 0.6)))
    checks.append(("stratified_max_distribution_shift within-stratum-equal rows == 0.0",
                   approx(gu.stratified_max_distribution_shift([[[5, 5], [5, 5]]]), 0.0)))
    # 3-row case where MAX over row pairs (0.6, rows 0 vs 2) != MEAN over pairs (0.4): pins the reducer
    # as a max, not a mean (a 2x2-only fixture cannot separate them -- there is a single row pair).
    checks.append(("stratified_max_distribution_shift 3-row == max over pairs (0.6, not mean 0.4)",
                   approx(gu.stratified_max_distribution_shift([[[80, 20], [50, 50], [20, 80]]]), 0.6)))

    # --- binary_rate_delta: rate0=8/10, rate1=3/10 -> 0.5; non-2x2 -> 0 ---
    checks.append(("binary_rate_delta [[8,2],[3,7]] == 0.5", approx(gu.binary_rate_delta([[8, 2], [3, 7]]), 0.5)))
    checks.append(("binary_rate_delta non-2x2 == 0.0", approx(gu.binary_rate_delta([[1, 2, 3], [4, 5, 6]]), 0.0)))

    # --- mh_common_odds_ratio: single [[8,2],[2,8]] -> 64/4 = 16; no discordant mass -> NaN ---
    checks.append(("mh_common_odds_ratio [[8,2],[2,8]] == 16.0", approx(gu.mh_common_odds_ratio([[[8, 2], [2, 8]]]), 16.0)))
    checks.append(("mh_common_odds_ratio mutually-exclusive -> NaN", np.isnan(gu.mh_common_odds_ratio([[[10, 0], [0, 10]]]))))
    checks.append(("mh_common_odds_ratio is stratum-invariant for identical strata (still 16)",
                   approx(gu.mh_common_odds_ratio([[[8, 2], [2, 8]], [[8, 2], [2, 8]]]), 16.0)))

    # --- mh_stratified_effect: single [[8,2],[2,8]] col-0 rate diff = 0.6 ---
    checks.append(("mh_stratified_effect [[8,2],[2,8]] == 0.6", approx(gu.mh_stratified_effect([[[8, 2], [2, 8]]]), 0.6)))

    # --- weighted_within_stratum_median_range: medians 2 vs 20 -> range 18; two-stratum weighted -> 11.5 ---
    checks.append(("weighted_within_stratum_median_range single == 18.0",
                   approx(gu.weighted_within_stratum_median_range([[[1, 2, 3], [10, 20, 30]]]), 18.0)))
    checks.append(("weighted_within_stratum_median_range 2 strata (6*18+6*5)/12 == 11.5",
                   approx(gu.weighted_within_stratum_median_range([[[1, 2, 3], [10, 20, 30]], [[0, 0, 0], [5, 5, 5]]]), 11.5)))

    # --- van_elteren_kw: single stratum == tie-corrected Kruskal-Wallis (scipy) ---
    g1, g2, g3 = [1., 2, 3, 4, 5], [3., 4, 5, 6, 7], [10., 11, 12, 13, 14]
    stat, p, used, df = gu.van_elteren_kw([[g1, g2, g3]])
    H, pk = ss.kruskal(g1, g2, g3)
    checks.append(("van_elteren_kw single stratum stat == scipy Kruskal-Wallis H", approx(stat, H) and approx(p, pk)))
    checks.append(("van_elteren_kw single stratum df == g-1 == 2", df == 2 and used == 1))
    # two strata with the SAME effect accumulate to a MORE significant p than one alone (added evidence)
    st2, p2, _, _ = gu.van_elteren_kw([[g1, g2, g3], [g1, g2, g3]])
    checks.append(("van_elteren_kw pools evidence: 2 identical strata -> smaller p than 1", p2 < p))

    # --- van_elteren_heterogeneity: homogeneous ordering -> not significant; reversed -> significant ---
    homo = [[g1, g3], [g1, g3]]
    rev = [[g1, g3], [g3, g1]]
    _, ph, _, _ = gu.van_elteren_heterogeneity(homo)
    _, pr, _, _ = gu.van_elteren_heterogeneity(rev)
    checks.append(("van_elteren_heterogeneity homogeneous strata -> p > 0.10", np.isfinite(ph) and ph > 0.10))
    checks.append(("van_elteren_heterogeneity reversed ordering -> p < 0.05", np.isfinite(pr) and pr < 0.05))

    # --- add_heterogeneity_flag: BH-adjust, True/False/<NA> semantics; NaN row -> <NA>, not False ---
    df_in = pd.DataFrame({"strata_heterogeneity_p": [0.0005, 0.5, np.nan, 0.02]})
    out = gu.add_heterogeneity_flag(df_in.copy(), alpha=0.05)
    padj = out["strata_heterogeneity_p_adj"].values
    het = out["strata_heterogeneous"]
    ref_adj = gu.benjamini_hochberg([0.0005, 0.5, np.nan, 0.02])
    checks.append(("add_heterogeneity_flag: p_adj == benjamini_hochberg", np.allclose(padj, ref_adj, atol=1e-9, equal_nan=True)))
    checks.append(("add_heterogeneity_flag: significant row (0.0005) -> True", het.iloc[0] is True or het.iloc[0] == True))  # noqa: E712
    checks.append(("add_heterogeneity_flag: homogeneous row (0.5) -> False", het.iloc[1] == False))  # noqa: E712
    checks.append(("add_heterogeneity_flag: untested (NaN p) row -> <NA>, NOT False", pd.isna(het.iloc[2])))
    checks.append(("add_heterogeneity_flag: strata_heterogeneous is nullable boolean", str(het.dtype) == "boolean"))

    # --- stratified_primary: mode routing (cmh / exact_single / none) ---
    def stub_exact(tbl):
        stub_exact.called_with = np.asarray(tbl)
        return ("stub_exact", "stub_stat", 1.23, 0.0456)
    two = [np.array([[8, 2], [2, 8]], float), np.array([[7, 3], [3, 7]], float)]
    r_cmh = gu.stratified_primary(two, stub_exact)
    checks.append(("stratified_primary >=2 strata -> mode 'cmh', n_informative==2", r_cmh[5] == "cmh" and r_cmh[4] == 2))
    one = [np.array([[6, 1], [1, 6]], float)]
    r_one = gu.stratified_primary(one, stub_exact)
    checks.append(("stratified_primary 1 stratum -> mode 'exact_single', uses the exact test (p==0.0456)",
                   r_one[5] == "exact_single" and approx(r_one[3], 0.0456) and r_one[0] == "stub_exact"))
    checks.append(("stratified_primary exact_single called on THAT stratum's table",
                   np.array_equal(getattr(stub_exact, "called_with", None), one[0])))
    r_none = gu.stratified_primary([], stub_exact)
    checks.append(("stratified_primary 0 strata -> mode 'none', NaN p", r_none[5] == "none" and np.isnan(r_none[3])))

    n_pass = 0
    for name, ok in checks:
        print(f"  [{'PASS' if ok else '**FAIL**'}] {name}")
        n_pass += bool(ok)
    print(f"\ngenotype_utils stat-primitive values: {n_pass}/{len(checks)} checks passed")
    sys.exit(0 if n_pass == len(checks) else 1)


if __name__ == "__main__":
    main()
