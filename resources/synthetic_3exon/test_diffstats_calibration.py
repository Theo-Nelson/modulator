#!/usr/bin/env python3
"""Null-calibration regression test for the beta-binomial between-condition test (guards C1).

The dispersion is estimated with a Cox-Reid adjusted profile likelihood (diffstats._fit_theta). Plain
ML biased the dispersion low and made the test ~2x anti-conservative on overdispersed replicates (null
type-I ~0.11 at phi>=0.005). This test simulates TRUE-NULL sites (both groups share one mu) at a few
dispersions and asserts the false-positive rate stays controlled -- it fails loudly if the bias returns.

Run:  <modulator-env>/bin/python resources/synthetic_3exon/test_diffstats_calibration.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "workflow", "scripts"))
import diffstats  # noqa: E402
import genotype_utils as gu  # noqa: E402


def mc_exact_rx2_null_type_i(depths, rate, seed0=0, M=800):
    """Type-I of the shared run_contingency_test single-stratum path on sparse r x 2 tables -- the M2
    regime that the asymptotic chi2 inflated to 0.15-0.25 and every single-sample site hits."""
    rng = np.random.default_rng(seed0)
    hit = used = 0
    for _ in range(M):
        tab = np.array([[(m := rng.binomial(d, rate)), d - m] for d in depths], dtype=float)
        if not ((tab.sum(0) > 0).all() and (tab.sum(1) > 0).all()):
            continue
        used += 1
        _, _, _, p = gu.run_contingency_test(tab)
        hit += (p < 0.05)
    return hit / max(1, used)


def null_type_i(phi, nrep, nsites=500, seed=20260822, mu_range=(0.2, 0.8), cov=(50, 300)):
    rng = np.random.default_rng(seed)
    theta = (1.0 / phi - 1.0) if phi > 0 else np.inf
    gidx = np.array([0] * nrep + [1] * nrep, dtype=int)
    sites = []
    for i in range(nsites):
        mu = rng.uniform(*mu_range)
        n = rng.integers(cov[0], cov[1], size=2 * nrep).astype(float)
        p = np.full(2 * nrep, mu) if np.isinf(theta) else rng.beta(mu * theta, (1 - mu) * theta, 2 * nrep)
        k = rng.binomial(n.astype(int), np.clip(p, 0, 1)).astype(float)
        sites.append((i, k, n, gidx))
    res = diffstats.beta_binomial_diff(sites, prior_weight=20.0, min_group_samples=2,
                                       ref_df=diffstats.REF_DF, site_weight="auto")
    return float(np.mean(np.array([r["p_value"] for r in res]) < 0.05))


def _low_stoich_universe(rng, m, nrep, depth=120):
    """A realistic modification-site universe: ~80% of sites are <1% modified, so most have an entire
    group at zero counts. Before the prior-collapse fix these all-zero-group sites (theta pinned at a
    Cox-Reid artifact ~0.37) set the across-site prior median and over-shrank every real site."""
    gidx = np.array([0] * nrep + [1] * nrep, dtype=int)
    sites = []
    for i in range(m):
        rate = rng.choice([0.002, 0.005, 0.008, 0.02], p=[0.45, 0.25, 0.13, 0.17])
        k = rng.poisson(rate * depth, 2 * nrep).astype(float)
        n = np.full(2 * nrep, float(depth))
        sites.append((("bg", i), np.minimum(k, n), n, gidx))
    return sites


def low_stoich_null_type_i(nrep=3, m=400, depth=120, seed=20260830):
    """Type-I at low stoichiometry: both groups share one low rate. Specificity must hold here too."""
    rng = np.random.default_rng(seed)
    res = diffstats.beta_binomial_diff(_low_stoich_universe(rng, m, nrep, depth), prior_weight=20.0)
    return float(np.mean(np.array([r["p_value"] for r in res]) < 0.05))


def low_stoich_sensitivity(nrep=3, depth=120, seed=20260830, p_thresh=1e-3):
    """POWER at realistic stoichiometry: clear, well-separated effect sites EMBEDDED in the
    low-stoichiometry universe must still be detected. This is the assertion the prior-collapse BLOCKER
    would have failed (embedded effects went p~1e-7 -> ~0.2): the old suite tested only mu 0.2-0.8 with
    no embedding and asserted only false positives, so it passed 4/4 while the headline analysis was
    dead on real data. Returns the fraction of effect sites reaching p < p_thresh."""
    rng = np.random.default_rng(seed)
    gidx = np.array([0] * nrep + [1] * nrep, dtype=int)
    effects = []
    for j, (mr, mt) in enumerate([(0.20, 0.72), (0.04, 0.42), (0.10, 0.55), (0.30, 0.80), (0.05, 0.45)]):
        kr = rng.binomial(depth, mr, nrep); kt = rng.binomial(depth, mt, nrep)
        k = np.concatenate([kr, kt]).astype(float); n = np.full(2 * nrep, float(depth))
        effects.append((("eff", j), k, n, gidx))
    sites = effects + _low_stoich_universe(rng, 400, nrep, depth)
    res = {r["key"]: r["p_value"] for r in diffstats.beta_binomial_diff(sites, prior_weight=20.0)}
    hit = sum(1 for j in range(len(effects)) if res.get(("eff", j), 1.0) < p_thresh)
    return hit / float(len(effects))


def all_zero_dropped_and_effect_kept():
    """M5: sites all-unmodified (or all-modified) in BOTH groups admit no between-group difference and
    must be DROPPED from the output entirely -- not emitted as "tested" p=1 rows that sit in the BH
    family, inflating m and crushing every adjusted p. A real effect embedded with them is still
    tested. Returns (n_all_zero_rows_emitted, effect_was_tested)."""
    g = np.array([0, 0, 0, 1, 1, 1], dtype=int)
    all_zero = [(("z", i), np.zeros(6), np.full(6, 100.0), g) for i in range(20)]
    effect = (("e",), np.array([20, 24, 26, 84, 90, 78], dtype=float), np.full(6, 120.0), g)
    res = diffstats.beta_binomial_diff(all_zero + [effect], prior_weight=20.0)
    keys = [r["key"] for r in res]
    n_zero = sum(1 for k in keys if isinstance(k, tuple) and k and k[0] == "z")
    return n_zero, (("e",) in keys)


def het_empty_last_col_null(reps=400, seed0=1):
    """MAJOR 1: heterogeneity null type-I when the last column is empty in every stratum (a phantom
    level). The two rows have DIFFERENT but stratum-CONSTANT proportions ([.8,.2] vs [.2,.8]) so the true
    risk difference is homogeneous (null holds) yet non-zero -- which is what makes the empty column bite:
    with BOTH column guards removed the rank-deficient theta vs Agresti-smoothed covariance inflates Q and
    this runs ~0.52 (an EQUAL-rows null gives theta~0 and stays low even when broken, so it could not
    catch a regression). informative_strata drops the globally-empty level and stratum_heterogeneity's
    min_col drops it again as a pooled near-empty column; together they hold this at ~nominal."""
    rng = np.random.default_rng(seed0)
    rows = ([0.8, 0.2], [0.2, 0.8])
    flags = used = 0
    for _ in range(reps):
        strata = []
        for _k in range(3):
            T = np.zeros((2, 3))
            for i in range(2):
                T[i, :2] = rng.multinomial(rng.integers(20, 60), rows[i])  # last col always 0
            strata.append(T)
        inf = gu.informative_strata(strata)
        if len(inf) < 2:
            continue
        used += 1
        _, hp, _, _ = gu.stratum_heterogeneity(inf)
        flags += (np.isfinite(hp) and hp < 0.05)
    return flags / max(1, used)


def het_near_empty_col_null(reps=500, seed0=3):
    """MAJOR 1 residual (min_col): heterogeneity null type-I when the last column is NEAR-empty -- a few
    reads pooled across strata, so informative_strata (which only drops levels empty in EVERY stratum)
    cannot help. Rows are homogeneous but unequal ([.8,.2]/[.2,.8]); the near-degenerate last column
    inflates Q to ~0.26 (~5x nominal) unless stratum_heterogeneity's pooled min_col guard drops it."""
    rng = np.random.default_rng(seed0)
    rows = ([0.8, 0.2], [0.2, 0.8])
    flags = used = 0
    for _ in range(reps):
        strata = []
        for _k in range(3):
            T = np.zeros((2, 3))
            for i in range(2):
                T[i, :2] = rng.multinomial(rng.integers(20, 60), rows[i])
                T[i, 2] = rng.integers(0, 2) if rng.random() < 0.3 else 0  # pooled < min_col
            strata.append(T)
        inf = gu.informative_strata(strata)
        if len(inf) < 2:
            continue
        used += 1
        _, hp, _, _ = gu.stratum_heterogeneity(inf)
        flags += (np.isfinite(hp) and hp < 0.05)
    return flags / max(1, used)


def het_informative_but_skipped_col_null(reps=600, seed0=11):
    """MAJOR 1 residual, one level deeper (min_col pooled over the min_row-passing strata): the last
    column's reads live entirely in a stratum the loop SKIPS (a row below min_row), so pooling min_col
    over ALL informative strata keeps that column even though it is empty in every USED stratum -> the
    near-empty degeneracy returns at full magnitude (~0.52). The two USED strata are identical unequal-row
    draws (true null). Guards that draw every row from integers(20,60) cannot reach this -- the skipped
    stratum here has a <min_row row on purpose."""
    rng = np.random.default_rng(seed0)
    rows = ([0.8, 0.2], [0.2, 0.8])
    flags = used = 0
    for _ in range(reps):
        strata = []
        for _k in range(2):                                   # 2 USED strata (pass min_row), last col empty
            T = np.zeros((2, 3))
            for i in range(2):
                T[i, :2] = rng.multinomial(rng.integers(20, 60), rows[i])
            strata.append(T)
        strata.append(np.array([[0., 0., 6.], [40., 40., 0.]]))  # row0 total 6 (<min_row) -> stratum skipped
        inf = gu.informative_strata(strata)
        if len(inf) < 2:
            continue
        used += 1
        _, hp, _, _ = gu.stratum_heterogeneity(inf)
        flags += (np.isfinite(hp) and hp < 0.05)
    return flags / max(1, used)


def het_col_drop_unmasks_skipped_col_null(reps=600, seed0=21):
    """MAJOR 1 residual, deepest level (fixed-point column reduction): the min_col drop itself lowers row
    sums, so a BOUNDARY stratum (row totals just >= min_row) that cleared min_col for a second column can
    fall under min_row once a FIRST near-empty column is dropped -- stranding that second column empty
    among the strata actually used. A single min_col pass (even pooled over the initially-min_row-passing
    strata) cannot see this; only recomputing the used set after each drop (iterating to a fixed point)
    closes it. Two identical main strata (homogeneous truth) + one boundary stratum (row totals 11)
    carrying col X (pooled 4 < min_col) and col Y (pooled 8 >= min_col): dropping X pushes the boundary
    under min_row, leaving Y empty among the used strata. Single-pass ~0.48; fixed point ~0.05."""
    rng = np.random.default_rng(seed0)
    rows = ([0.8, 0.2], [0.2, 0.8])
    flags = used = 0
    for _ in range(reps):
        strata = []
        for _k in range(2):                                   # main used strata; cols X,Y empty
            T = np.zeros((2, 4))
            for i in range(2):
                T[i, :2] = rng.multinomial(rng.integers(20, 60), rows[i])
            strata.append(T)
        strata.append(np.array([[3., 2., 2., 4.], [2., 3., 2., 4.]]))  # row totals 11; X pooled 4, Y pooled 8
        inf = gu.informative_strata(strata)
        if len(inf) < 2:
            continue
        used += 1
        _, hp, _, _ = gu.stratum_heterogeneity(inf)
        flags += (np.isfinite(hp) and hp < 0.05)
    return flags / max(1, used)


def main():
    # thresholds sit ~4 SE (SE~0.01 at n=500) above the observed Cox-Reid rates and WELL below the
    # ~0.11 pre-fix inflation, so this bites on a regression without flaking on Monte-Carlo noise.
    checks = [
        ("null type-I @ phi=0,    3v3 <= 0.08", null_type_i(0.0, 3) <= 0.08),
        ("null type-I @ phi=0.005,3v3 <= 0.09", null_type_i(0.005, 3) <= 0.09),
        ("null type-I @ phi=0.02, 3v3 <= 0.09", null_type_i(0.02, 3) <= 0.09),
        ("null type-I @ phi=0.05, 3v3 <= 0.09", null_type_i(0.05, 3) <= 0.09),
        # low-stoichiometry regime (the real m6A universe: ~80% of sites <1% modified). Specificity
        # must hold AND real effects embedded in it must still be detected (prior-collapse BLOCKER).
        ("low-stoich null type-I,  3v3 <= 0.09", low_stoich_null_type_i() <= 0.09),
        ("low-stoich SENSITIVITY,  3v3 >= 0.80", low_stoich_sensitivity() >= 0.80),
    ]
    _n_zero, _eff_kept = all_zero_dropped_and_effect_kept()
    checks += [
        ("M5: all-zero sites DROPPED (not in BH family)", _n_zero == 0),
        ("M5: a real effect among them is still tested", _eff_kept),
    ]
    # Shared single-stratum exact test (genotype_utils.run_contingency_test) -- used by 5 other scripts.
    _nm, _, _, _p2xc = gu.run_contingency_test(np.array([[50, 10, 40], [10, 50, 40]], dtype=float))
    _nm3, _, _, _pbig = gu.run_contingency_test(np.array([[1000, 0, 0], [0, 1000, 0], [0, 0, 1000]], dtype=float))
    # zero-row guard: a 3x2 stratum with an all-zero row (a transcript with no reads in that sample --
    # routine at 3v3, the design point) must DROP the zero row and test the rest, not return NaN and
    # silently drop a genuine hit from the BH family.
    _nmz, _, _, _pz = gu.montecarlo_exact_test(np.array([[0, 0], [10, 90], [85, 15]], dtype=float))
    _nmu, _, _, _pu = gu.montecarlo_exact_test(np.array([[0, 0], [0, 0], [50, 50]], dtype=float))
    checks += [
        ("shared exact: r x 2 null type-I [1,200,1] <= 0.06", mc_exact_rx2_null_type_i([1, 200, 1], 0.05) <= 0.06),
        ("shared exact: 2 x c handled (both margins), not floored/garbage", _nm.startswith("montecarlo_exact") and _p2xc < 1e-6),
        ("shared exact: p-floor removed (huge effect << 1e-4)", _pbig < 1e-4),
        ("shared exact: zero row DROPPED, not NaN'd out of the BH family", np.isfinite(_pz) and _pz < 1e-6),
        ("shared exact: <2x2 after dropping is untestable (NaN)", _nmu == "untestable" and not np.isfinite(_pu)),
    ]
    # run_contingency_test (the path all FIVE genotype scripts use) must ALSO reduce a zero-row/zero-col
    # table and test the survivors -- not screen the whole thing to untestable before dispatch, which
    # silently dropped a significant reduced table and disagreed with the direct-call path.
    _rc_name, _, _, _rc_p = gu.run_contingency_test(np.array([[0, 0], [10, 90], [85, 15]], dtype=float))
    _rc_col = gu.run_contingency_test(np.array([[50, 10, 0], [10, 50, 0]], dtype=float))
    _rc_u = gu.run_contingency_test(np.array([[0, 0], [0, 0], [50, 50]], dtype=float))
    checks += [
        ("run_contingency_test: zero row reduced + tested (not untestable)", np.isfinite(_rc_p) and _rc_p < 1e-6),
        ("run_contingency_test: zero column reduced + tested", np.isfinite(_rc_col[3]) and _rc_col[3] < 1e-6),
        ("run_contingency_test: genuine <2x2 collapse still untestable", _rc_u[0] == "untestable" and not np.isfinite(_rc_u[3])),
    ]
    # MAJOR 1+2: informative_strata drops levels zero in EVERY stratum, so CMH + heterogeneity get the
    # reduction (not just the exact_single path). A phantom empty column made heterogeneity anti-
    # conservative and cost CMH a spurious df.
    _phantom = [np.array([[80, 20, 0], [20, 80, 0]], dtype=float) for _ in range(3)]  # last col empty everywhere
    _red = gu.informative_strata(_phantom)
    _cmh_full = gu.cmh_stratified_test([np.array([[80, 20], [20, 80]], dtype=float) for _ in range(3)])
    _cmh_phantom = gu.cmh_stratified_test(_red)
    checks += [
        ("MAJOR1: heterogeneity null w/ empty last col <= 0.10 (unequal rows; ~0.52 if unguarded)", het_empty_last_col_null() <= 0.10),
        ("MAJOR1 residual: heterogeneity null w/ NEAR-empty col <= 0.12 (~0.26 without min_col)", het_near_empty_col_null() <= 0.12),
        ("MAJOR1 residual: heterogeneity null w/ INFORMATIVE-but-min_row-SKIPPED col <= 0.10 (~0.52 if pooled over all)", het_informative_but_skipped_col_null() <= 0.10),
        ("MAJOR1 residual: heterogeneity null w/ col-drop UNMASKING a skipped stratum <= 0.10 (~0.48 single-pass; needs fixed point)", het_col_drop_unmasks_skipped_col_null() <= 0.10),
        ("MAJOR1+2: informative_strata drops the globally-empty level", all(T.shape == (2, 2) for T in _red)),
        ("MAJOR2: phantom empty level does not inflate CMH p (df restored)", abs(_cmh_phantom[3] - _cmh_full[3]) < 1e-9),
    ]
    n_pass = 0
    for name, ok in checks:
        rate_ok = "PASS" if ok else "**FAIL**"
        print(f"  [{rate_ok}] {name}")
        n_pass += ok
    print(f"\ndiffstats null calibration: {n_pass}/{len(checks)} checks passed")
    sys.exit(0 if n_pass == len(checks) else 1)


if __name__ == "__main__":
    main()
