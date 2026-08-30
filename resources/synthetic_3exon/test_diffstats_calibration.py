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
    n_pass = 0
    for name, ok in checks:
        rate_ok = "PASS" if ok else "**FAIL**"
        print(f"  [{rate_ok}] {name}")
        n_pass += ok
    print(f"\ndiffstats null calibration: {n_pass}/{len(checks)} checks passed")
    sys.exit(0 if n_pass == len(checks) else 1)


if __name__ == "__main__":
    main()
