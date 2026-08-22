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


def main():
    # thresholds sit ~4 SE (SE~0.01 at n=500) above the observed Cox-Reid rates and WELL below the
    # ~0.11 pre-fix inflation, so this bites on a regression without flaking on Monte-Carlo noise.
    checks = [
        ("null type-I @ phi=0,    3v3 <= 0.08", null_type_i(0.0, 3) <= 0.08),
        ("null type-I @ phi=0.005,3v3 <= 0.09", null_type_i(0.005, 3) <= 0.09),
        ("null type-I @ phi=0.02, 3v3 <= 0.09", null_type_i(0.02, 3) <= 0.09),
        ("null type-I @ phi=0.05, 3v3 <= 0.09", null_type_i(0.05, 3) <= 0.09),
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
