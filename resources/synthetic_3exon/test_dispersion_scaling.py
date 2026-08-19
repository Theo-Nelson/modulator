#!/usr/bin/env python3
"""
Guard the between-condition dispersion-shrinkage scaling fix (diffstats.beta_binomial_diff).

The legacy fixed per-site weight (w=1 vs prior_weight=20) shrank every site's dispersion ~95%
toward the near-binomial global median, regardless of cohort size. On a large or heterogeneous
cohort that crushes genuinely OVERDISPERSED sites onto the bulk and reads their replicate scatter
as signal -> false positives (up to ~50-100% FPR at those sites). site_weight="auto" (= N_site-2)
makes the shrinkage fade as the cohort grows, so at cohort scale overdispersed true-null sites
return to a nominal false-positive rate.

This test simulates a TRUE-NULL cohort (same mean in both groups) with a minority of overdispersed
sites and asserts: (1) 'auto' drives the overdispersed-site FPR to ~nominal at large N, well below
the legacy w=1; (2) clean (near-binomial) sites stay nominal under both.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow" / "scripts"))
import diffstats  # noqa: E402


def make_null_cohort(rng, N, n_clean=450, n_over=50, phi_over=0.10, cov=60):
    g = np.array([0] * (N // 2) + [1] * (N // 2))
    sites, over = [], []
    for i in range(n_clean + n_over):
        mu = rng.uniform(0.2, 0.8)
        is_over = i >= n_clean
        if is_over:
            th = (1 - phi_over) / phi_over
            p = rng.beta(mu * th, (1 - mu) * th, size=N)
        else:
            p = np.full(N, mu)                        # binomial (near-zero dispersion)
        sites.append((f"s{i}", rng.binomial(cov, p).astype(float), np.full(N, cov, float), g))
        over.append(is_over)
    return sites, np.array(over)


def fpr(sites, over, site_weight):
    res = diffstats.beta_binomial_diff(sites, prior_weight=20.0, ref_df=10, site_weight=site_weight)
    p = np.array([r["p_value"] for r in res])
    return float(np.mean(p[over] < 0.05)), float(np.mean(p[~over] < 0.05))


def main():
    rng = np.random.default_rng(11)
    N = 200                                            # cohort scale (100 vs 100)
    sites, over = make_null_cohort(rng, N)
    over_legacy, clean_legacy = fpr(sites, over, 1.0)
    over_auto, clean_auto = fpr(sites, over, "auto")

    print(f"  N={N} ({N//2}v{N//2}), true null, 10% sites overdispersed (phi=0.10)")
    print(f"  {'':16}{'overdispersed FPR':>20}{'clean FPR':>14}")
    print(f"  {'legacy w=1':16}{over_legacy:>20.3f}{clean_legacy:>14.3f}")
    print(f"  {'auto (N-2)':16}{over_auto:>20.3f}{clean_auto:>14.3f}")

    checks = [
        ("auto restores near-nominal FPR at overdispersed sites (<0.12)", over_auto < 0.12),
        ("auto is a large improvement over legacy at overdispersed sites", over_auto < over_legacy - 0.15),
        ("clean sites stay nominal under auto (<0.10)", clean_auto < 0.10),
        ("legacy is inflated at overdispersed sites (>0.20) [documents the bug]", over_legacy > 0.20),
    ]
    n_fail = 0
    print()
    for msg, ok in checks:
        print(f"  [{'PASS' if ok else '**FAIL**'}] {msg}")
        n_fail += (not ok)
    print(f"\ndispersion scaling: {len(checks) - n_fail}/{len(checks)} checks passed")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
