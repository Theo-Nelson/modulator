#!/usr/bin/env python3
"""Guard finding L: TES clustering must never split one APA site into two partitions.

Finding L (AUDIT_2136a74): under the old single-linkage `cluster_positions`, two reads of the SAME
intron chain whose 3' ends were within `apa_window` could land in different clusters, becoming two
fragmentforms with distinct `zn_index` in the same metagene. `test_diffs` then compared an isoform
against ITSELF (measured: 2490 such pairs on real HG002 chr21-scale data). The mode-seeking rewrite
(4b181fb) fixes this because each fragmentform's TES IS a cluster representative, and the invariant
below guarantees representatives are pairwise MORE than `window` apart -- so two fragmentforms in one
(chrom, strand) have TES either identical (same cluster) or > window apart, never "within window but
distinct". This test would FAIL under the old single-linkage clustering.

Run: <modulator-env>/bin/python resources/synthetic_3exon/test_tes_clustering_invariant.py
"""
from __future__ import annotations
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
from assemble_transcripts import cluster_positions  # noqa: E402


def _lcg(seed):
    """Deterministic PRNG (no Math.random-style nondeterminism, so failures reproduce)."""
    x = seed & 0xFFFFFFFF
    while True:
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        yield x


def main():
    checks = []

    def check(name, ok, detail=""):
        checks.append(ok)
        print(f"  {'PASS' if ok else '**FAIL**'}  {name}{('  --  ' + detail) if detail else ''}")

    rng = _lcg(1337)
    windows = [5, 10, 20, 40, 60, 100]

    rep_gap_violations = 0        # two cluster reps within `window` -> the finding-L failure mode
    coverage_violations = 0       # a position not within `window` of its own cluster rep
    partition_violations = 0      # clusters not a partition of the input multiset
    monotone_violations = 0       # cluster count increases as window grows (single-linkage pathology)

    TRIALS = 3000
    for t in range(TRIALS):
        n = 1 + next(rng) % 40
        # tight, overlapping positions in a small span so within-window collisions are common
        positions = sorted((next(rng) % 200) for _ in range(n))

        prev_count = None
        for w in windows:
            clusters = cluster_positions(positions, w)
            reps = [c["rep"] for c in clusters]

            # (1) representatives strictly more than `window` apart
            sreps = sorted(reps)
            if any(sreps[i + 1] - sreps[i] <= w for i in range(len(sreps) - 1)):
                rep_gap_violations += 1

            # (2) every member within `window` of its cluster rep
            for c in clusters:
                if any(abs(p - c["rep"]) > w for p in c["positions"]):
                    coverage_violations += 1
                    break

            # (3) clusters partition the input multiset exactly
            got = Counter()
            for c in clusters:
                got.update(c["positions"])
            if got != Counter(positions):
                partition_violations += 1

            # (4) cluster count is non-increasing as the window grows
            if prev_count is not None and len(clusters) > prev_count:
                monotone_violations += 1
            prev_count = len(clusters)

    check(f"cluster reps pairwise > window apart ({TRIALS} random multisets x {len(windows)} windows)",
          rep_gap_violations == 0, f"{rep_gap_violations} violations")
    check("every position within window of its cluster rep", coverage_violations == 0,
          f"{coverage_violations} violations")
    check("clusters exactly partition the input multiset", partition_violations == 0,
          f"{partition_violations} violations")
    check("cluster count non-increasing as window grows", monotone_violations == 0,
          f"{monotone_violations} violations")

    # The concrete finding-L scenario: same chain, two 3' ends 10 nt apart (< apa_window=20), the
    # far one with far more support -> must collapse to ONE cluster (one fragmentform), not two.
    cl = cluster_positions(sorted([39342305] * 106 + [39342315] * 1888), 20)
    check("finding-L pair (TES 10nt apart, apa_window=20) collapses to ONE cluster",
          len(cl) == 1, f"got {len(cl)} clusters, reps={[c['rep'] for c in cl]}")

    n_ok = sum(checks)
    print(f"\nTES clustering invariant: {n_ok}/{len(checks)} checks passed")
    sys.exit(0 if n_ok == len(checks) else 1)


if __name__ == "__main__":
    main()
