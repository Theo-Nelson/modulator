#!/usr/bin/env python3
"""
Calibrate ``between_conditions.ref_df`` for a dataset.

WHY THIS EXISTS
---------------
The between-condition count tests (differential modification / isoform / APA /
junction usage) use a beta-binomial LRT whose reference distribution is
``F(1, ref_df)`` (see ``diffstats.py``). ``chi2(1)`` (i.e. ref_df -> inf) is
anti-conservative at small replicate counts because the dispersion is *estimated*;
``F(1, ref_df)`` absorbs that uncertainty. The right ``ref_df`` therefore depends
on how reproducible YOUR replicates are, so it must be checked per dataset.

THE CHECK (simple, empirical)
-----------------------------
Build a NULL contrast in which there is, by construction, no real biological
difference -- two arbitrary halves of the replicates OF THE SAME CONDITION -- and
run the real differential-modification engine on it across a grid of ``ref_df``.
Under the null the fraction of sites with ``p < 0.05`` should be ~0.05. Pick the
smallest ``ref_df`` whose null rate is not anti-conservative (<= ~0.06); a larger
value is safe but costs sensitivity.

    ref_df too small  -> null p<0.05 rate >> 0.05  (anti-conservative: false positives)
    ref_df too large  -> null p<0.05 rate << 0.05  (over-conservative: lost power)

USAGE
-----
    python workflow/scripts/calibrate_between_conditions.py \
        --zn-long   results/aggregate_zn/<prefix>_FILTERED_sites_long.tsv \
        --sample-metadata results/<prefix>_sample_metadata.tsv \
        [--column condition] [--min-cov 20] [--ref-df-grid 4,6,8,10,12,15] \
        [--target 0.05] [--tolerance 0.01]

Needs a condition with >=4 replicates (to split 2 vs 2). Reuses
``test_condition_mod_diffs.py`` unchanged, so the null uses the exact production
engine.
"""
from __future__ import annotations
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def build_null_metadata(meta: pd.DataFrame, column: str) -> tuple[pd.DataFrame, str, str, list[str]]:
    """Split the largest condition group into two arbitrary halves -> a null contrast."""
    counts = meta[column].value_counts()
    grp = counts.index[0]
    members = sorted(meta.loc[meta[column] == grp, "sample"].tolist())
    if len(members) < 4:
        sys.exit(f"[calibrate] largest group '{grp}' has only {len(members)} replicate(s); "
                 f"need >=4 to build a 2-vs-2 within-condition null. "
                 f"(Add replicates, or calibrate on a dataset that has them.)")
    half = len(members) // 2
    a, b = members[:half], members[half:]
    null = meta[meta[column] == grp].copy()
    null["null_group"] = ["A" if s in a else "B" for s in null["sample"]]
    return null, "A", "B", members


def run_one(zn_long: str, null_meta_path: str, out_path: str, ref_df: int,
            min_cov: int, prior_weight: float) -> pd.DataFrame:
    cmd = [sys.executable, str(HERE / "test_condition_mod_diffs.py"),
           "--in-tsv", zn_long, "--sample-metadata", null_meta_path,
           "--out-tsv", out_path, "--column", "null_group",
           "--test", "B", "--reference", "A", "--contrast-name", f"null_refdf{ref_df}",
           "--min-cov", str(min_cov), "--min-samples-per-group", "2",
           "--prior-weight", str(prior_weight), "--ref-df", str(ref_df)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return pd.read_csv(out_path, sep="\t")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zn-long", required=True, help="<prefix>_FILTERED_sites_long.tsv")
    ap.add_argument("--sample-metadata", required=True)
    ap.add_argument("--column", default="condition")
    ap.add_argument("--min-cov", type=int, default=20)
    ap.add_argument("--prior-weight", type=float, default=20.0)
    ap.add_argument("--ref-df-grid", default="4,6,8,10,12,15")
    ap.add_argument("--target", type=float, default=0.05, help="desired null p<0.05 rate")
    ap.add_argument("--tolerance", type=float, default=0.01,
                    help="accept ref_df whose null rate <= target+tolerance")
    args = ap.parse_args()

    meta = pd.read_csv(args.sample_metadata, sep="\t")
    if "sample" not in meta.columns or args.column not in meta.columns:
        sys.exit(f"[calibrate] metadata needs 'sample' and '{args.column}' columns "
                 f"(found {list(meta.columns)})")
    null_meta, test_lvl, ref_lvl, members = build_null_metadata(meta, args.column)
    grid = [int(x) for x in args.ref_df_grid.split(",") if x.strip()]

    print(f"[calibrate] null contrast = within-condition split of "
          f"{len(members)} replicates: {test_lvl} vs {ref_lvl}")
    print(f"[calibrate] {'ref_df':>7}  {'n_sites':>8}  {'null p<0.05':>12}  verdict")

    with tempfile.TemporaryDirectory() as td:
        nm = Path(td) / "null_metadata.tsv"
        null_meta.to_csv(nm, sep="\t", index=False)
        rows = []
        for rdf in grid:
            out = Path(td) / f"null_{rdf}.tsv"
            res = run_one(args.zn_long, str(nm), str(out), rdf, args.min_cov, args.prior_weight)
            n = len(res)
            rate = float((res["p_value"] < 0.05).mean()) if n else float("nan")
            calibrated = (rate == rate) and rate <= args.target + args.tolerance
            verdict = ("OK (calibrated)" if calibrated
                       else "anti-conservative" if (rate == rate and rate > args.target + args.tolerance)
                       else "no testable sites")
            print(f"[calibrate] {rdf:>7}  {n:>8}  {rate:>12.3f}  {verdict}")
            rows.append((rdf, n, rate, calibrated))

    ok = [r for r in rows if r[3]]
    if ok:
        rec = min(ok, key=lambda r: r[0])   # smallest ref_df that is calibrated
        print(f"\n[calibrate] RECOMMENDATION: set between_conditions.ref_df={rec[0]} "
              f"(null p<0.05 = {rec[2]:.3f}, ~target {args.target}).")
        print(f"[calibrate] A larger ref_df is safe but costs sensitivity; a smaller one "
              f"is anti-conservative on this data.")
    else:
        print(f"\n[calibrate] No ref_df in the grid reached the target null rate "
              f"(<= {args.target + args.tolerance}). Try a larger grid (e.g. 20,30,50) "
              f"or check that the null contrast has enough covered sites.")


if __name__ == "__main__":
    main()
