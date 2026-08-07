#!/usr/bin/env python3
"""Differential poly(A) TAIL LENGTH between conditions, replicate-aware.

Consumes the per-read tail table (build_read_polya_table.py) + the samplesheet metadata and asks,
per fragmentform (or per gene): is the poly(A) tail length different between two conditions?

Tail length is CONTINUOUS, so this is the one between-condition test that is not beta-binomial.
Each REPLICATE is reduced to a single summary (its median tail for that feature) and the two groups
of replicate summaries are compared with Welch's t-test -- see diffstats.continuous_diff.

Why not a per-read Mann-Whitney (as test_taillength_diffs.py uses)? Because that is the right test
for a different question. Comparing the fragmentforms *within* one pooled population treats reads as
the unit legitimately; comparing CONDITIONS does not -- there the biological unit is the replicate,
and pooling 27M reads across 3 replicates would be pseudoreplication (arbitrarily small p-values for
trivial differences). Hence: one number per replicate, then test across replicates.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

import diffstats
from genotype_utils import benjamini_hochberg

OUT_COLS = ["contrast", "level", "feature", "gene_name", "n_reference", "n_test",
            "reads_reference", "reads_test", "median_tail_reference", "median_tail_test",
            "delta_nt", "stat", "p_value", "p_adj_bh", "per_sample_json"]


def parse_args():
    ap = argparse.ArgumentParser(description="Replicate-aware differential poly(A) tail length between conditions.")
    ap.add_argument("--tail-tsv", required=True, help="*_read_tail_lengths.tsv (per-read, has sample)")
    ap.add_argument("--sample-metadata", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--level", choices=["fragmentform", "gene"], default="fragmentform")
    ap.add_argument("--contrast-name", default="")
    ap.add_argument("--column", default="condition")
    ap.add_argument("--test", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--min-reads-per-sample", type=int, default=10,
                    help="min reads for a feature in EVERY sample (a median needs support)")
    ap.add_argument("--min-samples-per-group", type=int, default=2)
    ap.add_argument("--min-tail", type=int, default=1, help="drop reads with tail_len < this (pt:i:0 = no estimate)")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def main():
    args = parse_args()
    name = args.contrast_name or f"{args.test}_vs_{args.reference}"

    meta = pd.read_csv(args.sample_metadata, sep="\t", low_memory=False)
    if "sample" not in meta.columns or args.column not in meta.columns:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return
    grp = dict(zip(meta["sample"].astype(str), meta[args.column].astype(str)))
    ref_s = {s for s, g in grp.items() if g == args.reference}
    test_s = {s for s, g in grp.items() if g == args.test}
    if len(ref_s) < args.min_samples_per_group or len(test_s) < args.min_samples_per_group:
        print(f"[condition_tail] {name}: need >={args.min_samples_per_group}/group "
              f"({len(ref_s)} vs {len(test_s)})", file=sys.stderr, flush=True)
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return

    want = ["sample", "tail_len", "ZT", "gene_name"]
    hdr = pd.read_csv(args.tail_tsv, sep="\t", nrows=0).columns
    df = pd.read_csv(args.tail_tsv, sep="\t", low_memory=False, usecols=[c for c in want if c in hdr])
    df["sample"] = df["sample"].astype(str)
    df = df[df["sample"].isin(ref_s | test_s)]
    df = df[pd.to_numeric(df["tail_len"], errors="coerce").fillna(0) >= int(args.min_tail)]
    if df.empty:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return
    feat_col = "ZT" if args.level == "fragmentform" else "gene_name"
    df = df[df[feat_col].astype(str).ne("")]

    # One summary per (feature, replicate): the median tail. The replicate is the unit of analysis.
    per = df.groupby([feat_col, "sample"], sort=False)["tail_len"].agg(["median", "size"]).reset_index()
    per = per[per["size"] >= args.min_reads_per_sample]
    gene_of = (df.groupby(feat_col, sort=False)["gene_name"].first()
               if "gene_name" in df.columns and args.level == "fragmentform" else None)

    rows = []
    for feature, sub in per.groupby(feat_col, sort=False):
        m = dict(zip(sub["sample"], sub["median"]))
        n = dict(zip(sub["sample"], sub["size"]))
        a = [m[s] for s in ref_s if s in m]
        b = [m[s] for s in test_s if s in m]
        if len(a) < args.min_samples_per_group or len(b) < args.min_samples_per_group:
            continue
        r = diffstats.continuous_diff(a, b)
        if not np.isfinite(r["p_value"]):
            continue
        rows.append({
            "contrast": name, "level": args.level, "feature": feature,
            "gene_name": (gene_of.get(feature, "") if gene_of is not None else feature),
            "n_reference": r["n_reference"], "n_test": r["n_test"],
            "reads_reference": int(sum(n[s] for s in ref_s if s in n)),
            "reads_test": int(sum(n[s] for s in test_s if s in n)),
            "median_tail_reference": round(r["mean_reference"], 2),
            "median_tail_test": round(r["mean_test"], 2),
            "delta_nt": round(r["delta"], 2), "stat": round(r["stat"], 4), "p_value": r["p_value"],
            "per_sample_json": json.dumps({s: [round(float(m[s]), 1), int(n[s])] for s in sorted(m)},
                                          separators=(",", ":")),
        })
    out = pd.DataFrame(rows)
    if out.empty:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return
    out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
    out["_abs"] = out["delta_nt"].abs()
    out = out.sort_values(["p_adj_bh", "_abs"], ascending=[True, False]).drop(columns="_abs")
    out = out[OUT_COLS].reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    if args.verbose:
        print(f"[condition_tail:{args.level}] {name}: {len(out):,} features tested, "
              f"{int((out['p_adj_bh'] < 0.05).sum()):,} at FDR<0.05 -> {args.out_tsv}", flush=True)


if __name__ == "__main__":
    main()
