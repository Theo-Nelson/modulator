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
            "reads_reference", "reads_test", "mean_tail_reference", "mean_tail_test",
            "delta_nt", "stat", "p_value", "p_adj_bh", "per_sample_json", "per_replicate_json"]


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
    # Ensure the output directory exists BEFORE any early-return writes an empty table -- the
    # pipeline does not pre-create between_conditions/, so a no-result contrast otherwise crashes.
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    name = args.contrast_name or f"{args.test}_vs_{args.reference}"

    meta = pd.read_csv(args.sample_metadata, sep="\t", low_memory=False, keep_default_na=False)
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
    # pt:i:0 = dorado "no estimate", not a 0-nt tail: require an actual estimate (tail_estimated).
    _tl = pd.to_numeric(df["tail_len"], errors="coerce").fillna(0)
    if "tail_estimated" in df.columns:
        _est = df["tail_estimated"].astype(str).str.strip().str.lower().isin(("true", "1", "yes"))
    else:
        _est = _tl > 0
    df = df[_est & (_tl >= int(args.min_tail))]
    if df.empty:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return
    feat_col = "ZT" if args.level == "fragmentform" else "gene_name"
    df = df[df[feat_col].astype(str).ne("")]

    # One summary per (feature, replicate); the replicate is the unit of analysis.
    gene_untestable = []       # gene level: genes dropped by the common-fragmentform restriction
    ff_all = None
    if args.level == "fragmentform":
        per = df.groupby([feat_col, "sample"], sort=False)["tail_len"].agg(["median", "size"]).reset_index()
        per = per[per["size"] >= args.min_reads_per_sample]
    else:
        # GENE level: a replicate's summary is the mean of its per-FRAGMENTFORM median tails (each
        # fragmentform weighted EQUALLY), NOT the pooled-read median. Pooling reads weights each
        # fragmentform by its usage, so a condition-linked isoform-usage shift moves the gene median even
        # when no isoform's tail changes -- reporting e.g. -35.5 nt "significant" on a gene whose isoforms
        # are individually flat. The fragmentform-averaged summary is usage-invariant: a pure usage shift
        # leaves each fragmentform's own median unchanged, so the gene delta is ~0; a real gene-wide tail
        # shift (all isoforms move) is preserved. Only fragmentforms with >= min_reads in a replicate
        # count toward that replicate's mean.
        _dff = df[df["ZT"].astype(str).ne("")] if "ZT" in df.columns else df.iloc[0:0]
        ff_all = _dff.groupby(["gene_name", "ZT", "sample"], sort=False)["tail_len"].agg(["median", "size"]).reset_index()
        ff_all = ff_all[ff_all["size"] >= args.min_reads_per_sample]
        # A fragmentform must clear the read threshold in EVERY sample that summarises the gene, else
        # WHICH forms enter a replicate's mean changes between conditions: a form dropping below the
        # threshold in one condition (a pure isoform-usage shift -- exactly the case this test targets)
        # silently leaves that condition's mean, turning the confound into a STEP at the threshold and
        # manufacturing a tail difference where no isoform's tail moved. Restrict to the fragmentforms
        # common to all contributing samples so every replicate's mean averages the SAME form set; a
        # genuine gene-wide shift (all shared forms move) still survives.
        if not ff_all.empty:
            n_samples_per_gene = ff_all.groupby("gene_name")["sample"].transform("nunique")
            n_samples_per_zt = ff_all.groupby(["gene_name", "ZT"])["sample"].transform("nunique")
            ff = ff_all[n_samples_per_zt == n_samples_per_gene]
        else:
            ff = ff_all
        # A gene with qualifying reads but NO fragmentform covered in all its samples cannot form a
        # composition-stable gene mean, so it is not tested here. It must NOT vanish silently (the loss
        # is biased toward multi-isoform genes with uneven coverage -- the population this analysis is
        # about): emit an explicit untestable row (NaN stats) for each so it stays visible and counted.
        gene_untestable = sorted(set(ff_all["gene_name"]) - set(ff["gene_name"])) if not ff_all.empty else []
        per = (ff.groupby(["gene_name", "sample"], sort=False)
                 .agg(median=("median", "mean"), size=("size", "sum")).reset_index())
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
            # these are the MEAN over replicates of each replicate's median tail (continuous_diff/Welch
            # works on replicate summaries), so name them mean_* -- delta_nt = mean_test - mean_reference.
            "mean_tail_reference": round(r["mean_reference"], 2),
            "mean_tail_test": round(r["mean_test"], 2),
            "delta_nt": round(r["delta"], 2), "stat": round(r["stat"], 4), "p_value": r["p_value"],
            "per_sample_json": json.dumps({s: [round(float(m[s]), 1), int(n[s])] for s in sorted(m)},
                                          separators=(",", ":")),
            # sort the (set-valued) ref_s/test_s so the JSON key order is deterministic across runs --
            # otherwise string-hash randomization (PYTHONHASHSEED) reorders these dicts and the output
            # TSV is not byte-reproducible.
            "per_replicate_json": json.dumps(
                {"reference": {s: round(float(m[s]), 1) for s in sorted(ref_s) if s in m},
                 "test": {s: round(float(m[s]), 1) for s in sorted(test_s) if s in m}},
                separators=(",", ":")),
        })
    # Surface genes the common-fragmentform restriction removed (only those that WOULD have been
    # testable -- qualifying coverage in >= min_samples_per_group replicates of BOTH conditions) as
    # explicit untestable rows (NaN stats), so composition-unstable multi-isoform genes are visible
    # and counted rather than silently absent. per_replicate_json carries the reason.
    n_untestable = 0
    for g in gene_untestable:
        sub = ff_all[ff_all["gene_name"] == g]
        ref_reps = sorted({s for s in sub["sample"] if s in ref_s})
        test_reps = sorted({s for s in sub["sample"] if s in test_s})
        if len(ref_reps) < args.min_samples_per_group or len(test_reps) < args.min_samples_per_group:
            continue
        cov = {}
        for _, rr in sub.iterrows():
            cov.setdefault(str(rr["sample"]), []).append(str(rr["ZT"]).split(".")[-1])
        rows.append({
            "contrast": name, "level": args.level, "feature": g, "gene_name": g,
            "n_reference": len(ref_reps), "n_test": len(test_reps),
            "reads_reference": int(sub[sub["sample"].isin(ref_s)]["size"].sum()),
            "reads_test": int(sub[sub["sample"].isin(test_s)]["size"].sum()),
            "mean_tail_reference": float("nan"), "mean_tail_test": float("nan"),
            "delta_nt": float("nan"), "stat": float("nan"), "p_value": float("nan"),
            "per_sample_json": json.dumps({s: sorted(cov[s]) for s in sorted(cov)}, separators=(",", ":")),
            "per_replicate_json": json.dumps(
                {"status": "untestable_no_common_fragmentform",
                 "note": "no fragmentform cleared min_reads_per_sample in every replicate; "
                         "a composition-stable gene mean cannot be formed"},
                separators=(",", ":")),
        })
        n_untestable += 1
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
        n_tested = int(np.isfinite(pd.to_numeric(out["p_value"], errors="coerce")).sum())
        msg = (f"[condition_tail:{args.level}] {name}: {n_tested:,} features tested, "
               f"{int((out['p_adj_bh'] < 0.05).sum()):,} at FDR<0.05")
        if gene_untestable:
            msg += (f"; {n_untestable:,} gene(s) NOT tested (no fragmentform covered in every "
                    f"replicate -> composition-unstable; emitted as untestable rows)")
        print(msg + f" -> {args.out_tsv}", flush=True)


if __name__ == "__main__":
    main()
