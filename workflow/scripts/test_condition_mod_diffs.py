#!/usr/bin/env python3
"""Differential modification BETWEEN CONDITIONS (e.g. zikv vs mock), replicate-aware.

Consumes the per-sample ZN long table (``*_FILTERED_sites_long.tsv``: one row per
site x transcript-partition x sample with Nvalid_cov / Nmod) plus the samplesheet-derived metadata
(sample -> condition), and asks, per modification site: is the modified fraction different between
the two conditions, ACCOUNTING FOR REPLICATE VARIABILITY?

This is deliberately NOT the same shape as test_stoichiometry_diffs.py, which pools samples to
compare transcript partitions within one population. Pooling reads across replicates for a
BETWEEN-CONDITION question is pseudoreplication: with millions of reads but n=3 per group it returns
p=1e-300 for trivial differences (measured: 62% of simulated NULL sites reach p<0.05). The biological
unit here is the replicate, so the test is a beta-binomial LRT with dispersion shrinkage -- see
diffstats.py for the model, the calibration, and the validation numbers.

Counts are summed over transcript partitions first: the question is about the SITE, not the isoform.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

import diffstats
from genotype_utils import benjamini_hochberg

_WANT = ["sample", "ZN_transcript_index", "chrom", "start0", "end0", "strand", "mod_code",
         "Nvalid_cov", "Nmod", "gene_name"]
# ZN_transcript_index is inserted after gene_name only in --by-transcript mode (see build_out_cols).
OUT_COLS = ["contrast", "chrom", "start0", "end0", "strand", "mod_code", "gene_name",
            "n_reference", "n_test", "reads_reference", "reads_test",
            "mu_reference", "mu_test", "delta", "dispersion", "lrt_stat", "p_value", "p_adj_bh",
            "per_replicate_json"]


def build_out_cols(by_transcript):
    if not by_transcript:
        return list(OUT_COLS)
    cols = list(OUT_COLS)
    cols.insert(cols.index("gene_name") + 1, "ZN_transcript_index")
    return cols


def parse_args():
    ap = argparse.ArgumentParser(description="Replicate-aware differential modification between two conditions.")
    ap.add_argument("--in-tsv", required=True, help="*_FILTERED_sites_long.tsv (per-sample site counts)")
    ap.add_argument("--sample-metadata", required=True, help="TSV with sample + condition columns")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--contrast-name", default="")
    ap.add_argument("--column", default="condition", help="metadata column holding the group")
    ap.add_argument("--test", required=True, help="the test level (e.g. zikv)")
    ap.add_argument("--reference", required=True, help="the reference level (e.g. mock)")
    ap.add_argument("--min-cov", type=int, default=20, help="min Nvalid_cov per sample for a site to be tested")
    ap.add_argument("--min-samples-per-group", type=int, default=2)
    ap.add_argument("--prior-weight", type=float, default=20.0, help="dispersion shrinkage strength (diffstats)")
    ap.add_argument("--ref-df", type=int, default=diffstats.REF_DF, help="F(1, df) reference; see diffstats.py")
    ap.add_argument("--site-weight", default="auto",
                    help="per-site weight in dispersion shrinkage: 'auto' = N_site-2 (scales with cohort "
                         "size, so heterogeneous large cohorts aren't over-shrunk), or a fixed number "
                         "(1 = legacy behaviour).")
    ap.add_argument("--mod-filter", nargs="*", default=None, help="restrict to these mod codes")
    ap.add_argument("--by-transcript", action="store_true",
                    help="Test each site PER TRANSCRIPT PARTITION (ZN) instead of summing partitions: "
                         "compares the SAME transcript between conditions. Lower per-test coverage but "
                         "resolves which fragmentform carries the change.")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def main():
    args = parse_args()
    name = args.contrast_name or f"{args.test}_vs_{args.reference}"

    meta = pd.read_csv(args.sample_metadata, sep="\t", low_memory=False, keep_default_na=False)
    if "sample" not in meta.columns or args.column not in meta.columns:
        print(f"[condition_mod] metadata needs 'sample' and {args.column!r}", file=sys.stderr, flush=True)
        pd.DataFrame(columns=build_out_cols(args.by_transcript)).to_csv(args.out_tsv, sep="\t", index=False)
        return
    grp = dict(zip(meta["sample"].astype(str), meta[args.column].astype(str)))
    ref_s = {s for s, g in grp.items() if g == args.reference}
    test_s = {s for s, g in grp.items() if g == args.test}
    if len(ref_s) < args.min_samples_per_group or len(test_s) < args.min_samples_per_group:
        print(f"[condition_mod] {name}: need >={args.min_samples_per_group} samples per group "
              f"(reference={len(ref_s)}, test={len(test_s)}); nothing to do", file=sys.stderr, flush=True)
        pd.DataFrame(columns=build_out_cols(args.by_transcript)).to_csv(args.out_tsv, sep="\t", index=False)
        return

    hdr = pd.read_csv(args.in_tsv, sep="\t", nrows=0).columns
    df = pd.read_csv(args.in_tsv, sep="\t", low_memory=False, usecols=[c for c in _WANT if c in hdr])
    df["sample"] = df["sample"].astype(str)
    df = df[df["sample"].isin(ref_s | test_s)]
    if args.mod_filter:
        df = df[df["mod_code"].isin(set(args.mod_filter))]
    if df.empty:
        pd.DataFrame(columns=build_out_cols(args.by_transcript)).to_csv(args.out_tsv, sep="\t", index=False)
        return

    # Genomic key. In --by-transcript mode the transcript partition (ZN) joins the key, so each
    # (site, ZN) is tested independently -- "the SAME transcript between conditions". Otherwise the
    # partitions are summed and the question is about the site.
    by_tx = args.by_transcript and "ZN_transcript_index" in df.columns
    base_key = ["chrom", "start0", "end0", "strand", "mod_code"]
    key = base_key + (["ZN_transcript_index"] if by_tx else [])
    out_cols = build_out_cols(by_tx)

    agg = df.groupby(key + ["sample"], sort=False, observed=True)[["Nvalid_cov", "Nmod"]].sum().reset_index()
    genes = df.groupby(key, sort=False, observed=True)["gene_name"].first() if "gene_name" in df.columns else None

    cov = agg.pivot_table(index=key, columns="sample", values="Nvalid_cov")
    mod = agg.pivot_table(index=key, columns="sample", values="Nmod")
    samples = [s for s in cov.columns if s in ref_s or s in test_s]
    cov, mod = cov[samples], mod[samples]
    gidx = np.array([0 if s in ref_s else 1 for s in samples], dtype=int)

    # Require the coverage floor in at least min_samples_per_group samples WITHIN EACH group, NOT in
    # every sample. `cov.min(axis=1) >= min_cov` demanded the floor in all samples, so one low-coverage
    # library in either arm killed a site that was well-covered in the rest -- and the beta-binomial
    # already models per-sample coverage and drops n<=0, so the floor's job is only to guarantee each
    # group has enough covered libraries to estimate. (This is coverage-robustness; it does NOT report
    # isoform-switching-driven site changes -- that stays with the isoform-usage test by design.)
    _covn = cov.to_numpy(dtype=float)
    _ok = (~np.isnan(_covn)) & (_covn >= args.min_cov)          # (site x sample) meets the floor
    _nref = _ok[:, gidx == 0].sum(axis=1)
    _ntest = _ok[:, gidx == 1].sum(axis=1)
    keep = pd.Series((_nref >= args.min_samples_per_group) & (_ntest >= args.min_samples_per_group),
                     index=cov.index)
    cov, mod = cov[keep], mod[keep]
    unit = "site x transcript" if by_tx else "site"
    if args.verbose:
        print(f"[condition_mod] {name}: {len(cov):,} {unit}s with >={args.min_cov}x in >="
              f"{args.min_samples_per_group} samples/group ({len(ref_s)} {args.reference} vs "
              f"{len(test_s)} {args.test})", flush=True)
    if cov.empty:
        pd.DataFrame(columns=out_cols).to_csv(args.out_tsv, sep="\t", index=False)
        return

    # A site may now include samples uncovered at this (site, ZN) (NaN) -- the per-group filter above
    # only guarantees ENOUGH covered samples per group. Map NaN -> 0 so beta_binomial_diff drops those
    # as n<=0 (rather than propagating NaN into the likelihood).
    K = np.nan_to_num(mod.to_numpy(dtype=float), nan=0.0)
    N = np.nan_to_num(cov.to_numpy(dtype=float), nan=0.0)
    K = np.minimum(K, N)   # defensive: Nmod>Nvalid_cov (upstream glitch) otherwise yields a garbage significant p
    ref_names = [s for s in samples if s in ref_s]
    test_names = [s for s in samples if s in test_s]
    sites = [(i, K[i], N[i], gidx) for i in range(K.shape[0])]
    res = diffstats.beta_binomial_diff(sites, prior_weight=args.prior_weight,
                                       min_group_samples=args.min_samples_per_group,
                                       ref_df=args.ref_df, calibrate=False,
                                       site_weight=diffstats.parse_site_weight(args.site_weight))
    if not res:
        pd.DataFrame(columns=out_cols).to_csv(args.out_tsv, sep="\t", index=False)
        return

    idx = cov.index
    rows = []
    for r in res:
        i = r["key"]
        kv = dict(zip(key, idx[i] if isinstance(idx[i], tuple) else (idx[i],)))
        # per-replicate observed modified fractions, so the report can plot the replicate spread
        with np.errstate(divide="ignore", invalid="ignore"):
            mu_i = np.where(N[i] > 0, K[i] / N[i], np.nan)
        # keep only finite per-replicate values -> always valid JSON (a NaN would serialize to the
        # bare token NaN, which is not valid JSON and breaks the report's json.loads)
        per_rep = {"reference": {s: round(float(mu_i[j]), 4)
                                 for j, s in enumerate(samples) if s in ref_s and np.isfinite(mu_i[j])},
                   "test": {s: round(float(mu_i[j]), 4)
                            for j, s in enumerate(samples) if s in test_s and np.isfinite(mu_i[j])}}
        row = {
            "contrast": name, "chrom": kv["chrom"], "start0": kv["start0"], "end0": kv["end0"],
            "strand": kv["strand"], "mod_code": kv["mod_code"],
            "gene_name": (genes.loc[idx[i]] if genes is not None else ""),
            "n_reference": r["n_reference"], "n_test": r["n_test"],
            "reads_reference": r["reads_reference"], "reads_test": r["reads_test"],
            "mu_reference": round(r["mu_reference"], 5), "mu_test": round(r["mu_test"], 5),
            "delta": round(r["delta"], 5), "dispersion": round(r["dispersion"], 6),
            "lrt_stat": round(r["lrt_stat"], 4), "p_value": r["p_value"],
            "per_replicate_json": json.dumps(per_rep, separators=(",", ":")),
        }
        if by_tx:
            row["ZN_transcript_index"] = kv["ZN_transcript_index"]
        rows.append(row)
    out = pd.DataFrame(rows)
    out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
    out["_abs"] = out["delta"].abs()
    out = out.sort_values(["p_adj_bh", "_abs"], ascending=[True, False]).drop(columns="_abs")
    out = out[out_cols].reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    if args.verbose:
        sig = int((out["p_adj_bh"] < 0.05).sum())
        print(f"[condition_mod] {name}: {len(out):,} sites tested, {sig:,} at FDR<0.05 -> {args.out_tsv}", flush=True)


if __name__ == "__main__":
    main()
