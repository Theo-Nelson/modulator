#!/usr/bin/env python3
"""
Test for stoichiometry differences between transcripts (ZN) at each genomic site.

Input
-----
A TSV like the one you showed (one row per sample x transcript x site), with columns:
  gene_name, gene_id, mod_code, chrom, start0, end0, strand,
  ZN_transcript_index, sample, Nvalid_cov, Nmod, ... , frac_modified

What it does
------------
- For each site (chrom, start0, end0, strand, mod_code, gene_name):
    * Aggregate counts across samples within each transcript (ZN):
        - Ncov_zn = sum(Nvalid_cov), Nmod_zn = sum(Nmod)
    * Keep only transcripts with Ncov_zn >= --min-cov (default 20).
    * If >= 2 transcripts pass:
        - Build a contingency table with rows=ZN, cols=[modified, unmodified].
        - Do Fisher's exact test if 2x2, else Chi-square test of independence.
        - Record max pairwise |p_i - p_j| across transcripts (effect size).
- Multiple-testing correction (Benjamini-Hochberg FDR) across all tested sites.
- Output a ranked TSV and PNG plots for the top K sites.

Outputs
-------
1) <out_prefix>__ZN_site_diff_results.tsv
2) <out_prefix>__figs/ (PNG plots for top --topk sites)

Requirements
------------
pandas, numpy, scipy, matplotlib
"""

import os, sys, argparse, json
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd

from scipy.stats import chi2_contingency, fisher_exact
import matplotlib.pyplot as plt


def parse_args():
    ap = argparse.ArgumentParser(description="Stoichiometry differences between transcripts at each site")
    ap.add_argument("--in-tsv", required=True,
                    help="Input TSV with columns incl. gene_name, mod_code, chrom, start0, end0, strand, ZN_transcript_index, sample, Nvalid_cov, Nmod")
    ap.add_argument("--out-prefix", required=True, help="Output prefix for results and figures")
    ap.add_argument("--min-cov", type=int, default=20,
                    help="Minimum TOTAL coverage per transcript (summing samples) to include in a site's test (default: 20)")
    ap.add_argument("--topk", type=int, default=10, help="How many top sites to plot (default: 10)")
    ap.add_argument("--gene-filter", nargs="*", default=None,
                    help="Optional: restrict to these gene_name(s)")
    ap.add_argument("--mod-filter", nargs="*", default=None,
                    help="Optional: restrict to these mod_code(s), e.g. a m 17596 ...")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def benjamini_hochberg(p):
    """Return BH-adjusted p-values for a 1D array-like of p-values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n+1)
    adj = p * n / ranks
    # enforce monotonicity
    adj_sorted = np.minimum.accumulate(adj[order][::-1])[::-1]
    adj_bh = np.empty_like(adj_sorted)
    adj_bh[order] = adj_sorted
    return np.minimum(adj_bh, 1.0)


def site_key(row):
    # include gene_name in the site (helpful for downstream splitting),
    # but the test itself is at the exact genomic locus + mod_code.
    return (row["gene_name"], row["mod_code"], row["chrom"], int(row["start0"]), int(row["end0"]), row["strand"])


def summarize_site(df_site, min_cov):
    """
    df_site: rows for a single genomic site across samples and transcripts.
    Returns dict with test results or None if not enough transcripts pass.
    """
    # aggregate per transcript across samples
    grp = df_site.groupby("ZN_transcript_index", as_index=False)[["Nvalid_cov","Nmod"]].sum()
    grp["Nunmod"] = grp["Nvalid_cov"] - grp["Nmod"]
    grp = grp.sort_values("ZN_transcript_index")

    # filter by coverage
    grp_f = grp[grp["Nvalid_cov"] >= min_cov].copy()
    if len(grp_f) < 2:
        return None

    # contingency matrix (rows=ZN, cols=[mod, unmod])
    table = grp_f[["Nmod","Nunmod"]].values.astype(int)

    # test
    if table.shape == (2, 2):
        # Fisher exact (two-sided)
        odds, pval = fisher_exact(table, alternative="two-sided")
        stat_name = "fisher_odds"
        stat_val = odds
        test_name = "fisher_exact_2x2"
    else:
        chi2, pval, dof, _ = chi2_contingency(table, correction=False)
        stat_name = "chi2"
        stat_val = chi2
        test_name = f"chi2_{table.shape[0]}x2"

    # effect size: max absolute difference between transcript fractions (pooled across samples)
    grp_f["frac"] = grp_f["Nmod"] / grp_f["Nvalid_cov"].replace(0, np.nan)
    max_diff = 0.0
    for a, b in combinations(grp_f["ZN_transcript_index"].tolist(), 2):
        pa = float(grp_f.loc[grp_f["ZN_transcript_index"]==a, "frac"])
        pb = float(grp_f.loc[grp_f["ZN_transcript_index"]==b, "frac"])
        max_diff = max(max_diff, abs(pa - pb))

    # keep per-transcript stats for output
    per_tx = []
    for _, r in grp_f.iterrows():
        per_tx.append({
            "ZN": int(r["ZN_transcript_index"]),
            "Ncov": int(r["Nvalid_cov"]),
            "Nmod": int(r["Nmod"]),
            "frac": round(float(r["frac"]), 6)
        })

    return {
        "n_tx_tested": int(len(grp_f)),
        "test_name": test_name,
        "stat_name": stat_name,
        "stat_value": float(stat_val),
        "p_value": float(pval),
        "effect_max_abs_frac_diff": round(float(max_diff), 6),
        "per_transcript": per_tx,
    }


def make_plot(df_site, per_tx, title, out_png):
    """
    Two panels in one figure:
      A) per-sample dots for each transcript (x=ZN, y=frac per sample)
      B) bars for pooled fractions per transcript (using aggregated counts)
    """
    # per-sample fractions
    df_samp = (df_site.groupby(["ZN_transcript_index","sample"], as_index=False)
                     [["Nvalid_cov","Nmod"]].sum())
    df_samp = df_samp[df_samp["Nvalid_cov"] > 0].copy()
    df_samp["frac"] = df_samp["Nmod"] / df_samp["Nvalid_cov"]

    # pooled per transcript (from per_tx)
    pooled = pd.DataFrame(per_tx)

    zn_order = sorted(pooled["ZN"].tolist())
    # Figure with two axes (stacked vertically)
    fig = plt.figure(figsize=(7, 7))

    # Panel A: per-sample dots
    ax1 = fig.add_subplot(2,1,1)
    for zn in zn_order:
        sub = df_samp[df_samp["ZN_transcript_index"]==zn]
        ax1.scatter([zn]*len(sub), sub["frac"], label=f"ZN {zn}")
    ax1.set_xlabel("Transcript index (ZN)")
    ax1.set_ylabel("Stoichiometry (per-sample)")
    ax1.set_title("Per-sample stoichiometries")
    # Put a light grid
    ax1.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

    # Panel B: pooled bars
    ax2 = fig.add_subplot(2,1,2)
    ax2.bar(pooled["ZN"], pooled["frac"])
    ax2.set_xlabel("Transcript index (ZN)")
    ax2.set_ylabel("Stoichiometry (pooled)")
    ax2.set_title("Pooled stoichiometries (across samples)")
    ax2.set_ylim(0, 1.0)
    ax2.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

    fig.suptitle(title, y=0.98, fontsize=12)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def main():
    args = parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t", low_memory=False)
    needed = {"gene_name","mod_code","chrom","start0","end0","strand",
              "ZN_transcript_index","sample","Nvalid_cov","Nmod"}
    missing = needed - set(df.columns)
    if missing:
        sys.exit(f"Missing required columns: {sorted(missing)}")

    # Filters
    if args.gene_filter:
        df = df[df["gene_name"].astype(str).isin(args.gene_filter)]
    if args.mod_filter:
        df = df[df["mod_code"].astype(str).isin(args.mod_filter)]

    # Ensure integer/numeric types
    for c in ["start0","end0","Nvalid_cov","Nmod","ZN_transcript_index"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["start0","end0","Nvalid_cov","Nmod","ZN_transcript_index"]).copy()
    df["start0"] = df["start0"].astype(int)
    df["end0"] = df["end0"].astype(int)
    df["ZN_transcript_index"] = df["ZN_transcript_index"].astype(int)

    # Define site key
    df["site_key"] = df.apply(site_key, axis=1)

    results = []
    site_groups = df.groupby("site_key", sort=False)

    if args.verbose:
        print(f"[info] evaluating {len(site_groups)} sites with min_cov={args.min_cov}", file=sys.stderr)

    for sk, df_site in site_groups:
        res = summarize_site(df_site, min_cov=args.min_cov)
        if res is None:
            continue
        gene_name, mod_code, chrom, start0, end0, strand = sk
        results.append({
            "gene_name": gene_name,
            "mod_code": mod_code,
            "chrom": chrom,
            "start0": int(start0),
            "end0": int(end0),
            "strand": strand,
            "n_tx_tested": res["n_tx_tested"],
            "test_name": res["test_name"],
            "stat_name": res["stat_name"],
            "stat_value": res["stat_value"],
            "p_value": res["p_value"],
            "effect_max_abs_frac_diff": res["effect_max_abs_frac_diff"],
            "per_transcript_json": json.dumps(res["per_transcript"], separators=(",",":")),
        })

    if not results:
        sys.exit("No sites had ≥2 transcripts meeting the coverage threshold; nothing to test.")

    res_df = pd.DataFrame(results)
    res_df = res_df.sort_values("p_value").reset_index(drop=True)
    res_df["p_adj_bh"] = benjamini_hochberg(res_df["p_value"].values)
    res_df = res_df.sort_values(["p_adj_bh","effect_max_abs_frac_diff"], ascending=[True, False])

    out_tsv = f"{args.out_prefix}__ZN_site_diff_results.tsv"
    res_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[ok] wrote {out_tsv}  (n_sites_tested={len(res_df)})")

    # Make plots for top K
    figs_dir = f"{args.out_prefix}__figs"
    os.makedirs(figs_dir, exist_ok=True)
    topk = min(args.topk, len(res_df))

    # build a lookup DataFrame keyed by site for plotting
    df_lookup = df.set_index(["gene_name","mod_code","chrom","start0","end0","strand"])

    for i in range(topk):
        r = res_df.iloc[i]
        key = (r["gene_name"], r["mod_code"], r["chrom"], int(r["start0"]), int(r["end0"]), r["strand"])
        df_site = df_lookup.loc[[key]].reset_index()
        per_tx = json.loads(r["per_transcript_json"])

        title = (f"{r['gene_name']} | {r['mod_code']} | {r['chrom']}:{r['start0']}-{r['end0']}({r['strand']})\n"
                 f"{r['test_name']} p={r['p_value']:.2e}, FDR={r['p_adj_bh']:.2e}, "
                 f"max|Δfrac|={r['effect_max_abs_frac_diff']:.3f}")
        out_png = os.path.join(figs_dir,
                               f"site_{i+1:02d}__{r['gene_name']}__{r['mod_code']}__{r['chrom']}_{r['start0']}_{r['end0']}_{r['strand']}.png")
        make_plot(df_site, per_tx, title, out_png)

    print(f"[ok] saved {topk} figure(s) under {figs_dir}")


if __name__ == "__main__":
    main()

