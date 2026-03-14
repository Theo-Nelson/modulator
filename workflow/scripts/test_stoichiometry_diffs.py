#!/usr/bin/env python3
"""
Test stoichiometry differences between transcripts (ZN) at each genomic site.

Inputs
------
A "long" TSV produced by the ZN aggregator (RAW or FILTERED) with at least:
  gene_name, mod_code, chrom, start0, end0, strand,
  ZN_transcript_index, sample, Nvalid_cov, Nmod

What it does
------------
• Collapses per-sample rows to per-transcript totals at each site
• Keeps only transcripts with total coverage >= --min-cov
• If >=2 transcripts remain:
    - test='auto'   : Fisher's exact for 2×2; otherwise chi-square with pseudocount
    - test='fisher' : Fisher for 2×2; otherwise falls back to chi-square (warn)
    - test='chi2'   : chi-square for r×2 with pseudocount
• Computes an effect size: max absolute difference among pooled stoichiometries
• Benjamini–Hochberg FDR correction across sites
• Writes a results table and plots top-K sites

Notes
-----
• Robust to accidental "null"/"None"/"NA" strings provided to filters
• Uses matplotlib only; no seaborn; single-figure per site with 2 panels
"""

import os
import sys
import json
import argparse
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
import matplotlib.pyplot as plt


# ------------------------------ CLI ------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Test stoichiometry differences between transcripts (ZN) at each site"
    )
    ap.add_argument(
        "--in-tsv", required=True,
        help="Input long TSV (from ZN aggregation)."
    )
    ap.add_argument(
        "--out-prefix", required=True,
        help="Output prefix (results TSV and figures dir)."
    )
    ap.add_argument(
        "--min-cov", type=int, default=20,
        help="Minimum TOTAL coverage per transcript at a site (sum over samples). Default: 20"
    )
    ap.add_argument(
        "--topk", type=int, default=10,
        help="Number of top sites (by FDR then effect size) to plot. Default: 10"
    )
    ap.add_argument(
        "--test", choices=["auto", "fisher", "chi2"], default="auto",
        help='Hypothesis test: "auto" (2×2 Fisher else chi-square), "fisher", or "chi2".'
    )
    ap.add_argument(
        "--pseudocount", type=float, default=0.5,
        help="Pseudocount added to each cell for chi-square to stabilize expected values. Default: 0.5"
    )
    ap.add_argument(
        "--alternative", choices=["two-sided", "greater", "less"], default="two-sided",
        help='Alternative for Fisher’s exact (2×2 only). Default: "two-sided"'
    )
    ap.add_argument(
        "--gene-filter", nargs="*", default=None,
        help='Optional subset of gene_name, e.g. --gene-filter ALCAM RIOK3'
    )
    ap.add_argument(
        "--mod-filter", nargs="*", default=None,
        help='Optional subset of mod_code, e.g. --mod-filter a m'
    )
    ap.add_argument(
        "--verbose", action="store_true",
        help="Verbose logging to stderr."
    )
    return ap.parse_args()


# --------------------------- Utilities ---------------------------

NULL_TOKENS = {"", "null", "none", "na", "nil", "Null", "None", "NA", "NIL"}

def is_nullish(x) -> bool:
    if x is None:
        return True
    if isinstance(x, str) and x.strip() in NULL_TOKENS:
        return True
    return False


def benjamini_hochberg(pvals):
    """Vectorized BH-FDR."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    adj = p * n / ranks
    adj_sorted = np.minimum.accumulate(adj[order][::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = adj_sorted
    return np.minimum(out, 1.0)


def site_key_tuple(row):
    """Stable tuple key for grouping/lookup."""
    return (
        str(row["gene_name"]),
        str(row["mod_code"]),
        str(row["chrom"]),
        int(row["start0"]),
        int(row["end0"]),
        str(row["strand"]),
    )


def site_key_str_from_tuple(tup):
    g, m, chrom, s0, e0, st = tup
    return f"{g}|{m}|{chrom}:{int(s0)}-{int(e0)}({st})"


def summarize_site(df_site, min_cov, which_test, pseudocount, alternative):
    """
    Collapse per-sample rows to per-transcript totals and run the appropriate test.
    Returns dict with stats or None if <2 transcripts pass coverage.
    """
    # per-transcript totals
    grp = (df_site.groupby("ZN_transcript_index", as_index=False)[["Nvalid_cov", "Nmod"]]
                 .sum())
    grp["Nunmod"] = grp["Nvalid_cov"] - grp["Nmod"]
    grp = grp.sort_values("ZN_transcript_index")

    # coverage filter
    grp_f = grp.loc[grp["Nvalid_cov"] >= min_cov].copy()
    if len(grp_f) < 2:
        return None

    # build contingency (rows=transcripts; cols=[mod, unmod])
    table = grp_f[["Nmod", "Nunmod"]].to_numpy(dtype=float)

    # choose test
    nrows = table.shape[0]
    used_test = which_test
    stat_name = None
    stat_value = None
    pval = None

    def do_fisher_2x2(tab):
        odds, p = fisher_exact(tab.astype(int), alternative=alternative)
        return "fisher_exact_2x2", "fisher_odds", float(odds), float(p)

    def do_chi2_rx2(tab, pc):
        tab_pc = tab + pc
        chi2, p, dof, _ = chi2_contingency(tab_pc, correction=False)
        return f"chi2_{tab.shape[0]}x{tab.shape[1]}_pc{pc:g}", "chi2", float(chi2), float(p)

    if which_test == "auto":
        if nrows == 2:
            used_test, stat_name, stat_value, pval = do_fisher_2x2(table)
        else:
            used_test, stat_name, stat_value, pval = do_chi2_rx2(table, pseudocount)
    elif which_test == "fisher":
        if nrows == 2:
            used_test, stat_name, stat_value, pval = do_fisher_2x2(table)
        else:
            # graceful fallback
            used_test, stat_name, stat_value, pval = do_chi2_rx2(table, pseudocount)
    else:  # "chi2"
        used_test, stat_name, stat_value, pval = do_chi2_rx2(table, pseudocount)

    # effect size: max absolute Δ among pooled stoichiometries
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = grp_f["Nmod"] / grp_f["Nvalid_cov"].replace(0, np.nan)
    frac = frac.fillna(0.0).to_numpy()
    max_diff = 0.0
    for i in range(len(frac)):
        for j in range(i + 1, len(frac)):
            max_diff = max(max_diff, float(abs(frac[i] - frac[j])))

    # serialize per-transcript for convenience
    per_tx = []
    for _, r in grp_f.iterrows():
        ncv = int(r["Nvalid_cov"])
        nmd = int(r["Nmod"])
        per_tx.append({
            "ZN": int(r["ZN_transcript_index"]),
            "Ncov": ncv,
            "Nmod": nmd,
            "frac": 0.0 if ncv == 0 else float(nmd / ncv),
        })

    return {
        "n_tx_tested": int(len(grp_f)),
        "test_name": used_test,
        "stat_name": stat_name,
        "stat_value": stat_value,
        "p_value": pval,
        "effect_max_abs_frac_diff": round(max_diff, 6),
        "per_transcript": per_tx,
    }


def make_plot(df_site, per_tx, title, out_png):
    """
    Two-panel figure:
      (1) per-sample stoichiometries by transcript
      (2) pooled stoichiometry per transcript
    """
    # per-sample
    samp = (df_site.groupby(["ZN_transcript_index", "sample"], as_index=False)[["Nvalid_cov", "Nmod"]]
                  .sum())
    samp = samp.loc[samp["Nvalid_cov"] > 0].copy()
    samp["frac"] = samp["Nmod"] / samp["Nvalid_cov"]

    pooled = pd.DataFrame(per_tx)
    zn_order = sorted(pooled["ZN"].tolist())

    fig = plt.figure(figsize=(7, 7))

    ax1 = fig.add_subplot(2, 1, 1)
    for zn in zn_order:
        sub = samp.loc[samp["ZN_transcript_index"] == zn]
        if len(sub):
            ax1.scatter([zn] * len(sub), sub["frac"])
    ax1.set_xlabel("Transcript index (ZN)")
    ax1.set_ylabel("Stoichiometry (per-sample)")
    ax1.set_title("Per-sample stoichiometries")
    ax1.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

    ax2 = fig.add_subplot(2, 1, 2)
    ax2.bar(pooled["ZN"], pooled["frac"])
    ax2.set_xlabel("Transcript index (ZN)")
    ax2.set_ylabel("Stoichiometry (pooled)")
    ax2.set_title("Pooled stoichiometries (across samples)")
    ax2.set_ylim(0, 1.0)
    ax2.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

    fig.suptitle(title, y=0.98, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


# ----------------------------- Main ------------------------------

def main():
    args = parse_args()

    # Load
    df = pd.read_csv(args.in_tsv, sep="\t", low_memory=False)

    required = {
        "gene_name", "mod_code", "chrom", "start0", "end0", "strand",
        "ZN_transcript_index", "sample", "Nvalid_cov", "Nmod",
    }
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"Missing columns in --in-tsv: {sorted(missing)}")

    # Optional filters (tolerate nullish strings)
    if args.gene_filter and not (len(args.gene_filter) == 1 and is_nullish(args.gene_filter[0])):
        df = df[df["gene_name"].astype(str).isin([str(g) for g in args.gene_filter])]
    if args.mod_filter and not (len(args.mod_filter) == 1 and is_nullish(args.mod_filter[0])):
        df = df[df["mod_code"].astype(str).isin([str(m) for m in args.mod_filter])]

    # Coerce numeric columns
    for c in ["start0", "end0", "Nvalid_cov", "Nmod", "ZN_transcript_index"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["start0", "end0", "Nvalid_cov", "Nmod", "ZN_transcript_index"]).copy()
    df["start0"] = df["start0"].astype(int)
    df["end0"] = df["end0"].astype(int)
    df["ZN_transcript_index"] = df["ZN_transcript_index"].astype(int)

    # Build keys (fix for "Cannot set a DataFrame with multiple columns...")
    df["site_key"] = df.apply(site_key_tuple, axis=1)
    df["site_key_str"] = df["site_key"].map(site_key_str_from_tuple)

    # Group by site
    site_groups = df.groupby("site_key", sort=False)
    if args.verbose:
        print(f"[info] evaluating {len(site_groups)} sites with min_cov={args.min_cov}, test={args.test}", file=sys.stderr)

    results = []
    for sk, df_site in site_groups:
        res = summarize_site(
            df_site,
            min_cov=args.min_cov,
            which_test=args.test,
            pseudocount=args.pseudocount,
            alternative=args.alternative,
        )
        if res is None:
            continue

        g, m, chrom, s0, e0, st = sk
        results.append({
            "gene_name": g,
            "mod_code": m,
            "chrom": chrom,
            "start0": int(s0),
            "end0": int(e0),
            "strand": st,
            "n_tx_tested": res["n_tx_tested"],
            "test_name": res["test_name"],
            "stat_name": res["stat_name"],
            "stat_value": res["stat_value"],
            "p_value": res["p_value"],
            "effect_max_abs_frac_diff": res["effect_max_abs_frac_diff"],
            "per_transcript_json": json.dumps(res["per_transcript"], separators=(",", ":")),
        })

    if not results:
        sys.exit("No sites had ≥2 transcripts meeting the coverage threshold; nothing to test.")

    res_df = pd.DataFrame(results)
    res_df["p_adj_bh"] = benjamini_hochberg(res_df["p_value"].values)
    res_df = res_df.sort_values(["p_adj_bh", "effect_max_abs_frac_diff"], ascending=[True, False]).reset_index(drop=True)

    # Write results
    out_tsv = f"{args.out_prefix}__ZN_site_diff_results.tsv"
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    res_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[ok] wrote {out_tsv}  (n_sites_tested={len(res_df)})")

    # Figures
    figs_dir = f"{args.out_prefix}__figs"
    os.makedirs(figs_dir, exist_ok=True)
    topk = int(min(args.topk, len(res_df)))

    # Fast index by tuple key for re-slice
    df_site_index = df.set_index("site_key")

    for i in range(topk):
        r = res_df.iloc[i]
        key_tuple = (
            str(r["gene_name"]), str(r["mod_code"]), str(r["chrom"]),
            int(r["start0"]), int(r["end0"]), str(r["strand"])
        )
        try:
            df_site = df_site_index.loc[[key_tuple]].reset_index(drop=True)
        except KeyError:
            # should not happen, but be defensive
            continue

        per_tx = json.loads(r["per_transcript_json"])
        title = (
            f"{r['gene_name']} | {r['mod_code']} | "
            f"{r['chrom']}:{int(r['start0'])}-{int(r['end0'])}({r['strand']})\n"
            f"{r['test_name']} p={r['p_value']:.2e}, FDR={r['p_adj_bh']:.2e}, "
            f"max|Δfrac|={r['effect_max_abs_frac_diff']:.3f}"
        )
        out_png = os.path.join(
            figs_dir,
            f"site_{i+1:02d}__{r['gene_name']}__{r['mod_code']}__{r['chrom']}_{int(r['start0'])}_{int(r['end0'])}_{r['strand']}.png"
        )
        make_plot(df_site, per_tx, title, out_png)

    print(f"[ok] saved {topk} figure(s) under {figs_dir}")


if __name__ == "__main__":
    main()

