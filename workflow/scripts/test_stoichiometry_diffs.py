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
from scipy.stats import chi2 as _chi2_dist, chi2_contingency, fisher_exact
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from plot_utils import save_figure
from genotype_utils import add_heterogeneity_flag, informative_strata, stratum_heterogeneity

# Canonical output schema, so the header is identical whether or not any site was testable (a
# zero-sites early-exit previously wrote a 14-column header while a real run emitted 21+).
RESULT_COLS = [
    "gene_name", "mod_code", "chrom", "start0", "end0", "strand",
    "n_tx_tested", "test_name", "stat_name", "stat_value", "p_value",
    "n_strata_informative", "strata_heterogeneous", "strata_heterogeneity_p", "strata_heterogeneity_p_adj",
    "effect_max_abs_frac_diff", "test_name_pooled", "stat_value_pooled", "p_value_pooled", "effect_max_abs_frac_diff_pooled",
    "per_transcript_json", "p_adj_bh",
]


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
    """Vectorized BH-FDR. Ranks only finite p-values (statsmodels semantics); a NaN would
    otherwise poison every adjusted p-value via the reversed minimum.accumulate. NaN stays NaN."""
    p = np.asarray(pvals, dtype=float)
    if p.size == 0:
        return p
    out = np.full(p.size, np.nan, dtype=float)
    idx = np.flatnonzero(np.isfinite(p))
    m = idx.size
    if m == 0:
        return out
    pf = p[idx]
    order = np.argsort(pf)
    ranks = np.empty(m, dtype=int)
    ranks[order] = np.arange(1, m + 1)
    adj = pf * m / ranks
    adj_sorted = np.minimum.accumulate(adj[order][::-1])[::-1]
    adj_final = np.empty(m, dtype=float)
    adj_final[order] = adj_sorted
    out[idx] = np.clip(adj_final, 0.0, 1.0)
    return out


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


def cmh_general_association(strata):
    """Generalized Cochran-Mantel-Haenszel GENERAL-ASSOCIATION test across a list of r x 2 strata
    (rows = transcripts in a FIXED order, cols = [Nmod, Nunmod]); one stratum per sample.

    BLOCKER-3: the previous test summed every sample's counts into ONE 2x2 and tested that. With
    per-sample rate heterogeneity AND per-sample transcript-composition imbalance, Simpson's paradox
    manufactures a large "effect" at astronomical significance where the within-sample difference is
    exactly zero (measured p=1.2e-23 on a constructed null). Stratifying by sample removes the confound.

    Reduces to the standard 2x2 CMH for r == 2. Returns (statistic, p_value, df, n_informative_strata);
    a stratum with no modified OR no unmodified calls, or <2 covered transcripts, carries no information
    and is dropped (this is the intended power cost -- a sample contributing to only one partition cannot
    separate transcript from sample). NaN p if no stratum is informative."""
    r = strata[0].shape[0]
    if r < 2:
        return float("nan"), float("nan"), 0, 0
    A = np.zeros(r - 1)
    V = np.zeros((r - 1, r - 1))
    used = 0
    for T in strata:
        T = np.asarray(T, dtype=float)
        R = T.sum(axis=1)            # per-transcript coverage in this stratum
        C = T.sum(axis=0)            # [total_mod, total_unmod]
        N = T.sum()
        if N < 2 or C[0] <= 0 or C[1] <= 0 or (R > 0).sum() < 2:
            continue
        a = T[:r - 1, 0]
        Rr = R[:r - 1]
        A += a - Rr * C[0] / N
        f = C[0] * C[1] / (N * N * (N - 1))
        V += f * (N * np.diag(Rr) - np.outer(Rr, Rr))
        used += 1
    if used == 0:
        return float("nan"), float("nan"), r - 1, 0
    try:
        Q = float(A @ np.linalg.solve(V, A))
    except np.linalg.LinAlgError:
        Q = float(A @ np.linalg.pinv(V) @ A)   # transcript with no across-stratum variation -> singular
    if not np.isfinite(Q) or Q < 0:
        return float("nan"), float("nan"), r - 1, used
    return Q, float(_chi2_dist.sf(Q, r - 1)), r - 1, used


def mh_max_abs_rate_diff(strata):
    """Mantel-Haenszel coverage-weighted rate difference, max over transcript pairs -- the effect size
    consistent with the stratified test (weights w_k = R_i R_j / N_k, only strata covering both)."""
    r = strata[0].shape[0]
    best = 0.0
    for i in range(r):
        for j in range(i + 1, r):
            num = den = 0.0
            for T in strata:
                T = np.asarray(T, dtype=float)
                Ri, Rj, N = T[i].sum(), T[j].sum(), T.sum()
                if N <= 0 or Ri <= 0 or Rj <= 0:
                    continue
                w = Ri * Rj / N
                num += w * (T[i, 0] / Ri - T[j, 0] / Rj)
                den += w
            if den > 0:
                best = max(best, abs(num / den))
    return best


def summarize_site(df_site, min_cov, which_test, pseudocount, alternative):
    """
    Collapse per-sample rows to per-transcript totals and run the appropriate test.
    Returns dict with stats or None if <2 transcripts pass coverage.
    """
    # per-transcript totals
    grp = (df_site.groupby("ZN_transcript_index", as_index=False)[["Nvalid_cov", "Nmod"]]
                 .sum())
    # clip Nmod to coverage: an upstream Nmod > Nvalid_cov glitch would otherwise make Nunmod negative
    # and blow up the contingency test (ValueError) -- the sibling test_condition_mod_diffs guards this too.
    grp["Nmod"] = grp[["Nmod", "Nvalid_cov"]].min(axis=1)
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
        # single zero cell -> odds genuinely infinite (keep inf); nan margins handled below.
        odds = float(odds) if np.isfinite(odds) else float("inf")
        return "fisher_exact_2x2", "fisher_odds", odds, float(p)

    def do_chi2_rx2(tab, pc):
        tab_pc = tab + pc
        chi2, p, dof, _ = chi2_contingency(tab_pc, correction=False)
        return f"chi2_{tab.shape[0]}x{tab.shape[1]}_pc{pc:g}", "chi2", float(chi2), float(p)

    def do_mc_exact_rx2(tab, n_resamples=9999, seed=12345):
        """Monte-Carlo EXACT test for an r x 2 table (cols = [modified, unmodified]).

        M2: the single-informative-stratum PRIMARY previously used an asymptotic chi2 here, which is
        ~2x anti-conservative on the realistic m6A regime -- sparse, unequal coverage, low rate (null
        P(p<0.05) = 0.103 at per-transcript depths [1,200,1], rate 0.15). This is R's
        chisq.test(simulate.p.value=TRUE): resample tables with the SAME margins under independence
        (multivariate-hypergeometric: distribute the total modified reads across transcripts by their
        coverage) and compare chi2 statistics. Deterministic (fixed seed) so the output stays
        byte-reproducible; fully vectorised. Falls back to asymptotic chi2 only if resampling is
        impossible (a degenerate margin, already screened as untestable upstream)."""
        tab = np.asarray(tab, dtype=float)
        row_tot = tab.sum(axis=1)                     # per-transcript coverage (fixed row margins)
        n_mod = float(tab[:, 0].sum()); N = float(tab.sum()); n_unmod = N - n_mod
        if n_mod <= 0 or n_unmod <= 0 or (row_tot <= 0).any():
            return do_chi2_rx2(tab, pseudocount)
        exp_mod = row_tot * n_mod / N
        exp_unmod = row_tot * n_unmod / N
        obs = float(((tab[:, 0] - exp_mod) ** 2 / exp_mod + (tab[:, 1] - exp_unmod) ** 2 / exp_unmod).sum())
        rng = np.random.default_rng(seed)
        draws = rng.multivariate_hypergeometric(row_tot.astype(int), int(n_mod), size=n_resamples)
        mod = draws.astype(float); unmod = row_tot[None, :] - mod
        chi2s = ((mod - exp_mod) ** 2 / exp_mod + (unmod - exp_unmod) ** 2 / exp_unmod).sum(axis=1)
        p = (int(np.sum(chi2s >= obs - 1e-9)) + 1) / (n_resamples + 1)
        return f"montecarlo_exact_{tab.shape[0]}x2", "chi2", obs, float(p)

    # A zero marginal (all transcripts 100%- or 0%-modified, i.e. no variation in the modified column)
    # is UNTESTABLE: fisher returns a nan odds ratio and chi2_contingency RAISES at pseudocount=0.
    # Set NaN p (excluded from the BH family -- it can never be significant) rather than crashing the
    # stage or writing a nan/inf statistic that later coerces to a misleading value. The effect size and
    # per-transcript block below are still computed so the row is otherwise complete.
    if not ((table.sum(axis=0) > 0).all() and (table.sum(axis=1) > 0).all()):
        used_test, stat_name, stat_value, pval = "untestable", "none", float("nan"), float("nan")
    elif which_test == "auto":
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

    # pooled effect size: max absolute Δ among POOLED stoichiometries (kept for transparency)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = grp_f["Nmod"] / grp_f["Nvalid_cov"].replace(0, np.nan)
    frac = frac.fillna(0.0).to_numpy()
    max_diff_pooled = 0.0
    for i in range(len(frac)):
        for j in range(i + 1, len(frac)):
            max_diff_pooled = max(max_diff_pooled, float(abs(frac[i] - frac[j])))

    # PRIMARY test (BLOCKER-3): stratify by SAMPLE. Build one r x 2 table per sample over the tested
    # transcripts (fixed row order; a transcript uncovered in a sample is a zero row) and run the
    # generalized CMH. This replaces the sample-pooled Fisher/chi2 as the reported statistic; the pooled
    # values are retained as *_pooled so the confound is visible rather than hidden.
    tested_zn = [int(z) for z in grp_f["ZN_transcript_index"].tolist()]
    zn_pos = {z: i for i, z in enumerate(tested_zn)}
    strata = []
    if "sample" in df_site.columns:
        sample_iter = df_site.groupby("sample")
    else:
        sample_iter = [("_all", df_site)]
    for _s, sub in sample_iter:
        by_zn = sub.groupby("ZN_transcript_index")[["Nvalid_cov", "Nmod"]].sum()
        T = np.zeros((len(tested_zn), 2))
        for z, row in by_zn.iterrows():
            if int(z) in zn_pos:
                cov = float(row["Nvalid_cov"]); mod = min(float(row["Nmod"]), cov)
                T[zn_pos[int(z)]] = (mod, cov - mod)
        strata.append(T)
    # PRIMARY from the INFORMATIVE strata only: >=2 -> generalized CMH; exactly 1 -> the EXACT test on
    # THAT stratum alone; 0 -> NaN (leaves the BH family). The old <2 fallback used the fully-POOLED
    # table -- but that pools the non-informative samples' reads back in, which is the very Simpson
    # statistic the stratification removed (e.g. one exactly-independent informative sample + two
    # constant samples -> pooled p=1.5e-126, effect 0.72, where the honest answer is p=1.0, effect 0).
    inf = informative_strata(strata)
    n_strata = len(inf)
    if n_strata >= 2:
        cmh_stat, cmh_p, _cmh_df, _ = cmh_general_association(inf)
        primary_test = "cmh_2x2" if len(tested_zn) == 2 else f"cmh_general_{len(tested_zn)}x2"
        primary_stat_name, primary_stat, primary_p = "cmh_chi2", cmh_stat, cmh_p
        primary_eff = mh_max_abs_rate_diff(inf)
    elif n_strata == 1:
        T1 = inf[0]
        # 2x2 -> Fisher (exact); r x 2 -> Monte-Carlo EXACT (M2: the old asymptotic chi2 here was ~2x
        # anti-conservative on sparse/unequal/low-rate strata, the realistic m6A regime).
        primary_test, primary_stat_name, primary_stat, primary_p = (
            do_fisher_2x2(T1) if T1.shape[0] == 2 else do_mc_exact_rx2(T1))
        primary_eff = mh_max_abs_rate_diff(inf)
    else:
        primary_test, primary_stat_name = "untestable", "none"
        primary_stat, primary_p, primary_eff = float("nan"), float("nan"), float("nan")
    # count-aware heterogeneity across the informative strata (generalizes to this r x 2 headline test)
    _hstat, het_p, _hdf, _ = stratum_heterogeneity(inf)
    strata_heterogeneous = bool(np.isfinite(het_p) and het_p < 0.05)

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
        # primary = CMH (>=2 informative strata) / exact on the single informative stratum / NaN
        "test_name": primary_test,
        "stat_name": primary_stat_name,
        "stat_value": primary_stat,
        "p_value": primary_p,
        "n_strata_informative": int(n_strata),
        "strata_heterogeneous": bool(strata_heterogeneous),
        "strata_heterogeneity_p": round(het_p, 6) if np.isfinite(het_p) else float("nan"),
        "effect_max_abs_frac_diff": round(primary_eff, 6),
        # pooled companions (the OLD sample-pooled statistic -- do not rank on these)
        "test_name_pooled": used_test,
        "stat_value_pooled": stat_value,
        "p_value_pooled": pval,
        "effect_max_abs_frac_diff_pooled": round(max_diff_pooled, 6),
        "per_transcript": per_tx,
    }


def make_plot(df_site, per_tx, title, out_png):
    """
    Two-panel figure:
      (1) per-sample stoichiometries by transcript
      (2) pooled coverages with modified-call overlay
    """
    samp = (df_site.groupby(["ZN_transcript_index", "sample"], as_index=False)[["Nvalid_cov", "Nmod"]]
                  .sum())
    samp = samp.loc[samp["Nvalid_cov"] > 0].copy()
    samp["frac"] = samp["Nmod"] / samp["Nvalid_cov"]

    pooled = pd.DataFrame(per_tx).sort_values("ZN").reset_index(drop=True)
    zn_order = sorted(pooled["ZN"].tolist())
    sample_names = sorted(samp["sample"].astype(str).unique())
    offsets = np.linspace(-0.18, 0.18, num=max(1, len(sample_names)))
    palette = plt.get_cmap("tab10")
    sample_colors = {
        sample: palette(idx % palette.N)
        for idx, sample in enumerate(sample_names)
    }

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(8.5, 7.8),
        sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1.0]},
        layout="constrained",   # reflows at draw time -> absorbs the enlarged house-style fonts
    )

    handles = []
    for offset, sample in zip(offsets, sample_names):
        sub = samp.loc[samp["sample"].astype(str) == sample]
        if sub.empty:
            continue
        handle = ax1.scatter(
            sub["ZN_transcript_index"].astype(float) + offset,
            sub["frac"],
            s=48,
            color=sample_colors[sample],
            edgecolors="white",
            linewidths=0.7,
            alpha=0.92,
            label=sample,
            zorder=3,
        )
        handles.append(handle)
    ax1.set_ylabel("Modified fraction")
    ax1.set_title("Per-sample stoichiometries")
    ax1.set_ylim(0.0, 1.02)
    ax1.grid(True, axis="y", linestyle="--", linewidth=0.6, alpha=0.4)

    ax2.bar(
        pooled["ZN"],
        pooled["Ncov"],
        width=0.72,
        color="#d9e1ea",
        edgecolor="#6a7785",
        linewidth=1.0,
        label="Total coverage",
        zorder=1,
    )
    ax2.bar(
        pooled["ZN"],
        pooled["Nmod"],
        width=0.72,
        color="#d46a5d",
        edgecolor="#983f33",
        linewidth=1.0,
        label="Modified calls",
        zorder=2,
    )
    max_cov = float(pooled["Ncov"].max()) if not pooled.empty else 0.0
    label_pad = max(1.0, max_cov * 0.02)
    for row in pooled.itertuples(index=False):
        frac_pct = 100.0 * float(row.frac)
        ax2.text(
            row.ZN,
            float(row.Ncov) + label_pad,
            f"{frac_pct:.0f}%",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#4a3f35",
        )
    ax2.set_xlabel("Transcript index (ZN)")
    ax2.set_ylabel("Pooled coverage")
    ax2.set_title("Pooled coverage with modified-call overlay")
    ax2.set_ylim(0.0, max(1.0, max_cov * 1.18))
    ax2.grid(True, axis="y", linestyle="--", linewidth=0.6, alpha=0.4)

    for ax in (ax1, ax2):
        ax.set_xticks(zn_order)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#5e554c")
        ax.spines["bottom"].set_color("#5e554c")
        ax.tick_params(colors="#3f3933")

    if handles:
        legend_cols = min(3, len(handles))
        # 'outside lower center' -> constrained_layout reserves a band for the legend below the
        # axes, so it never overlaps the x-axis label (even at the enlarged house font sizes).
        fig.legend(
            handles=handles,
            labels=sample_names,
            loc="outside lower center",
            ncol=legend_cols,
            frameon=False,
            fontsize=8,
            handletextpad=0.4,
            columnspacing=1.2,
        )

    fig.suptitle(title, fontsize=12)
    save_figure(fig, out_png, dpi=300, bbox_inches="tight")   # PNG + PDF + SVG
    plt.close(fig)


# ----------------------------- Main ------------------------------

def main():
    args = parse_args()

    if args.alternative != "two-sided":
        # A one-sided alternative reaches only sites tested by the single-stratum EXACT test
        # (do_fisher_2x2); the multi-stratum primary is the generalized CMH, which is inherently
        # two-sided. So a directional --alternative silently mixes one- and two-sided p-values in one
        # BH family. Surface it rather than let the family be quietly inconsistent.
        print(f"[test_diffs] WARNING: --alternative={args.alternative} is one-sided, but the multi-stratum "
              f"primary test (generalized CMH) is inherently TWO-SIDED. The directional alternative applies "
              f"ONLY to single-informative-stratum sites (the exact test); CMH sites stay two-sided, so the "
              f"BH family mixes sidedness. Use --alternative two-sided for a homogeneous family.",
              file=sys.stderr, flush=True)

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

    # An empty input (0 data rows after filtering -- e.g. a gene set where no site passed the
    # aggregation filter, or an aggressive gene/mod filter) must not reach df.apply below:
    # on an empty frame df.apply(axis=1) returns an empty DataFrame (not a Series), which
    # crashes the site_key assignment ("Cannot set a DataFrame with multiple columns...").
    # Emit the empty result table + figs dir and exit 0 so the pipeline continues.
    if df.empty:
        out_tsv = f"{args.out_prefix}__ZN_site_diff_results.tsv"
        os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
        pd.DataFrame(columns=RESULT_COLS).to_csv(out_tsv, sep="\t", index=False)
        os.makedirs(f"{args.out_prefix}__figs", exist_ok=True)
        print(f"[info] --in-tsv has no usable data rows; wrote empty {out_tsv} and continuing.")
        return

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
            "n_strata_informative": res["n_strata_informative"],
            "strata_heterogeneous": res["strata_heterogeneous"],
            "strata_heterogeneity_p": res["strata_heterogeneity_p"],
            "effect_max_abs_frac_diff": res["effect_max_abs_frac_diff"],
            "test_name_pooled": res["test_name_pooled"],
            "stat_value_pooled": res["stat_value_pooled"],
            "p_value_pooled": res["p_value_pooled"],
            "effect_max_abs_frac_diff_pooled": res["effect_max_abs_frac_diff_pooled"],
            "per_transcript_json": json.dumps(res["per_transcript"], separators=(",", ":")),
        })

    if not results:
        # Single-isoform loci (e.g. mitochondrial chrM genes) have no site with >=2
        # transcript partitions to contrast -> nothing to test. Emit an empty results
        # table (header only) + empty figs dir and exit 0 so the pipeline continues to
        # genotype/report instead of aborting.
        out_tsv = f"{args.out_prefix}__ZN_site_diff_results.tsv"
        os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
        pd.DataFrame(columns=RESULT_COLS).to_csv(out_tsv, sep="\t", index=False)
        os.makedirs(f"{args.out_prefix}__figs", exist_ok=True)
        print("[info] No sites had >=2 transcripts meeting the coverage threshold; "
              f"wrote empty {out_tsv} and continuing (nothing to test).")
        return

    res_df = pd.DataFrame(results)
    res_df["p_adj_bh"] = benjamini_hochberg(res_df["p_value"].values)
    res_df = add_heterogeneity_flag(res_df)        # BH-adjust the heterogeneity flag like every other p
    res_df = res_df.sort_values(["p_adj_bh", "effect_max_abs_frac_diff"], ascending=[True, False]).reset_index(drop=True)
    for _c in RESULT_COLS:                          # stable schema regardless of what was found
        if _c not in res_df.columns:
            res_df[_c] = pd.NA
    res_df = res_df[RESULT_COLS]

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
