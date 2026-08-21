#!/usr/bin/env python3
"""Poly(A) tail-length distributions per fragmentform, and differential tail length between the
fragmentforms of a gene.

Consumes the per-read tail table (build_read_polya_table.py). Emits two TSVs:
  {out_prefix}_polya_fragmentform.tsv   one row per fragmentform (ZT): tail-length distribution
  {out_prefix}_taillength_diffs.tsv     one row per metagene (gene group): are the tail-length
                                        distributions of its competing fragmentforms different?

The differential test mirrors test_stoichiometry_diffs.py's group -> filter -> test -> results-row ->
per_*_json -> BH -> topK shape, but tail length is continuous, so it uses a continuous test
(Mann-Whitney U for 2 fragmentforms, Kruskal-Wallis for >2) instead of a contingency test. Effect
size is the spread of per-fragmentform median tail lengths (nt). Empty-input safe (header-only out).
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

from genotype_utils import benjamini_hochberg
from plot_utils import save_figure

FRAG_COLS = ["ZT", "gene_name", "metagene_index", "classification", "n_reads", "n_samples",
             "median_tail", "mean_tail", "std_tail", "q25_tail", "q75_tail", "min_tail", "max_tail"]
DIFF_COLS = ["metagene_index", "gene_name", "n_fragmentforms_tested", "n_reads", "test_name",
             "stat_name", "stat_value", "p_value", "effect_median_range_nt",
             "min_median_tail", "max_median_tail", "per_fragmentform_json", "p_adj_bh"]


def parse_args():
    ap = argparse.ArgumentParser(description="Per-fragmentform poly(A) tail distributions + between-fragmentform differential tail length.")
    ap.add_argument("--in-tsv", required=True, help="Per-read tail table from build_read_polya_table.py")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-reads", type=int, default=10, help="Min reads for a fragmentform to enter the differential test")
    ap.add_argument("--min-total-reads", type=int, default=20, help="Min total reads across tested fragmentforms")
    ap.add_argument("--min-tail", type=int, default=1, help="Drop reads with tail_len < this (pt:i:0 = no estimate)")
    ap.add_argument("--test", choices=["auto"], default="auto")
    ap.add_argument("--gene-filter", nargs="*", default=None, help="Restrict to these gene names")
    ap.add_argument("--figs-dir", default="", help="If set, write top-K per-gene tail-distribution PNGs here")
    ap.add_argument("--top-k", type=int, default=10, help="Number of top genes (by p_adj) to plot")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def _plot_gene(df, mg, gene, p_adj, min_reads, path, max_forms=20):
    """Boxplot of poly(A) tail length per fragmentform for one gene/metagene. The p_adj is computed
    over ALL fragmentforms with >= min_reads; when there are more than max_forms the boxplot shows only
    the most-supported ones (and says so), since dozens of boxes are unreadable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    sub = df[df["metagene_index"] == mg]
    groups = [(zt, g["tail_len"].values) for zt, g in sub.groupby("ZT", sort=False)
              if g.shape[0] >= int(min_reads)]
    if len(groups) < 2:
        return False
    n_total = len(groups)                    # number of fragmentforms in the TEST
    if n_total > max_forms:
        # keep the extremes of the distribution (shortest- and longest-median-tail fragmentforms) so
        # the plotted range is honest, PLUS the top-10 most-supported by read count; dedup by ZT.
        by_median = sorted(groups, key=lambda kv: np.median(kv[1]))
        by_reads = sorted(groups, key=lambda kv: -kv[1].size)
        picked = {}
        for grp in [by_median[-1], by_median[0]] + by_reads[:10]:   # max, min, then top-10 by reads
            picked[grp[0]] = grp
        groups = list(picked.values())
    groups.sort(key=lambda kv: np.median(kv[1]))
    labels = [zt.split(".")[-1] for zt, _ in groups]          # T1, T2, ... (transcript index)
    data = [v for _, v in groups]
    fig, ax = plt.subplots(figsize=(max(5.0, 0.7 * len(groups) + 2.0), 4.0))
    bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.62,
                    medianprops=dict(color="#c1121f", linewidth=1.4))
    for patch in bp["boxes"]:
        patch.set(facecolor="#dce7f1", edgecolor="#3b6ea5", linewidth=0.8)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color="#3b6ea5", linewidth=0.8)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("poly(A) tail length (nt)")
    cap_note = (f"\nKruskal over all {n_total} fragmentforms; showing the min- & max-tail forms + top 10 by reads"
                if n_total > len(groups)
                else f"\nKruskal over all {n_total} fragmentforms")
    ax.set_title(f"{gene} — tail length by fragmentform  (p_adj={p_adj:.1e}){cap_note}", fontsize=9)
    # add headroom above the boxes so the n= labels sit clear of the top whisker/box
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + 0.10 * (ymax - ymin))
    top = ax.get_ylim()[1]
    nfs = 5.5 if len(groups) <= 10 else 4.5   # shrink n= labels when many boxes are packed together
    for i, (_, v) in enumerate(groups):
        ax.text(i + 1, top, f"n={v.size}", ha="center", va="top", fontsize=nfs, color="#5b6773", rotation=90)
    fig.tight_layout()
    save_figure(fig, path, dpi=130, bbox_inches="tight")   # PNG + SVG
    plt.close(fig)
    return True


def _dist(vals):
    v = np.asarray(vals, dtype=float)
    return {
        "n": int(v.size),
        "median": round(float(np.median(v)), 2),
        "mean": round(float(v.mean()), 2),
        "std": round(float(v.std(ddof=1)), 2) if v.size > 1 else 0.0,
        "q25": round(float(np.percentile(v, 25)), 2),
        "q75": round(float(np.percentile(v, 75)), 2),
        "min": int(v.min()),
        "max": int(v.max()),
    }


def main():
    args = parse_args()
    usecols = ["sample", "qname", "tail_len", "ZT", "metagene_index", "gene_name", "classification"]
    try:
        df = pd.read_csv(args.in_tsv, sep="\t", low_memory=False,
                         usecols=lambda c: c in usecols)
    except Exception:
        df = pd.DataFrame(columns=usecols)

    frag_out = args.out_prefix + "_polya_fragmentform.tsv"
    diff_out = args.out_prefix + "_taillength_diffs.tsv"

    if not df.empty:
        df = df[pd.to_numeric(df["tail_len"], errors="coerce").fillna(0) >= int(args.min_tail)].copy()
        df["tail_len"] = df["tail_len"].astype(float)
        df = df[df["ZT"].astype(str).ne("")]
        if args.gene_filter:
            df = df[df["gene_name"].isin(set(args.gene_filter))]

    if df.empty:
        pd.DataFrame(columns=FRAG_COLS).to_csv(frag_out, sep="\t", index=False)
        pd.DataFrame(columns=DIFF_COLS).to_csv(diff_out, sep="\t", index=False)
        if args.verbose:
            print("[taillength] no tail rows after filtering; wrote header-only outputs", flush=True)
        return

    # ---- per-fragmentform distribution table ----
    frag_rows = []
    for zt, g in df.groupby("ZT", sort=False):
        d = _dist(g["tail_len"].values)
        frag_rows.append({
            "ZT": zt, "gene_name": g["gene_name"].iloc[0], "metagene_index": g["metagene_index"].iloc[0],
            "classification": g["classification"].iloc[0] if "classification" in g else "",
            "n_reads": d["n"], "n_samples": int(g["sample"].nunique()),
            "median_tail": d["median"], "mean_tail": d["mean"], "std_tail": d["std"],
            "q25_tail": d["q25"], "q75_tail": d["q75"], "min_tail": d["min"], "max_tail": d["max"],
        })
    frag_df = pd.DataFrame(frag_rows).sort_values(["metagene_index", "median_tail"], ascending=[True, False])
    frag_df.to_csv(frag_out, sep="\t", index=False)

    # ---- between-fragmentform differential tail length, per metagene ----
    diff_rows = []
    for mg, g in df.groupby("metagene_index", sort=False):
        by_zt = {zt: sub["tail_len"].values for zt, sub in g.groupby("ZT", sort=False)}
        kept = {zt: v for zt, v in by_zt.items() if v.size >= int(args.min_reads)}
        if len(kept) < 2:
            continue
        n_total = int(sum(v.size for v in kept.values()))
        if n_total < int(args.min_total_reads):
            continue
        groups = list(kept.values())
        if len(groups) == 2:
            stat, p = mannwhitneyu(groups[0], groups[1], alternative="two-sided")
            test_name, stat_name = "mannwhitneyu", "U"
        else:
            # kruskal raises ValueError("All numbers are identical") when every read across all
            # groups has the same tail length -- guard so a degenerate gene doesn't abort the run.
            try:
                stat, p = kruskal(*groups)
            except ValueError:
                stat, p = float("nan"), 1.0
            test_name, stat_name = "kruskal", "H"
        per_frag = []
        for zt, v in sorted(kept.items(), key=lambda kv: -np.median(kv[1])):
            d = _dist(v)
            per_frag.append({"ZT": zt, "n": d["n"], "median": d["median"], "mean": d["mean"],
                             "q25": d["q25"], "q75": d["q75"]})
        medians = [pf["median"] for pf in per_frag]
        diff_rows.append({
            "metagene_index": mg, "gene_name": g["gene_name"].iloc[0],
            "n_fragmentforms_tested": len(kept), "n_reads": n_total,
            "test_name": test_name, "stat_name": stat_name,
            "stat_value": round(float(stat), 4), "p_value": float(p),
            "effect_median_range_nt": round(float(max(medians) - min(medians)), 2),
            "min_median_tail": float(min(medians)), "max_median_tail": float(max(medians)),
            "per_fragmentform_json": json.dumps(per_frag, separators=(",", ":")),
        })

    diff_df = pd.DataFrame(diff_rows)
    if not diff_df.empty:
        diff_df["p_adj_bh"] = benjamini_hochberg(diff_df["p_value"].values)
        diff_df = diff_df.sort_values(["p_adj_bh", "effect_median_range_nt"],
                                      ascending=[True, False]).reset_index(drop=True)
        diff_df = diff_df[DIFF_COLS]
    else:
        diff_df = pd.DataFrame(columns=DIFF_COLS)
    diff_df.to_csv(diff_out, sep="\t", index=False)

    # Top-K per-gene tail-distribution figures (from the raw per-read tails still in df).
    if args.figs_dir and int(args.top_k) > 0 and not diff_df.empty:
        os.makedirs(args.figs_dir, exist_ok=True)
        n_fig = 0
        for i, r in enumerate(diff_df.head(int(args.top_k)).itertuples(index=False), start=1):
            safe = re.sub(r"[^A-Za-z0-9._-]", "_", str(r.gene_name))
            path = os.path.join(args.figs_dir, f"rank{i:02d}_{safe}.png")
            if _plot_gene(df, r.metagene_index, r.gene_name, r.p_adj_bh, args.min_reads, path):
                n_fig += 1
        if args.verbose:
            print(f"[taillength] wrote {n_fig} per-gene tail figures -> {args.figs_dir}", flush=True)

    if args.verbose:
        print(f"[taillength] fragmentforms={len(frag_df)} metagenes_tested={len(diff_df)} -> {diff_out}", flush=True)


if __name__ == "__main__":
    main()
