#!/usr/bin/env python3
"""Poly(A) tail length x modification.

Per modification site, do reads that are MODIFIED carry a different poly(A) tail length than reads
that are UNMODIFIED at that site? Joins the per-read tail table (build_read_polya_table.py) to the
per-read mod-call table (molecule_mod_calls.tsv) on (sample, qname), then per mod_site_id compares
tail lengths of modified vs unmodified reads with a Mann-Whitney U test (continuous). This is the
tail x m6A/mod mechanistic readout (tail length is linked to modification and stability).

The mod table is genotype-scale (tens of GB), so it is sharded per chromosome (shard_tsv_by_chrom,
O(1) RAM router) and streamed with usecols; the per-read tail table is small and held as a
(sample,qname) -> tail_len dict. Empty-input safe (header-only output).
"""
import argparse
import json
import os
import re
import shutil
import tempfile

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from genotype_utils import benjamini_hochberg, shard_tsv_by_chrom, tsv_header
from plot_utils import save_figure

_MOD_WANT = ["sample", "qname", "mod_site_id", "chrom", "target_mod_code", "target_modified",
             "gene_name", "state_detail", "usable", "fail", "within_alignment", "ZN"]
OUT_COLS = ["mod_site_id", "chrom", "gene_name", "target_mod_code",
            "n_reads", "n_modified", "n_unmodified", "median_tail_modified", "median_tail_unmodified",
            "effect_median_diff_nt", "test_name", "stat_name", "stat_value", "p_value", "p_adj_bh",
            "n_forms_comparable", "n_forms_concordant", "per_fragmentform_json"]


def parse_args():
    ap = argparse.ArgumentParser(description="Poly(A) tail length vs modification state, per mod site.")
    ap.add_argument("--tail-tsv", required=True, help="Per-read tail table from build_read_polya_table.py")
    ap.add_argument("--mod-tsv", required=True, help="Per-read molecule_mod_calls.tsv")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-state-reads", type=int, default=10, help="Min reads per state (modified / unmodified)")
    ap.add_argument("--min-tail", type=int, default=1, help="Drop reads with tail_len < this (pt:i:0 = no estimate)")
    ap.add_argument("--figs-dir", default="", help="If set, write top-K per-site tail-distribution PNGs here")
    ap.add_argument("--top-k", type=int, default=10, help="Number of top sites (by p_adj) to plot")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def _plot_site(mod_t, unmod_t, meta, path, per_ff=None):
    """Modified-vs-unmodified poly(A) tail for one mod site: (left) pooled histogram across all
    fragmentforms, and (right, when available) a per-fragmentform dumbbell showing the modified vs
    unmodified MEDIAN tail WITHIN each fragmentform -- so a pooled shift can be told apart from the
    modification merely tracking a differently-tailed fragmentform."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    # fragmentforms with BOTH states well-sampled (>=3 reads each)
    ff = [(zn, d) for zn, d in (per_ff or {}).items()
          if d.get("median_tail_mod") is not None and d.get("median_tail_unmod") is not None
          and d["n_mod"] >= 3 and d["n_unmod"] >= 3]
    n_ff_total = len(per_ff or {})               # all fragmentforms with ANY reads at the site
    ff.sort(key=lambda t: t[1]["delta_nt"])       # by within-form Δ, so a lone strong-effect form stands out
    two = len(ff) >= 1
    if two:
        fig, (ax, a2) = plt.subplots(1, 2, figsize=(11.2, max(4.2, 0.32 * len(ff) + 2.2)), layout="constrained")
    else:
        fig, ax = plt.subplots(figsize=(6.4, 4.0), layout="constrained")
    hi = float(np.percentile(np.concatenate([mod_t, unmod_t]), 99))
    bins = np.linspace(0, max(hi, 10), 40)
    ax.hist(unmod_t, bins=bins, density=True, color="#8aa0b5", alpha=0.6, label=f"unmodified (n={unmod_t.size})")
    ax.hist(mod_t, bins=bins, density=True, color="#c1121f", alpha=0.55, label=f"modified (n={mod_t.size})")
    ax.axvline(np.median(unmod_t), color="#41576d", ls="--", lw=1.2)
    ax.axvline(np.median(mod_t), color="#7a0c15", ls="--", lw=1.2)
    ax.set_xlabel("poly(A) tail length (nt)"); ax.set_ylabel("density")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("pooled across fragmentforms", fontsize=9)
    if two:
        xmax = max(max(d["median_tail_unmod"], d["median_tail_mod"]) for _, d in ff)
        for y, (zn, d) in enumerate(ff):
            um, mm = d["median_tail_unmod"], d["median_tail_mod"]
            a2.plot([um, mm], [y, y], color="#b8c2cc", lw=1.6, zorder=1)
            a2.scatter([um], [y], color="#8aa0b5", s=34, zorder=2)
            a2.scatter([mm], [y], color="#c1121f", s=34, zorder=2)
            a2.text(xmax * 1.04, y, f"Δ={d['delta_nt']:+.0f}", va="center", ha="left", fontsize=6.5,
                    color="#7a0c15" if d["delta_nt"] < 0 else "#1d6b4f")
        a2.set_yticks(range(len(ff)))
        a2.set_yticklabels([f"ZN{zn} (n {d['n_mod']}/{d['n_unmod']})" for zn, d in ff], fontsize=7.5)
        a2.set_xlim(right=xmax * 1.22)
        a2.set_xlabel("median poly(A) tail (nt)")
        # be explicit that only forms with BOTH states well-sampled are comparable here
        a2.set_title(f"per fragmentform  (● unmod, ● mod, Δ = mod−unmod)\n"
                     f"{len(ff)} of {n_ff_total} forms with ≥3 reads in each state", fontsize=8)
        a2.grid(True, axis="x", ls="--", lw=0.5, alpha=0.4)
        for sp in ("top", "right"):
            a2.spines[sp].set_visible(False)
    d = float(np.median(mod_t) - np.median(unmod_t))
    # warn when the pooled effect is NOT reproduced within fragmentforms (single-form-driven / confounded)
    n_comp = len(ff)
    # a zero-Δ form is neither concordant nor discordant -- excluding it keeps the count independent
    # of the pooled sign (else Δ=0 forms counted concordant only when pooled<0: False==False).
    n_conc = sum(1 for _, dd in ff if dd["delta_nt"] != 0 and (dd["delta_nt"] > 0) == (d > 0)) if d != 0 else 0
    warn = ""
    if n_comp <= 1:
        warn = (f"CAUTION: only {n_comp} fragmentform comparable — pooled effect cannot be separated "
                f"from fragmentform identity")
    elif n_conc < (n_comp + 1) // 2:
        warn = (f"CAUTION: not consistent across fragmentforms — only {n_conc}/{n_comp} shift in the "
                f"pooled direction (effect may be driven by a subset)")
    sup = (f"{meta['gene_name']}  {meta['mod_site_id']} — modified median {np.median(mod_t):.0f} vs "
           f"unmodified {np.median(unmod_t):.0f} nt  (Δ={d:+.0f}, p_adj={meta['p_adj_bh']:.1e})")
    # a dark-red title (with the caveat as a second line) flags "read carefully"; kept in the
    # layout-managed suptitle so it never overlaps the axes, and ASCII-only for font safety.
    fig.suptitle(sup + ("\n" + warn if warn else ""), fontsize=9, color=("#b00020" if warn else "black"))
    save_figure(fig, path, dpi=130, bbox_inches="tight")   # PNG + SVG
    plt.close(fig)
    return True


def _load_tail_map(path, min_tail):
    """(sample,qname) -> tail_len for reads with an estimated tail."""
    cols = ["sample", "qname", "tail_len"]
    df = pd.read_csv(path, sep="\t", low_memory=False, usecols=lambda c: c in cols)
    if df.empty:
        return {}
    df = df[pd.to_numeric(df["tail_len"], errors="coerce").fillna(0) >= int(min_tail)]
    key = df["sample"].astype(str) + "\x00" + df["qname"].astype(str)
    return dict(zip(key, df["tail_len"].astype(float)))


def _site_rows_for_chrom(mod_path, tail_map, args):
    hdr = tsv_header(mod_path)
    mod = pd.read_csv(mod_path, sep="\t", low_memory=False, usecols=[c for c in _MOD_WANT if c in hdr])
    if mod.empty:
        return [], []
    if "usable" in mod.columns:
        mod = mod[mod["usable"].fillna(False)].copy()
    elif "fail" in mod.columns and "within_alignment" in mod.columns:
        mod = mod[(~mod["fail"].fillna(True)) & mod["within_alignment"].fillna(False)].copy()
    if mod.empty:
        return [], []
    key = mod["sample"].astype(str) + "\x00" + mod["qname"].astype(str)
    mod["tail_len"] = key.map(tail_map).to_numpy()
    mod = mod[mod["tail_len"].notna()]
    if mod.empty:
        return [], []
    mod["target_modified"] = pd.to_numeric(mod["target_modified"], errors="coerce").fillna(0).astype(int)
    collect = bool(args.figs_dir) and int(args.top_k) > 0
    rows, cands = [], []
    for sid, g in mod.groupby("mod_site_id", sort=False):
        mod_t = g.loc[g["target_modified"] == 1, "tail_len"].values
        unmod_t = g.loc[g["target_modified"] == 0, "tail_len"].values
        if mod_t.size < int(args.min_state_reads) or unmod_t.size < int(args.min_state_reads):
            continue
        stat, p = mannwhitneyu(mod_t, unmod_t, alternative="two-sided")
        f = g.iloc[0]
        # per-FRAGMENTFORM (ZN) breakdown: does the modified-vs-unmodified tail gap hold WITHIN a
        # fragmentform, or is it just the modification tracking a differently-tailed fragmentform?
        per_ff = {}
        if "ZN" in g.columns:
            for zn, zg in g.groupby("ZN"):
                try:
                    zi = int(float(zn))
                except (TypeError, ValueError):
                    continue  # unassigned reads (blank ZN) belong to no fragmentform
                mt = zg.loc[zg["target_modified"] == 1, "tail_len"].values
                ut = zg.loc[zg["target_modified"] == 0, "tail_len"].values
                if mt.size == 0 and ut.size == 0:
                    continue
                per_ff[zi] = {
                    "n_mod": int(mt.size), "n_unmod": int(ut.size),
                    "median_tail_mod": round(float(np.median(mt)), 1) if mt.size else None,
                    "median_tail_unmod": round(float(np.median(ut)), 1) if ut.size else None,
                    "delta_nt": (round(float(np.median(mt) - np.median(ut)), 1)
                                 if mt.size and ut.size else None),
                }
        row = {
            "mod_site_id": sid, "chrom": f.get("chrom", ""), "gene_name": f.get("gene_name", ""),
            "target_mod_code": f.get("target_mod_code", ""),
            "n_reads": int(mod_t.size + unmod_t.size), "n_modified": int(mod_t.size), "n_unmodified": int(unmod_t.size),
            "median_tail_modified": round(float(np.median(mod_t)), 2),
            "median_tail_unmodified": round(float(np.median(unmod_t)), 2),
            "effect_median_diff_nt": round(float(np.median(mod_t) - np.median(unmod_t)), 2),
            "test_name": "mannwhitneyu", "stat_name": "U", "stat_value": round(float(stat), 4), "p_value": float(p),
            "per_fragmentform_json": json.dumps(per_ff),
        }
        # Is the pooled shift reproduced WITHIN fragmentforms, or driven by a single form (which would
        # make it inseparable from fragmentform identity)? Count comparable (>=3 reads/state) forms and
        # how many shift in the SAME direction as the pooled effect.
        _pooled = float(np.median(mod_t) - np.median(unmod_t))
        _comp = [d for d in per_ff.values()
                 if d["n_mod"] >= 3 and d["n_unmod"] >= 3 and d["delta_nt"] is not None]
        row["n_forms_comparable"] = len(_comp)
        row["n_forms_concordant"] = (sum(1 for d in _comp
                                         if d["delta_nt"] != 0 and (d["delta_nt"] > 0) == (_pooled > 0))
                                     if _pooled != 0 else 0)
        rows.append(row)
        if collect:
            cands.append((float(p), row, mod_t, unmod_t, per_ff))
    return rows, cands


def main():
    args = parse_args()
    if not (os.path.exists(args.mod_tsv) and os.path.getsize(args.mod_tsv)
            and os.path.exists(args.tail_tsv) and os.path.getsize(args.tail_tsv)):
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return
    tail_map = _load_tail_map(args.tail_tsv, args.min_tail)
    if args.verbose:
        print(f"[taillength_mod] {len(tail_map)} reads with a tail estimate", flush=True)
    if not tail_map:
        pd.DataFrame(columns=OUT_COLS).to_csv(args.out_tsv, sep="\t", index=False)
        return

    keep = max(1, 4 * int(args.top_k))  # bound the figure-candidate stash to ~4*top_k sites' raw arrays
    tmp = tempfile.mkdtemp(prefix=".tailmod_", dir=os.path.dirname(args.out_tsv) or ".")
    rows, fig_cands = [], []
    try:
        shards = shard_tsv_by_chrom(args.mod_tsv, os.path.join(tmp, "mod"))
        for chrom in sorted(shards):
            r, c = _site_rows_for_chrom(shards[chrom], tail_map, args)
            rows.extend(r)
            if c:
                fig_cands.extend(c)
                fig_cands.sort(key=lambda x: x[0])   # smallest p first (== smallest p_adj; BH is monotonic)
                del fig_cands[keep:]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_adj_bh"] = benjamini_hochberg(out["p_value"].values)
        out["_abs"] = out["effect_median_diff_nt"].abs()
        out = out.sort_values(["p_adj_bh", "_abs"], ascending=[True, False]).drop(columns="_abs").reset_index(drop=True)
        out = out[OUT_COLS]
    else:
        out = pd.DataFrame(columns=OUT_COLS)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    # Top-K per-site modified-vs-unmodified tail figures.
    if args.figs_dir and int(args.top_k) > 0 and fig_cands:
        os.makedirs(args.figs_dir, exist_ok=True)
        padj = dict(zip(out["mod_site_id"], out["p_adj_bh"])) if not out.empty else {}
        n_fig = 0
        for i, (p, row, mod_t, unmod_t, per_ff) in enumerate(sorted(fig_cands, key=lambda x: x[0])[:int(args.top_k)], start=1):
            meta = dict(row); meta["p_adj_bh"] = float(padj.get(row["mod_site_id"], p))
            safe = re.sub(r"[^A-Za-z0-9._-]", "_", f"{row['gene_name']}_{row['mod_site_id']}")
            if _plot_site(mod_t, unmod_t, meta, os.path.join(args.figs_dir, f"rank{i:02d}_{safe}.png"), per_ff):
                n_fig += 1
        if args.verbose:
            print(f"[taillength_mod] wrote {n_fig} per-site tail figures -> {args.figs_dir}", flush=True)
    if args.verbose:
        print(f"[taillength_mod] wrote {len(out)} site rows -> {args.out_tsv}", flush=True)


if __name__ == "__main__":
    main()
