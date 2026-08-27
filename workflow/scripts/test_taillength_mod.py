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
from scipy.stats import mannwhitneyu, norm as _norm, rankdata as _rankdata


def van_elteren_stratified(strata):
    """Van Elteren stratified Wilcoxon rank-sum test over a list of (mod_tails, unmod_tails) strata,
    one per (sample, fragmentform). Design-free weight w_k = 1/(N_k + 1).

    BLOCKER-2: the pooled Mann-Whitney compares modified vs unmodified reads across ALL fragmentforms,
    samples and conditions at once. When isoform usage and modification are both condition-linked and
    the isoforms have different tails, modification status proxies for isoform and the isoform's tail
    difference is attributed to the modification -- a +33 nt "effect" at FDR 1.3e-7 on the fixture where
    every within-stratum difference is ~0. Stratifying by (sample, fragmentform) removes the confound;
    a stratum with only one modification state carries no within-stratum information and is dropped.

    Returns (z, p, n_informative_strata); NaN when no stratum has both states."""
    num = 0.0
    var = 0.0
    used = 0
    for mt, ut in strata:
        n1, n2 = len(mt), len(ut)
        N = n1 + n2
        if n1 < 2 or n2 < 2:
            continue  # a 1-vs-1 (or 1-vs-n) stratum contributes a near-deterministic rank with tiny
                      # variance -- statistical noise, not signal; require >=2 reads in EACH state
        allv = np.concatenate([np.asarray(mt, float), np.asarray(ut, float)])
        ranks = _rankdata(allv)
        W = float(ranks[:n1].sum())            # rank sum of the MODIFIED group
        EW = n1 * (N + 1) / 2.0
        _, cnt = np.unique(allv, return_counts=True)
        tie = float((cnt ** 3 - cnt).sum())
        VW = n1 * n2 / 12.0 * ((N + 1) - tie / (N * (N - 1)))
        w = 1.0 / (N + 1)
        num += w * (W - EW)
        var += w * w * VW
        used += 1
    if used == 0 or var <= 0:
        return float("nan"), float("nan"), used
    z = num / np.sqrt(var)
    return float(z), float(2.0 * _norm.sf(abs(z))), used


def weighted_within_stratum_median_diff(strata):
    """Effect size consistent with the stratified test: within each (sample, fragmentform) stratum,
    median(mod) - median(unmod), averaged over strata weighted by stratum size."""
    num = den = 0.0
    for mt, ut in strata:
        if len(mt) < 1 or len(ut) < 1:
            continue
        w = len(mt) + len(ut)
        num += w * (float(np.median(mt)) - float(np.median(ut)))
        den += w
    return num / den if den > 0 else 0.0

from genotype_utils import benjamini_hochberg, shard_tsv_by_chrom, tsv_header
from plot_utils import save_figure

_MOD_WANT = ["sample", "qname", "mod_site_id", "chrom", "target_mod_code", "target_modified",
             "gene_name", "state_detail", "usable", "fail", "within_alignment", "ZN", "ZT"]
OUT_COLS = ["mod_site_id", "chrom", "gene_name", "target_mod_code",
            "n_reads", "n_modified", "n_unmodified", "median_tail_modified", "median_tail_unmodified",
            "effect_median_diff_nt", "test_name", "stat_name", "stat_value", "p_value", "p_adj_bh",
            "n_strata_informative", "n_forms_comparable", "n_forms_concordant", "tailmod_confounded",
            "test_name_pooled", "p_value_pooled", "effect_median_diff_nt_pooled", "per_fragmentform_json"]


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
    """(sample,qname) -> tail_len for reads with an ESTIMATED tail (pt:i:0 = no estimate, excluded)."""
    cols = ["sample", "qname", "tail_len", "tail_estimated"]
    df = pd.read_csv(path, sep="\t", low_memory=False, usecols=lambda c: c in cols)
    if df.empty:
        return {}
    # pt:i:0 is dorado's "no estimate" sentinel, not a 0-nt tail: require an actual estimate.
    _tl = pd.to_numeric(df["tail_len"], errors="coerce").fillna(0)
    if "tail_estimated" in df.columns:
        _est = df["tail_estimated"].astype(str).str.strip().str.lower().isin(("true", "1", "yes"))
    else:
        _est = _tl > 0
    df = df[_est & (_tl >= int(min_tail))]
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
    # a read at a site claimed by >1 overlapping gene appears once PER gene (mod_site_id omits gene_id),
    # so dedup by (sample, qname, mod_site_id) -- this test, unlike its siblings, does not otherwise
    # dedup by read, so duplicates were double-counted and flipped borderline significance.
    if {"sample", "qname", "mod_site_id"}.issubset(mod.columns):
        mod = mod.drop_duplicates(["sample", "qname", "mod_site_id"])
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
        # pooled Mann-Whitney (kept as *_pooled -- it is the confounded statistic, see BLOCKER-2)
        try:
            stat, p_pooled = mannwhitneyu(mod_t, unmod_t, alternative="two-sided")
        except ValueError:
            stat, p_pooled = float("nan"), float("nan")
        # PRIMARY: van Elteren stratified Wilcoxon over (sample, fragmentform). Strata remove the
        # isoform/condition confound; a site with no stratum carrying BOTH states is untestable (its
        # entire mod-vs-unmod gap is between fragmentforms/samples) and is dropped.
        _scol = "sample" if "sample" in g.columns else None
        _ffc = "ZT" if "ZT" in g.columns else ("ZN" if "ZN" in g.columns else None)
        strata = []
        if _scol and _ffc:
            for _sk, sg in g.groupby([_scol, _ffc]):
                if str(_sk[1]).strip() in ("", "nan"):
                    continue  # unassigned fragmentform
                smt = sg.loc[sg["target_modified"] == 1, "tail_len"].values
                sut = sg.loc[sg["target_modified"] == 0, "tail_len"].values
                if smt.size >= 2 and sut.size >= 2:   # match van_elteren_stratified's >=2/state guard
                    strata.append((smt, sut))
        z_strat, p, n_strata = van_elteren_stratified(strata)
        if not np.isfinite(p):
            continue  # no (sample, fragmentform) stratum has both states -> confounded, untestable
        f = g.iloc[0]
        # per-FRAGMENTFORM (ZN) breakdown: does the modified-vs-unmodified tail gap hold WITHIN a
        # fragmentform, or is it just the modification tracking a differently-tailed fragmentform?
        per_ff = {}
        # group by the FRAGMENTFORM (ZT), NOT ZN: ZN is a graph colour that non-overlapping
        # fragmentforms of a gene SHARE, so grouping by it pools them and computes the per-fragmentform
        # confounder guard on gene-pooled reads. Fall back to ZN only if ZT is absent.
        _ffcol = "ZT" if "ZT" in g.columns else ("ZN" if "ZN" in g.columns else None)
        if _ffcol:
            for ffval, zg in g.groupby(_ffcol):
                if _ffcol == "ZN":
                    try:
                        ffkey = int(float(ffval))
                    except (TypeError, ValueError):
                        continue  # unassigned reads (blank ZN) belong to no fragmentform
                else:
                    ffkey = str(ffval).strip()
                    if not ffkey:
                        continue  # unassigned reads (blank ZT) belong to no fragmentform
                mt = zg.loc[zg["target_modified"] == 1, "tail_len"].values
                ut = zg.loc[zg["target_modified"] == 0, "tail_len"].values
                if mt.size == 0 and ut.size == 0:
                    continue
                per_ff[ffkey] = {
                    "n_mod": int(mt.size), "n_unmod": int(ut.size),
                    "median_tail_mod": round(float(np.median(mt)), 1) if mt.size else None,
                    "median_tail_unmod": round(float(np.median(ut)), 1) if ut.size else None,
                    "delta_nt": (round(float(np.median(mt) - np.median(ut)), 1)
                                 if mt.size and ut.size else None),
                }
        row = {
            # a mod site can sit in a metagene spanning >1 gene; label with ALL of them (sorted, unique)
            # rather than the arbitrary first row -- the bug already fixed in test_taillength_diffs.
            "mod_site_id": sid, "chrom": f.get("chrom", ""),
            "gene_name": "+".join(sorted(g["gene_name"].dropna().astype(str).unique())) or "",
            "target_mod_code": f.get("target_mod_code", ""),
            "n_reads": int(mod_t.size + unmod_t.size), "n_modified": int(mod_t.size), "n_unmodified": int(unmod_t.size),
            "median_tail_modified": round(float(np.median(mod_t)), 2),
            "median_tail_unmodified": round(float(np.median(unmod_t)), 2),
            # primary effect = within-(sample, fragmentform) median difference, weighted by stratum size
            "effect_median_diff_nt": round(weighted_within_stratum_median_diff(strata), 2),
            "test_name": "van_elteren_stratified_wilcoxon", "stat_name": "z",
            "stat_value": round(float(z_strat), 4), "p_value": float(p),
            "n_strata_informative": int(n_strata),
            # pooled companions (the OLD confounded statistic -- do not rank on these)
            "test_name_pooled": "mannwhitneyu",
            "p_value_pooled": float(p_pooled),
            "effect_median_diff_nt_pooled": round(float(np.median(mod_t) - np.median(unmod_t)), 2),
            "per_fragmentform_json": json.dumps(per_ff),
        }
        # Secondary confounder flag (the stratified p is the primary guard now). Is the effect reproduced
        # WITHIN fragmentforms, or driven by a single form? Count comparable (>=3 reads/state) forms and
        # how many shift in the SAME direction as the STRATIFIED effect. A form is discordant if it moves
        # the other way OR its magnitude is <1/3 of the effect (near-flat). tailmod_confounded flags any
        # site where a strict majority of comparable forms are NOT concordant -- excluded from the report
        # headline count.
        _eff = row["effect_median_diff_nt"]
        _comp = [d for d in per_ff.values()
                 if d["n_mod"] >= 3 and d["n_unmod"] >= 3 and d["delta_nt"] is not None]
        row["n_forms_comparable"] = len(_comp)
        row["n_forms_concordant"] = (sum(1 for d in _comp
                                         if d["delta_nt"] != 0 and (d["delta_nt"] > 0) == (_eff > 0)
                                         and abs(d["delta_nt"]) >= abs(_eff) / 3.0)
                                     if _eff != 0 else 0)
        # confounded unless a STRICT MAJORITY of comparable forms are concordant. `2*n_conc <= len`
        # (concordant <= half) is the correct predicate: `n_conc < len/2` swung too far -- for a 2-form
        # site it needs n_conc < 1, i.e. 0 concordant, so the 1-of-2-discordant case the original finding
        # called out could never flag. `2*n_conc <= len` flags 1-of-2 (2*1<=2) while still passing a
        # genuine 2-of-3 or 3-of-4 concordant majority.
        row["tailmod_confounded"] = bool(_comp and 2 * row["n_forms_concordant"] <= len(_comp))
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
