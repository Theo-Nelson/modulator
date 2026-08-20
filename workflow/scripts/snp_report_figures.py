"""Per-example figures for the genotype/SNP sections of the modulator HTML report.

For each top row of the SNP->transcript, SNP->mod, SNP x transcript x mod
dependency, and haplotype association tables, render one small figure showing the
actual relationship for that individual SNP / haplotype block (e.g. how a SNP's
alleles partition across transcripts, or how methylation differs by allele within
each transcript). Figures are written to disk under <figs_dir>/<section>/ and also
returned as inline base64 data-URIs for embedding in the report.

All builders are defensive: any parse/plot failure for a row is skipped, and if
matplotlib is unavailable the gallery is simply omitted.
"""
from __future__ import annotations

import base64
import html
import json
import os
import re
from io import BytesIO

import pandas as pd

_EFFECT_COLS = [
    "effect_max_abs_tx_frac_diff", "effect_abs_delta_mod_frac",
    "weighted_within_tx_effect", "overall_effect_abs_delta_mod_frac",
    "effect_max_abs_mod_rate_diff",
]
_PADJ_COLS = ["p_adj_bh"]
_READS_COLS = ["n_reads", "support_reads", "n_alt_reads", "complete_reads"]


def _safe(s, n=70):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))[:n].strip("_") or "fig"


def _alleles(snp_id):
    # snp_id like "chr11:124746393:T>G"
    try:
        last = str(snp_id).split(":")[-1]
        ref, alt = last.split(">")
        return ref, alt
    except Exception:
        return "ref", "alt"


def _short_zt(zt):
    zt = str(zt)
    return zt.split(".")[-1] if "." in zt else zt  # e.g. NRGN.NRGN.G347.T1 -> T1


def _load_json(val):
    try:
        return json.loads(val)
    except Exception:
        return None


def _frac(mod, unmod):
    tot = float(mod or 0) + float(unmod or 0)
    return (float(mod or 0) / tot if tot else 0.0), tot


def _rank_pick(df, max_figs, min_reads):
    work = df.copy()
    reads_col = next((c for c in _READS_COLS if c in work.columns), None)
    if reads_col is not None:
        floored = work[pd.to_numeric(work[reads_col], errors="coerce").fillna(0) >= min_reads]
        if not floored.empty:
            work = floored
    padj = next((c for c in _PADJ_COLS if c in work.columns), None)
    eff = next((c for c in _EFFECT_COLS if c in work.columns), None)
    cols, asc = [], []
    if padj:
        cols.append(padj); asc.append(True)
    if eff:
        cols.append(eff); asc.append(False)
    if cols:
        work = work.sort_values(cols, ascending=asc, kind="mergesort")
    return work.head(max_figs)


# --------------------------- per-row figure builders ---------------------------

def _fig_snp_transcript(row, plt):
    import numpy as np
    data = _load_json(row.get("per_transcript_json"))
    if not data:
        return None
    txs = [_short_zt(d.get("ZT", "?")) for d in data]
    ref = [float(d.get("ref_reads", 0)) for d in data]
    alt = [float(d.get("alt_reads", 0)) for d in data]
    rs, as_ = (sum(ref) or 1.0), (sum(alt) or 1.0)
    ref_f = [r / rs for r in ref]
    alt_f = [a / as_ for a in alt]
    refb, altb = _alleles(row.get("snp_id"))
    x = np.arange(len(txs)); w = 0.38
    fig, ax = plt.subplots(figsize=(max(4.4, 1.15 * len(txs) + 2.0), 3.4))
    ax.bar(x - w / 2, ref_f, w, label=f"ref ({refb})", color="#4878a8")
    ax.bar(x + w / 2, alt_f, w, label=f"alt ({altb})", color="#c0552f")
    for xi, (rf, af, rc, ac) in enumerate(zip(ref_f, alt_f, ref, alt)):
        ax.text(xi - w / 2, rf + 0.01, f"{int(rc)}", ha="center", va="bottom", fontsize=7)
        ax.text(xi + w / 2, af + 0.01, f"{int(ac)}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(txs, fontsize=8, rotation=20, ha="right")
    ax.set_ylabel("fraction of allele's reads"); ax.set_ylim(0, 1.14)
    gene = str(row.get("gene_names", "")).split(",")[0]
    ax.set_title(f"{row.get('snp_id','')}  {gene}\nallele → transcript usage", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    return fig, f"{row.get('snp_id','')} {gene} (allele→transcript)"


def _fig_snp_mod(row, plt):
    data = _load_json(row.get("per_state_json"))
    if not data:
        return None
    rf, rn = _frac(data.get("ref_modified"), data.get("ref_not_target"))
    af, an = _frac(data.get("alt_modified"), data.get("alt_not_target"))
    refb, altb = _alleles(row.get("snp_id"))
    fig, ax = plt.subplots(figsize=(4.0, 3.4))
    bars = ax.bar([f"ref\n{refb}", f"alt\n{altb}"], [rf, af],
                  color=["#4878a8", "#c0552f"], width=0.6)
    for b, fr, tot in zip(bars, [rf, af], [rn, an]):
        ax.text(b.get_x() + b.get_width() / 2, fr + 0.01,
                f"{fr*100:.0f}%\n(n={int(tot)})", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, 1.16); ax.set_ylabel(f"fraction modified ({row.get('target_mod_code','')})")
    gene = str(row.get("gene_names", "")).split(",")[0]
    ax.set_title(f"{row.get('snp_id','')}  {gene}\n{row.get('mod_site_id','')}", fontsize=9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    return fig, f"{row.get('snp_id','')} {gene} {row.get('target_mod_code','')} (allele→mod)"


def _fig_hap(row, plt):
    import numpy as np
    data = _load_json(row.get("per_table_json"))
    if not data or "haplotypes" not in data or "counts" not in data:
        return None
    haps = [str(h) for h in data["haplotypes"]]
    states = [str(s) for s in data.get("states", [])]
    counts = np.asarray(data["counts"], dtype=float)
    if counts.ndim != 2 or counts.shape[0] != len(haps):
        return None
    block = row.get("block_id", "")
    is_mod = states[:2] == ["modified", "not_target"] or "mod_site_id" in row.index
    fig, ax = plt.subplots(figsize=(max(4.2, 0.95 * len(haps) + 2.0), 3.5))
    x = np.arange(len(haps))
    if is_mod and counts.shape[1] >= 2:
        fr = [counts[i, 0] / (counts[i].sum() or 1) for i in range(len(haps))]
        tot = [int(counts[i].sum()) for i in range(len(haps))]
        bars = ax.bar(x, fr, 0.62, color="#6a4c93")
        for b, f, t in zip(bars, fr, tot):
            ax.text(b.get_x() + b.get_width() / 2, f + 0.01, f"{f*100:.0f}%\n(n={t})",
                    ha="center", va="bottom", fontsize=7.5)
        ax.set_ylim(0, 1.16); ax.set_ylabel(f"fraction modified ({row.get('target_mod_code','')})")
        title = f"block {block}  {row.get('mod_site_id','')}\nmethylation by haplotype"
    else:
        # fragmentform usage: stacked proportion of each haplotype's reads across the gene's fragmentforms
        rowsum = counts.sum(axis=1, keepdims=True); rowsum[rowsum == 0] = 1
        prop = counts / rowsum
        bottom = np.zeros(len(haps))
        cmap = plt.get_cmap("tab20")
        for j in range(prop.shape[1]):
            lbl = _short_zt(states[j]) if j < len(states) else f"s{j}"
            if lbl.startswith("s") and lbl[1:].isdigit():   # generic state -> fragmentform label T1, T2, ...
                lbl = f"T{int(lbl[1:]) + 1}"
            ax.bar(x, prop[:, j], 0.62, bottom=bottom, label=lbl, color=cmap(j % 20))
            bottom += prop[:, j]
        ax.set_ylim(0, 1.0); ax.set_ylabel("fragmentform usage fraction")
        ax.legend(fontsize=6.5, frameon=False, ncol=1, bbox_to_anchor=(1.01, 1), loc="upper left",
                  title="fragmentform")
        gene = str(row.get("gene_names", "") or "").split(",")[0]
        title = (f"block {block} — {gene}\nfragmentform usage by haplotype" if gene
                 else f"block {block}\nfragmentform usage by haplotype")
    ax.set_xticks(x); ax.set_xticklabels(haps, fontsize=8, rotation=20, ha="right")
    ax.set_title(title, fontsize=8.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    return fig, f"block_{block}"


_KINDS = {
    "snp_tx":     ("_fig_snp_transcript", "SNP→transcript example(s)"),
    "snp_mod":    ("_fig_snp_mod",        "SNP→modification example(s)"),
    "hap_tx":     ("_fig_hap",            "haplotype→transcript example(s)"),
    "hap_mod":    ("_fig_hap",            "haplotype→modification example(s)"),
}


def _gallery(df, kind, figs_dir, max_figs, min_reads):
    if df is None or df.empty or max_figs <= 0:
        return ""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""
    builder = globals()[_KINDS[kind][0]]
    picked = _rank_pick(df, max_figs, min_reads)
    if picked.empty:
        return ""
    out_dir = os.path.join(figs_dir, kind) if figs_dir else None
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    pieces = []
    for rank, (_, row) in enumerate(picked.iterrows(), 1):
        try:
            built = builder(row, plt)
        except Exception:
            built = None
        if not built:
            continue
        fig, caption = built
        try:
            from plot_utils import setup_matplotlib_style, bump_fonts
            setup_matplotlib_style(); bump_fonts(fig)   # Arial-like + enlarged fonts
        except Exception:
            pass
        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
        if out_dir:  # vector PDF + SVG companions, written from the live figure before it is closed
            for _fmt in ("pdf", "svg"):
                try:
                    fig.savefig(os.path.join(out_dir, f"rank{rank:02d}__{_safe(caption)}.{_fmt}"),
                                format=_fmt, bbox_inches="tight")
                except Exception:
                    pass
        plt.close(fig)
        png = buf.getvalue()
        if out_dir:
            try:
                with open(os.path.join(out_dir, f"rank{rank:02d}__{_safe(caption)}.png"), "wb") as fh:
                    fh.write(png)
            except Exception:
                pass
        uri = "data:image/png;base64," + base64.b64encode(png).decode("ascii")
        cap = html.escape(caption)
        pieces.append(
            f"<figure class='image-card'>"
            f"<a class='image-link' href='{uri}' target='_blank' rel='noopener noreferrer' title='Open full-size image in a new tab'>"
            f"<img src='{uri}' alt='{cap}' /><span class='expand-badge' aria-hidden='true'>↗</span></a>"
            f"<figcaption>{cap}</figcaption></figure>"
        )
    if not pieces:
        return ""
    label = _KINDS[kind][1]
    return (
        f"<details class='definitions' open><summary>Top {len(pieces)} {label}</summary>"
        f"<div class='gallery'>" + "".join(pieces) + "</div></details>"
    )


def build_snp_galleries(*, snp_tx=None, snp_mod=None,
                        hap_tx=None, hap_mod=None, figs_dir="", max_figs=12, min_reads=10):
    """Return {section_key: gallery_html} for the per-example SNP figures.

    section_key in {snp_tx, snp_mod, hap_tx, hap_mod}. Missing/empty
    tables yield "" for that key.
    """
    return {
        "snp_tx":     _gallery(snp_tx,     "snp_tx",     figs_dir, max_figs, min_reads),
        "snp_mod":    _gallery(snp_mod,    "snp_mod",    figs_dir, max_figs, min_reads),
        "hap_tx":     _gallery(hap_tx,     "hap_tx",     figs_dir, max_figs, min_reads),
        "hap_mod":    _gallery(hap_mod,    "hap_mod",    figs_dir, max_figs, min_reads),
    }
