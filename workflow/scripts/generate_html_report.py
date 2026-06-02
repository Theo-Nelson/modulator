#!/usr/bin/env python3

import argparse
import base64
import glob
import html
import os
from pathlib import Path

import pandas as pd


CARD_DEFINITIONS = {
    "Assembled transcripts": "Final transcript models retained after transcript assembly, support, and poly(A)-related filters.",
    "Gene buckets": "Distinct gene-level transcript groups represented in the final classification summary.",
    "Metagenes": "Overlapping gene groups used for shared ZN partitioning.",
    "Max ZN partitions / metagene": "Largest number of non-overlapping ZN partitions required by any one metagene.",
    "Assigned reads": "Total reads assigned to retained transcript models in the classification summary.",
    "Exact-chain reads": "Assigned reads that exactly match the retained canonical intron chain.",
    "Truncation-assigned reads": "Assigned reads absorbed into a retained transcript from suffix-compatible shorter chains.",
    "Segregating SNPs": "Candidate SNP loci retained for genotype-aware association testing.",
    "Haplotype blocks": "Local read-backed SNP blocks retained for haplotype association testing.",
}


COLUMN_DEFINITIONS = {
    "zt_label": "Human-readable transcript code used in transcript-level outputs and BAM tags.",
    "code": "Transcript code or grouping identifier written by the relevant upstream stage.",
    "gtf_gene_name": "Reference gene name matched to the assembled transcript when available.",
    "gene_name": "Gene label assigned to the reported row or site.",
    "gene_names": "Gene labels associated with the reported SNP or molecule row.",
    "gene_ids": "Reference or assembled gene identifiers associated with the reported SNP or molecule row.",
    "gene_index": "Assembly-local gene bucket index assigned during transcript summarization.",
    "transcript_index": "Within-gene transcript rank assigned after transcript sorting.",
    "metagene_index": "Overlapping gene group index used for ZN partitioning.",
    "metagene_indices": "Metagene labels associated with the reported SNP or molecule row.",
    "zn_index": "ZN partition index assigned within a metagene.",
    "metagene_partition_count": "Total number of ZN partitions required for the corresponding metagene.",
    "read_support": "Total reads assigned to the transcript or partition.",
    "exact_chain_reads": "Assigned reads whose intron chain exactly matches the retained transcript model.",
    "trunc_assigned_reads": "Reads assigned to the transcript after suffix-compatible truncation absorption.",
    "anchor_reads": "Exact reads supporting the transcript-specific distal anchor used in support-first assignment.",
    "anchor_frac": "Fraction of reachable suffix-family reads that exactly support the retained canonical transcript.",
    "assignment_mode": "Recorded suffix-family assignment policy used during transcript assembly.",
    "sample": "Sample identifier derived from the BAM filename.",
    "chrom": "Reference chromosome or contig containing the reported feature.",
    "pos1": "1-based genomic coordinate of the reported SNP locus.",
    "mod_start0": "0-based inclusive start coordinate of the reported modification site.",
    "mod_end0": "0-based exclusive end coordinate of the reported modification site.",
    "total_reads": "Total assigned transcript reads for the sample.",
    "n_transcripts": "Number of transcripts detected in the sample.",
    "median_reads_per_tx": "Median assigned read count across transcripts detected in the sample.",
    "input_total_reads": "Total BAM alignments encountered before transcript assignment filtering.",
    "primary_reads": "Reads retained after removing secondary and supplementary alignments.",
    "mapq_reads": "Reads retained after applying the minimum MAPQ filter.",
    "intronic_reads": "Reads retained after the minimum intron-count filter.",
    "tagged_reads": "Reads written to the ZT-tagged BAM after transcript assignment.",
    "assigned_fraction": "Fraction of qualifying reads assigned to a retained transcript.",
    "assigned_reads": "Number of reads contributing to the reported transcript-length summary.",
    "mean_read_length": "Mean assigned read length in nucleotides.",
    "median_read_length": "Median assigned read length in nucleotides.",
    "min_read_length": "Shortest assigned read length in nucleotides.",
    "max_read_length": "Longest assigned read length in nucleotides.",
    "tes": "Transcript end site reported for the retained transcript model.",
    "mod_code": "Modification code reported by modkit or downstream aggregation.",
    "n_sites": "Number of unique genomic modification sites observed for the gene and modification code.",
    "p_value": "Nominal p-value from the reported hypothesis test.",
    "p_adj_bh": "Benjamini-Hochberg false-discovery-rate adjusted p-value.",
    "effect_max_abs_frac_diff": "Maximum absolute difference in pooled modified fraction across tested transcript partitions.",
    "effect_max_abs_tx_frac_diff": "Maximum absolute difference in transcript usage or stoichiometry between tested groups.",
    "effect_abs_delta_mod_frac": "Absolute difference in modified-site rate between the tested allele groups.",
    "weighted_within_tx_effect": "Coverage-weighted within-transcript SNP/mod effect size after transcript conditioning.",
    "classification": "Stage-specific label summarizing the inferred outcome for the reported row.",
    "snp_id": "Canonical SNP identifier in `chrom:pos:ref>alt` format.",
    "mod_site_id": "Canonical modification-site identifier in `chrom:start-end:strand:mod` format.",
    "target_mod_code": "Modification code tested at the reported site.",
    "n_alt_reads": "Number of reads carrying the alternative allele in the tested contingency table.",
    "n_ref_reads": "Number of reads carrying the reference allele in the tested contingency table.",
    "n_reads": "Total reads contributing to the reported test.",
    "n_modified": "Reads classified as modified for the target mod code.",
    "n_not_target": "Reads classified as canonical or as another modification state for the target site.",
    "n_transcripts_tested": "Number of transcript partitions retained in the reported test.",
    "cmh_p_value": "Cochran-Mantel-Haenszel p-value after stratifying by transcript.",
    "cmh_p_adj_bh": "Benjamini-Hochberg adjusted CMH p-value.",
    "complete_reads": "Reads covering every SNP in the reported haplotype block.",
    "support_reads": "Reads overlapping at least one SNP in the reported haplotype block.",
    "n_snps": "Number of SNPs represented in the reported haplotype block.",
    "block_id": "Identifier for the reported haplotype block.",
    "context_key": "Gene- or metagene-aware context key used to restrict genotype and modification joins to the same local feature family.",
    "haplotypes": "Observed allele strings for the retained read-backed haplotypes in the reported block.",
    "alt_frac": "Alternative-allele fraction across all supporting reads at the candidate SNP locus.",
    "total_cov": "Total coverage accumulated across samples at the reported SNP or mod site.",
    "ref_count": "Reference-base support across all reads for the candidate SNP locus.",
    "alt_count": "Alternative-base support across all reads for the candidate SNP locus.",
    "samples_with_alt": "Samples contributing at least one alternative-allele read at the candidate SNP locus.",
    "effect_max_abs_mod_rate_diff": "Maximum absolute difference in target modified fraction across the haplotypes retained in the reported test.",
    "category": "Granular structural category explaining why the isoforms differ in modification at the site (see category definitions).",
    "start0": "0-based inclusive start coordinate of the reported modification site.",
    "end0": "0-based exclusive end coordinate of the reported modification site.",
    "strand": "Genomic strand of the reported feature.",
    "n_tx_tested": "Number of transcript (ZN) partitions retained in the differential test at the site.",
    "hi_ZN": "ZN partition index of the higher-stoichiometry isoform in the classified contrast.",
    "hi_arch": "3' architecture of the higher-stoichiometry isoform versus the gene's longest-3'UTR anchor (IPA / TANDEM_APA / FULL_LENGTH / DISTAL_EXT / REFERENCE / AMBIGUOUS).",
    "hi_frac": "Pooled modified fraction (stoichiometry) of the higher isoform at the site.",
    "lo_ZN": "ZN partition index of the lower-stoichiometry isoform in the classified contrast.",
    "lo_arch": "3' architecture of the lower-stoichiometry isoform versus the gene's longest-3'UTR anchor.",
    "lo_frac": "Pooled modified fraction (stoichiometry) of the lower isoform at the site.",
    "anchor_ZN": "ZN partition index of the anchor (longest-3'UTR) isoform used as the structural reference.",
    "status_hi": "Position status of the site within the high isoform (exonic_terminal / exonic_internal / intronic / absent).",
    "status_lo": "Position status of the site within the low isoform.",
    "status_anchor": "Position status of the site within the anchor (longest) isoform.",
    "jd_hi": "Distance (nt) from the site to the nearest spliced junction in the high isoform.",
    "jd_lo": "Distance (nt) from the site to the nearest spliced junction in the low isoform.",
}


CATEGORY_DEFINITIONS = {
    "IPA_UNIQUE": "High isoform is an IPA (intronic polyadenylation) isoform; the A is exonic-terminal in it but intronic/absent in the longer anchor — the modified A only exists in the mature IPA transcript. Cleavage-dependent, IPA-private.",
    "SPLICED_EXON_UNIQUE": "The A sits in an internal/cassette exon present in the high isoform but spliced out (intronic) or absent in the comparator/anchor.",
    "LAST_EXON_DISTAL_ONLY": "The A is in the anchor's distal 3'UTR but the low (proximal) isoform's cleavage site is upstream of it — a distal-3'UTR-private A.",
    "IPA_SHARED_EJC": "High isoform is IPA; the A is shared (exonic in both) but exonic-terminal in IPA versus exonic-internal in the long anchor — the A gains m6A in IPA because the downstream exon-junction complex is removed. Cleavage-dependent.",
    "SPLICING_EJC": "Shared A, non-IPA: terminalized, or a junction within the EJC window in the low/anchor is removed in the high isoform — EJC relief.",
    "LAST_EXON_PROXIMAL_APA_FAVORED": "Same terminal exon (same acceptor), different cleavage site; the PROXIMAL (shorter 3'UTR) isoform carries more m6A. Tandem 3'UTR APA, cleavage-independent.",
    "LAST_EXON_DISTAL_APA_FAVORED": "Same terminal-exon geometry; the DISTAL (longer 3'UTR) isoform carries more m6A. Tandem 3'UTR APA.",
    "ALTERNATIVE_LAST_EXON": "High and low isoforms are both terminal at the site, but their terminal exons begin at different (nearby/overlapping) acceptors — alternative last exons.",
    "INTERGENIC_TERMINAL_EXON": "The site is exonic-terminal in a non-IPA high isoform whose terminal exon is genomically disjoint from the comparator's and separated by a large gap (>= --intergenic-gap) — spatially separated / read-through / intergenic-scale alternative last exons.",
    "SHARED_TERMINAL_EXON": "High and low isoforms share the same terminal exon AND the same cleavage site; m6A tracks isoform identity, not APA or EJC.",
    "SHARED_INTERNAL_EXON": "The site is in a constitutive internal exon with no junction asymmetry — not attributable to 3' architecture.",
    "UNEXPLAINED_SHARED": "Rare residual: terminal in the high isoform, internal in the low, with no nearby differential junction.",
    "HI_INTRONIC_ARTIFACT": "The high-m6A isoform does not structurally contain the A (intronic/absent) — the 'high' stoichiometry is intron-read noise, not a real isoform-specific site.",
    "UNCLASSIFIED": "Fewer than two covered isoform models at the site, or no anchor isoform — cannot be assigned a structural category.",
}


def parse_args():
    ap = argparse.ArgumentParser(description="Generate a lightweight HTML report for modulator outputs.")
    ap.add_argument("--classification", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--tx-counts", required=True)
    ap.add_argument("--pca-png", required=True)
    ap.add_argument("--sample-stats", required=True)
    ap.add_argument("--read-stats", required=True)
    ap.add_argument("--tx-lengths", required=True)
    ap.add_argument("--partition-map", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--zn-long", default="")
    ap.add_argument("--zt-long", default="")
    ap.add_argument("--diff-results", default="")
    ap.add_argument("--diff-figs-dir", default="")
    ap.add_argument("--classified-sites", default="",
                    help="{prefix}__ZN_site_classified.tsv from classify_diff_sites.py")
    ap.add_argument("--class-figs-dir", default="",
                    help="{prefix}__figs_by_category directory (one subdir per category) "
                         "holding the 2-panel per-sample stoichiometry / pooled-coverage figures")
    ap.add_argument("--arch-figs-dir", default="",
                    help="{prefix}__figs_by_category_arch directory (one subdir per category) "
                         "holding the isoform architecture-map (exon/intron locus-track) figures, "
                         "featured as the primary per-category figure")
    ap.add_argument("--max-class-figs-per-category", type=int, default=10,
                    help="max per-category figures to embed in the report. Default 10.")
    ap.add_argument("--multigene-summary-glob", default="")
    ap.add_argument("--candidate-snps", default="")
    ap.add_argument("--snp-tx-assoc", default="")
    ap.add_argument("--snp-mod-assoc", default="")
    ap.add_argument("--snp-tx-mod-assoc", default="")
    ap.add_argument("--hap-blocks", default="")
    ap.add_argument("--hap-tx-assoc", default="")
    ap.add_argument("--hap-mod-assoc", default="")
    ap.add_argument("--max-diff-figs", type=int, default=6)
    ap.add_argument("--top-transcripts", type=int, default=20)
    ap.add_argument("--top-genes", type=int, default=20)
    return ap.parse_args()


def clean_columns(df):
    df.columns = [str(c).lstrip("#") for c in df.columns]
    return df


def read_tsv(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    return clean_columns(pd.read_csv(path, sep="\t"))


def read_summary_metrics(paths):
    rows = []
    for path in sorted(paths):
        sample = os.path.basename(path).split(".")[0]
        metrics = {}
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    break
                if line == "metric\tvalue":
                    continue
                key, value = line.split("\t", 1)
                metrics[key] = value
        metrics["sample"] = sample
        rows.append(metrics)
    return pd.DataFrame(rows)


def embed_png(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return ""
    with open(path, "rb") as fh:
        encoded = base64.b64encode(fh.read()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def category_distribution_png(counts):
    """Horizontal bar chart of classified-site counts per category. Returns a
    base64 data URI (or "" if matplotlib/data unavailable)."""
    if not counts:
        return ""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from io import BytesIO
    except Exception:
        return ""
    items = sorted(counts.items(), key=lambda kv: kv[1])  # ascending -> largest on top
    labels = [k for k, _ in items]
    values = [v for _, v in items]
    total = sum(values) or 1
    height = max(2.4, 0.42 * len(labels) + 1.1)
    fig, ax = plt.subplots(figsize=(8.6, height))
    ypos = range(len(labels))
    bars = ax.barh(list(ypos), values, color="#c98a5e", edgecolor="#7d3c1f", linewidth=0.9)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Significant m6A sites")
    ax.set_title("Differential m6A sites by structural category")
    pad = max(values) * 0.01 if values else 0.1
    for rect, val in zip(bars, values):
        ax.text(rect.get_width() + pad, rect.get_y() + rect.get_height() / 2.0,
                f"{val:,} ({100.0 * val / total:.1f}%)", va="center", ha="left",
                fontsize=8, color="#4a3f35")
    ax.set_xlim(0, max(values) * 1.18 if values else 1)
    ax.grid(True, axis="x", linestyle="--", linewidth=0.6, alpha=0.4)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def clickable_image_html(src, alt, *, caption="", figure_class="image-card"):
    escaped_alt = html.escape(str(alt))
    caption_html = f"<figcaption>{html.escape(str(caption))}</figcaption>" if caption else ""
    return (
        f"<figure class='{figure_class}'>"
        f"<a class='image-link' href='{src}' target='_blank' rel='noopener noreferrer' title='Open full-size image in a new tab'>"
        f"<img src='{src}' alt='{escaped_alt}' />"
        "<span class='expand-badge' aria-hidden='true'>↗</span>"
        "</a>"
        f"{caption_html}"
        "</figure>"
    )


def fmt_int(val):
    try:
        return f"{int(val):,}"
    except Exception:
        return "0"


def visible_definitions(labels):
    return [(label, CARD_DEFINITIONS[label]) for label in labels if label in CARD_DEFINITIONS]


def column_definitions(columns):
    defs = []
    seen = set()
    for col in columns:
        if col in seen:
            continue
        seen.add(col)
        defs.append((col, COLUMN_DEFINITIONS.get(col, f"Value carried through from the upstream `{col}` field.")))
    return defs


def definitions_html(items, *, summary="Definitions", open_by_default=False):
    if not items:
        return ""
    open_attr = " open" if open_by_default else ""
    defs = "".join(
        f"<dt>{html.escape(str(term))}</dt><dd>{html.escape(str(desc))}</dd>"
        for term, desc in items
    )
    return (
        f"<details class='definitions'{open_attr}>"
        f"<summary>{html.escape(summary)}</summary>"
        f"<dl>{defs}</dl>"
        "</details>"
    )


def df_to_html(df, max_rows=25):
    if df is None or df.empty:
        return "<p class='muted'>No data available.</p>"
    table = df.head(max_rows).to_html(index=False, escape=True, classes="datatable", border=0)
    return f"<div class='table-wrap'>{table}</div>"


def section(title, body, *, intro="", definitions=""):
    intro_html = f"<p class='section-intro'>{html.escape(intro)}</p>" if intro else ""
    return f"<section><h2>{html.escape(title)}</h2>{intro_html}{definitions}{body}</section>"


def subsection(title, body, *, definitions=""):
    return f"<div class='subsection'><h3>{html.escape(title)}</h3>{definitions}{body}</div>"


def category_figure_gallery(figs_dir, category, max_figs,
                            summary_label="site figure(s) — per-sample stoichiometry &amp; pooled coverage",
                            open_by_default=False):
    """Collapsible gallery of the per-category figures under figs_dir/<category>/*.png."""
    if not figs_dir or max_figs <= 0:
        return ""
    cat_dir = os.path.join(figs_dir, category)
    if not os.path.isdir(cat_dir):
        return ""
    fig_paths = sorted(glob.glob(os.path.join(cat_dir, "*.png")))[:max_figs]
    if not fig_paths:
        return ""
    pieces = []
    for path in fig_paths:
        img = embed_png(path)
        if img:
            pieces.append(clickable_image_html(img, os.path.basename(path), caption=os.path.basename(path)))
    if not pieces:
        return ""
    gallery = "<div class='gallery'>" + "".join(pieces) + "</div>"
    open_attr = " open" if open_by_default else ""
    return (
        f"<details class='definitions'{open_attr}>"
        f"<summary>Top {len(pieces)} {summary_label}</summary>"
        f"{gallery}</details>"
    )


def build_classification_section(class_df, class_figs_dir, arch_figs_dir, max_figs_per_category, top_n):
    """Site-classification overview: per-category counts, the distribution graph,
    and per-category top sites. Each category features the isoform architecture-map
    (locus-track) figures as the primary visual, with the 2-panel per-sample
    stoichiometry figures kept as a secondary collapsible gallery."""
    if class_df is None or class_df.empty or "category" not in class_df.columns:
        return section(
            "Site Classification",
            "<p class='muted'>No classified differential sites available.</p>",
            intro="Granular structural classification of significant between-transcript modification sites.",
        )

    counts = class_df["category"].value_counts().to_dict()
    total = int(sum(counts.values())) or 1

    dist_uri = category_distribution_png(counts)
    dist_html = clickable_image_html(dist_uri, "Category distribution", figure_class="hero-figure") if dist_uri else ""

    summary_df = (
        pd.DataFrame({"category": list(counts.keys()), "n_sites": list(counts.values())})
        .sort_values("n_sites", ascending=False)
        .reset_index(drop=True)
    )
    summary_df["pct"] = (100.0 * summary_df["n_sites"] / total).map(lambda v: f"{v:.1f}%")

    present = list(summary_df["category"])
    cat_defs = [(c, CATEGORY_DEFINITIONS.get(c, "Structural category assigned by classify_diff_sites.")) for c in present]

    overview = (
        "<div class='overview-layout'>"
        f"<div>{df_to_html(summary_df, max_rows=len(summary_df))}"
        f"{definitions_html(cat_defs, summary='Category definitions', open_by_default=False)}</div>"
        f"<div class='hero'>{dist_html or '<p class=\"muted\">Distribution graph unavailable.</p>'}</div>"
        "</div>"
    )

    detail_cols = [
        c for c in [
            "gene_name", "mod_code", "chrom", "start0", "strand", "hi_ZN", "hi_arch", "hi_frac",
            "lo_ZN", "lo_arch", "lo_frac", "effect_max_abs_frac_diff", "p_adj_bh",
        ] if c in class_df.columns
    ]
    sort_cols = [c for c in ["effect_max_abs_frac_diff"] if c in class_df.columns]

    blocks = []
    for cat in present:
        sub = class_df[class_df["category"] == cat]
        if sort_cols:
            sub = sub.sort_values(sort_cols, ascending=False)
        n = len(sub)
        table_html = df_to_html(sub[detail_cols] if detail_cols else sub, max_rows=top_n)
        # Architecture-map (locus-track) figures are the PRIMARY visual (open by default);
        # the 2-panel per-sample stoichiometry figures are kept as a secondary gallery.
        arch_gallery_html = category_figure_gallery(
            arch_figs_dir, cat, max_figs_per_category,
            summary_label="isoform architecture map(s) — exon/intron locus tracks with the site marked",
            open_by_default=True,
        )
        stoich_gallery_html = category_figure_gallery(class_figs_dir, cat, max_figs_per_category)
        defn = CATEGORY_DEFINITIONS.get(cat, "")
        defn_html = f"<p class='section-intro'>{html.escape(defn)}</p>" if defn else ""
        blocks.append(
            subsection(
                f"{cat} — top {min(top_n, n)} of {n} site(s)",
                defn_html + table_html + arch_gallery_html + stoich_gallery_html,
            )
        )

    body = overview + "".join(blocks)
    return section(
        "Site Classification",
        body,
        intro=(
            "Every significant between-transcript modification site (across all detected "
            "mod_codes; BH-FDR and the &gt;10% absolute stoichiometry rule) is assigned one "
            "MECE structural category explaining why the isoforms differ, anchored to the "
            "gene's longest-3'UTR isoform. Each category lists its top sites by effect size, "
            "featuring an isoform architecture map (every isoform drawn as exon/intron tracks "
            "with the modified site marked) plus the per-sample stoichiometry / pooled-coverage "
            "figures."
        ),
        definitions=definitions_html(column_definitions(detail_cols), summary="Column definitions") if detail_cols else "",
    )


def main():
    args = parse_args()

    class_df = read_tsv(args.classification)
    metrics_df = read_tsv(args.metrics)
    tx_counts_df = read_tsv(args.tx_counts)
    sample_stats_df = read_tsv(args.sample_stats)
    read_stats_df = read_tsv(args.read_stats)
    tx_lengths_df = read_tsv(args.tx_lengths)
    partition_map_df = read_tsv(args.partition_map)
    zn_long_df = read_tsv(args.zn_long)
    zt_long_df = read_tsv(args.zt_long)
    diff_df = read_tsv(args.diff_results)
    classified_df = read_tsv(args.classified_sites)
    overlap_df = read_summary_metrics(glob.glob(args.multigene_summary_glob)) if args.multigene_summary_glob else pd.DataFrame()
    candidate_snps_df = read_tsv(args.candidate_snps)
    snp_tx_assoc_df = read_tsv(args.snp_tx_assoc)
    snp_mod_assoc_df = read_tsv(args.snp_mod_assoc)
    snp_tx_mod_assoc_df = read_tsv(args.snp_tx_mod_assoc)
    hap_blocks_df = read_tsv(args.hap_blocks)
    hap_tx_assoc_df = read_tsv(args.hap_tx_assoc)
    hap_mod_assoc_df = read_tsv(args.hap_mod_assoc)

    overview_cards = []
    n_tx = len(class_df)
    n_genes = class_df["gene_index"].nunique() if "gene_index" in class_df.columns and not class_df.empty else 0
    n_metagenes = class_df["metagene_index"].nunique() if "metagene_index" in class_df.columns and not class_df.empty else 0
    max_partitions = class_df["metagene_partition_count"].max() if "metagene_partition_count" in class_df.columns and not class_df.empty else 0
    total_reads = class_df["read_support"].astype(float).sum() if "read_support" in class_df.columns and not class_df.empty else 0
    trunc_reads = class_df["trunc_assigned_reads"].astype(float).sum() if "trunc_assigned_reads" in class_df.columns and not class_df.empty else 0
    exact_reads = class_df["exact_chain_reads"].astype(float).sum() if "exact_chain_reads" in class_df.columns and not class_df.empty else 0
    overview_cards.extend([
        ("Assembled transcripts", fmt_int(n_tx)),
        ("Gene buckets", fmt_int(n_genes)),
        ("Metagenes", fmt_int(n_metagenes)),
        ("Max ZN partitions / metagene", fmt_int(max_partitions)),
        ("Assigned reads", fmt_int(total_reads)),
        ("Exact-chain reads", fmt_int(exact_reads)),
        ("Truncation-assigned reads", fmt_int(trunc_reads)),
    ])
    if not candidate_snps_df.empty:
        overview_cards.extend([
            ("Segregating SNPs", fmt_int(len(candidate_snps_df))),
            ("Haplotype blocks", fmt_int(len(hap_blocks_df) if not hap_blocks_df.empty else 0)),
        ])

    pca_img = embed_png(args.pca_png)
    pca_html = clickable_image_html(pca_img, "PCA plot", figure_class="hero-figure") if pca_img else "<p class='muted'>PCA plot unavailable.</p>"

    top_tx_cols = [
        c for c in [
            "zt_label", "gtf_gene_name", "read_support", "exact_chain_reads",
            "trunc_assigned_reads", "anchor_reads", "anchor_frac",
            "metagene_index", "zn_index", "classification"
        ] if c in class_df.columns
    ]
    top_tx_df = class_df.sort_values(["read_support", "exact_chain_reads"], ascending=False) if not class_df.empty else class_df

    top_gene_sites_df = pd.DataFrame()
    if not zn_long_df.empty and "gene_name" in zn_long_df.columns:
        top_gene_sites_df = (
            zn_long_df.assign(site_key=zn_long_df[["chrom", "start0", "end0", "strand", "mod_code"]].astype(str).agg(":".join, axis=1))
            .groupby(["gene_name", "mod_code"], as_index=False)["site_key"].nunique()
            .rename(columns={"site_key": "n_sites"})
            .sort_values("n_sites", ascending=False)
        )

    diff_html = "<p class='muted'>No differential-site results available.</p>"
    diff_cols = []
    if not diff_df.empty:
        diff_cols = [
            c for c in [
                "gene_name", "mod_code", "chrom", "start0", "end0", "strand",
                "p_value", "p_adj_bh", "effect_max_abs_frac_diff"
            ] if c in diff_df.columns
        ]
        diff_html = df_to_html(
            diff_df[diff_cols].sort_values(["p_adj_bh", "effect_max_abs_frac_diff"], ascending=[True, False]),
            max_rows=20,
        )

    diff_fig_html = ""
    if args.diff_figs_dir and os.path.isdir(args.diff_figs_dir):
        fig_paths = sorted(glob.glob(os.path.join(args.diff_figs_dir, "*.png")))[: args.max_diff_figs]
        if fig_paths:
            pieces = []
            for path in fig_paths:
                img = embed_png(path)
                if img:
                    pieces.append(clickable_image_html(img, os.path.basename(path), caption=os.path.basename(path)))
            diff_fig_html = "<div class='gallery'>" + "".join(pieces) + "</div>"

    read_stats_counts_cols = [
        c for c in [
            "sample", "total_reads_bam", "total_mapped", "total_unmapped", "considered_reads",
            "failed_unmapped", "failed_secondary_or_supp", "failed_low_mapq", "failed_low_introns",
            "failed_low_softclip3p", "zt_tagged_exists", "zt_total_records", "zt_unmapped_records",
            "zt_mapped_records", "assigned_reads", "zt_mapped_unassigned_reads"
        ] if c in read_stats_df.columns
    ]
    read_stats_length_cols = ["sample"] + [c for c in read_stats_df.columns if "_len_" in c]
    read_stats_length_cols = [c for c in read_stats_length_cols if c in read_stats_df.columns]

    candidate_snps_df_view = pd.DataFrame()
    if not candidate_snps_df.empty:
        keep = [
            c for c in [
                "snp_id", "gene_names", "metagene_indices", "total_cov",
                "ref_count", "alt_count", "alt_frac", "samples_with_alt"
            ] if c in candidate_snps_df.columns
        ]
        candidate_snps_df_view = candidate_snps_df[keep].sort_values(["alt_count", "alt_frac"], ascending=False)

    snp_tx_df_view = pd.DataFrame()
    if not snp_tx_assoc_df.empty:
        keep = [
            c for c in [
                "snp_id", "gene_names", "metagene_indices", "n_alt_reads",
                "n_transcripts_tested", "p_value", "p_adj_bh", "effect_max_abs_tx_frac_diff"
            ] if c in snp_tx_assoc_df.columns
        ]
        snp_tx_df_view = snp_tx_assoc_df[keep].sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False])

    snp_mod_df_view = pd.DataFrame()
    if not snp_mod_assoc_df.empty:
        keep = [
            c for c in [
                "snp_id", "mod_site_id", "target_mod_code", "gene_names", "n_alt_reads",
                "n_modified", "p_value", "p_adj_bh", "effect_abs_delta_mod_frac"
            ] if c in snp_mod_assoc_df.columns
        ]
        snp_mod_df_view = snp_mod_assoc_df[keep].sort_values(["p_adj_bh", "effect_abs_delta_mod_frac"], ascending=[True, False])

    joint_df_view = pd.DataFrame()
    if not snp_tx_mod_assoc_df.empty:
        keep = [
            c for c in [
                "snp_id", "mod_site_id", "target_mod_code", "n_transcripts_tested",
                "cmh_p_value", "cmh_p_adj_bh", "weighted_within_tx_effect", "classification"
            ] if c in snp_tx_mod_assoc_df.columns
        ]
        joint_df_view = snp_tx_mod_assoc_df[keep].sort_values(["cmh_p_adj_bh", "weighted_within_tx_effect"], ascending=[True, False])

    hap_blocks_df_view = pd.DataFrame()
    if not hap_blocks_df.empty:
        keep = [c for c in ["block_id", "context_key", "chrom", "n_snps", "support_reads", "complete_reads", "haplotypes"] if c in hap_blocks_df.columns]
        hap_blocks_df_view = hap_blocks_df[keep].sort_values(["complete_reads", "n_snps"], ascending=False)

    hap_tx_df_view = pd.DataFrame()
    if not hap_tx_assoc_df.empty:
        keep = [c for c in ["block_id", "n_transcripts_tested", "p_value", "p_adj_bh", "effect_max_abs_tx_frac_diff"] if c in hap_tx_assoc_df.columns]
        hap_tx_df_view = hap_tx_assoc_df[keep].sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False])

    hap_mod_df_view = pd.DataFrame()
    if not hap_mod_assoc_df.empty:
        keep = [c for c in ["block_id", "mod_site_id", "target_mod_code", "p_value", "p_adj_bh", "effect_max_abs_mod_rate_diff"] if c in hap_mod_assoc_df.columns]
        hap_mod_df_view = hap_mod_assoc_df[keep].sort_values(["p_adj_bh", "effect_max_abs_mod_rate_diff"], ascending=[True, False])

    cards_html = "".join(
        f"<div class='card'><div class='label'>{html.escape(label)}</div><div class='value'>{html.escape(value)}</div></div>"
        for label, value in overview_cards
    )

    overview_defs = definitions_html(
        visible_definitions([label for label, _ in overview_cards]),
        summary="Definitions for the overview metrics",
        open_by_default=True,
    )

    body = []
    body.append(
        section(
            "Overview",
            (
                "<div class='overview-layout'>"
                "<div>"
                f"<div class='cards'>{cards_html}</div>{overview_defs}"
                "</div>"
                f"<div class='hero'>{pca_html}</div>"
                "</div>"
            ),
            intro="Top-level counts summarize the final assembled and downstream-tested objects written by the current modulator run.",
        )
    )
    body.append(
        section(
            "Top Transcript Partitions",
            df_to_html(top_tx_df[top_tx_cols], max_rows=args.top_transcripts) if top_tx_cols else "<p class='muted'>No transcript summary columns available.</p>",
            intro="Highest-support transcript models retained after assembly and reference classification.",
            definitions=definitions_html(column_definitions(top_tx_cols), summary="Column definitions"),
        )
    )
    body.append(
        section(
            "Sample Stats",
            df_to_html(sample_stats_df, max_rows=100),
            intro="Per-sample transcript assignment totals and detection breadth derived from the final retained transcript models.",
            definitions=definitions_html(column_definitions(list(sample_stats_df.columns)), summary="Column definitions"),
        )
    )
    body.append(
        section(
            "Read Funnel",
            (
                subsection(
                    "Read Retention Counts",
                    df_to_html(read_stats_df[read_stats_counts_cols], max_rows=100) if read_stats_counts_cols else "<p class='muted'>No read-count funnel columns available.</p>",
                    definitions=definitions_html(column_definitions(read_stats_counts_cols), summary="Column definitions") if read_stats_counts_cols else "",
                ) +
                subsection(
                    "Read Length Summaries",
                    df_to_html(read_stats_df[read_stats_length_cols], max_rows=100) if len(read_stats_length_cols) > 1 else "<p class='muted'>No read-length summary columns available.</p>",
                    definitions=definitions_html(column_definitions(read_stats_length_cols), summary="Column definitions") if len(read_stats_length_cols) > 1 else "",
                )
            ),
            intro="Sample-level retention across the primary read filters and transcript-assignment workflow.",
        )
    )
    body.append(
        section(
            "Partition Map",
            df_to_html(partition_map_df, max_rows=args.top_transcripts),
            intro="Mapping between human-readable transcript labels and the gene, metagene, and ZN identifiers used downstream.",
            definitions=definitions_html(column_definitions(list(partition_map_df.columns)), summary="Column definitions"),
        )
    )
    body.append(
        section(
            "Assigned Read Lengths",
            df_to_html(tx_lengths_df, max_rows=args.top_transcripts),
            intro="Assigned-read length summaries for the most supported retained transcript models.",
            definitions=definitions_html(column_definitions(list(tx_lengths_df.columns)), summary="Column definitions"),
        )
    )
    body.append(
        section(
            "Number of Modified Sites Per Gene (ZN)",
            df_to_html(top_gene_sites_df, max_rows=args.top_genes),
            intro="Counts of unique genomic modification sites observed per gene and modification code in the ZN aggregation.",
            definitions=definitions_html(column_definitions(list(top_gene_sites_df.columns)), summary="Column definitions"),
        )
    )
    body.append(
        section(
            "Overlap Resolution",
            df_to_html(overlap_df.fillna("0"), max_rows=100) if not overlap_df.empty else "<p class='muted'>No overlap-resolution summaries available.</p>",
            intro="Multigene-overlap outcomes from the optional read-resolution stage.",
            definitions=definitions_html(column_definitions(list(overlap_df.columns)), summary="Column definitions") if not overlap_df.empty else "",
        )
    )
    body.append(
        section(
            "Differential Sites",
            diff_html + diff_fig_html,
            intro="ZN sites where pooled transcript partitions differ in modified fraction after applying the ZN site filter.",
            definitions=definitions_html(column_definitions(diff_cols), summary="Result-column definitions") if diff_cols else "",
        )
    )
    body.append(
        build_classification_section(
            classified_df,
            args.class_figs_dir,
            args.arch_figs_dir,
            args.max_class_figs_per_category,
            args.max_class_figs_per_category,
        )
    )
    body.append(
        section(
            "Segregating SNP Candidates",
            df_to_html(candidate_snps_df_view, max_rows=args.top_genes) if not candidate_snps_df_view.empty else "<p class='muted'>No segregating SNP candidates available.</p>",
            intro="Read-supported non-reference loci discovered in the cleaned tagged BAMs.",
            definitions=definitions_html(column_definitions(list(candidate_snps_df_view.columns)), summary="Column definitions") if not candidate_snps_df_view.empty else "",
        )
    )
    body.append(
        section(
            "SNP to Transcript Associations",
            df_to_html(snp_tx_df_view, max_rows=args.top_genes) if not snp_tx_df_view.empty else "<p class='muted'>No SNP to transcript associations available.</p>",
            intro="Associations between segregating SNP alleles and transcript-partition usage.",
            definitions=definitions_html(column_definitions(list(snp_tx_df_view.columns)), summary="Column definitions") if not snp_tx_df_view.empty else "",
        )
    )
    body.append(
        section(
            "SNP to Epitranscriptome Associations",
            df_to_html(snp_mod_df_view, max_rows=args.top_genes) if not snp_mod_df_view.empty else "<p class='muted'>No SNP to epitranscriptome associations available.</p>",
            intro="Associations between segregating SNP alleles and target modification states on the same molecules.",
            definitions=definitions_html(column_definitions(list(snp_mod_df_view.columns)), summary="Column definitions") if not snp_mod_df_view.empty else "",
        )
    )
    body.append(
        section(
            "SNP Transcript Epitranscriptome Dependency",
            df_to_html(joint_df_view, max_rows=args.top_genes) if not joint_df_view.empty else "<p class='muted'>No joint SNP-transcript-epitranscriptome dependency results available.</p>",
            intro="Transcript-conditioned SNP/mod tests that distinguish direct epitranscriptome effects from transcript-composition shifts.",
            definitions=definitions_html(column_definitions(list(joint_df_view.columns)), summary="Column definitions") if not joint_df_view.empty else "",
        )
    )
    body.append(
        section(
            "Haplotype Blocks",
            df_to_html(hap_blocks_df_view, max_rows=args.top_genes) if not hap_blocks_df_view.empty else "<p class='muted'>No haplotype blocks available.</p>",
            intro="Local read-backed SNP blocks retained for haplotype association testing.",
            definitions=definitions_html(column_definitions(list(hap_blocks_df_view.columns)), summary="Column definitions") if not hap_blocks_df_view.empty else "",
        )
    )

    hap_sections = []
    if not hap_tx_df_view.empty:
        hap_sections.append(
            subsection(
                "Haplotype to Transcript",
                df_to_html(hap_tx_df_view, max_rows=args.top_genes),
                definitions=definitions_html(column_definitions(list(hap_tx_df_view.columns)), summary="Column definitions"),
            )
        )
    if not hap_mod_df_view.empty:
        hap_sections.append(
            subsection(
                "Haplotype to Epitranscriptome",
                df_to_html(hap_mod_df_view, max_rows=args.top_genes),
                definitions=definitions_html(column_definitions(list(hap_mod_df_view.columns)), summary="Column definitions"),
            )
        )
    body.append(
        section(
            "Haplotype Associations",
            "".join(hap_sections) if hap_sections else "<p class='muted'>No haplotype associations available.</p>",
            intro="Associations between local haplotype blocks and transcript or modification outcomes.",
        )
    )

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(args.title)}</title>
  <style>
    :root {{
      --bg: #f4efe7;
      --ink: #1d1b19;
      --card: #fffaf2;
      --line: #d8c9b3;
      --accent: #7d3c1f;
      --accent-soft: #f1dcc5;
      --panel: #fffdf8;
    }}
    body {{
      margin: 0;
      font-family: Georgia, "Times New Roman", serif;
      background: linear-gradient(180deg, #efe5d8 0%, var(--bg) 100%);
      color: var(--ink);
    }}
    main {{
      max-width: 1380px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2, h3 {{
      font-weight: 700;
      letter-spacing: 0.01em;
      margin-top: 0;
    }}
    header {{
      margin-bottom: 20px;
    }}
    section {{
      background: var(--card);
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 20px;
      margin-bottom: 20px;
      box-shadow: 0 10px 30px rgba(70, 44, 20, 0.08);
      overflow: hidden;
    }}
    .section-intro {{
      margin: 0 0 14px 0;
      color: #4f453a;
      line-height: 1.45;
    }}
    .overview-layout {{
      display: grid;
      grid-template-columns: minmax(0, 1.2fr) minmax(320px, 1fr);
      gap: 20px;
      align-items: start;
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
      gap: 12px;
      margin-bottom: 16px;
    }}
    .card {{
      background: var(--accent-soft);
      border-radius: 14px;
      padding: 14px;
      min-height: 92px;
    }}
    .label {{
      font-size: 0.85rem;
      color: #5c4d3d;
      margin-bottom: 6px;
      line-height: 1.3;
    }}
    .value {{
      font-size: 1.5rem;
      color: var(--accent);
      font-weight: 700;
    }}
    .hero {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 14px;
    }}
    .hero-figure, .image-card {{
      margin: 0;
    }}
    .image-link {{
      position: relative;
      display: block;
      text-decoration: none;
      color: inherit;
      cursor: zoom-in;
    }}
    img {{
      max-width: 100%;
      border-radius: 12px;
      border: 1px solid var(--line);
      background: white;
      display: block;
    }}
    .image-link img {{
      transition: transform 0.16s ease, box-shadow 0.16s ease;
    }}
    .image-link:hover img,
    .image-link:focus-visible img {{
      transform: translateY(-2px);
      box-shadow: 0 14px 28px rgba(45, 29, 16, 0.16);
    }}
    .expand-badge {{
      position: absolute;
      top: 12px;
      right: 12px;
      width: 34px;
      height: 34px;
      border-radius: 999px;
      display: grid;
      place-items: center;
      background: rgba(29, 27, 25, 0.78);
      color: white;
      font-size: 1rem;
      font-weight: 700;
      opacity: 0;
      transform: scale(0.92);
      transition: opacity 0.16s ease, transform 0.16s ease;
      pointer-events: none;
    }}
    .image-link:hover .expand-badge,
    .image-link:focus-visible .expand-badge {{
      opacity: 1;
      transform: scale(1);
    }}
    .definitions {{
      margin: 0 0 14px 0;
      border: 1px solid var(--line);
      border-radius: 12px;
      background: rgba(255, 255, 255, 0.72);
    }}
    .definitions summary {{
      cursor: pointer;
      padding: 10px 14px;
      font-weight: 600;
    }}
    .definitions dl {{
      margin: 0;
      padding: 0 14px 14px 14px;
      display: grid;
      grid-template-columns: minmax(180px, 240px) minmax(0, 1fr);
      gap: 8px 16px;
    }}
    .definitions dt {{
      font-weight: 700;
      color: #4a2e1d;
    }}
    .definitions dd {{
      margin: 0;
      color: #4f453a;
      line-height: 1.4;
    }}
    .table-wrap {{
      overflow: auto;
      max-height: 28rem;
      border: 1px solid var(--line);
      border-radius: 14px;
      background: white;
    }}
    .datatable {{
      width: 100%;
      min-width: 760px;
      border-collapse: collapse;
      font-size: 0.88rem;
    }}
    .datatable th, .datatable td {{
      border-bottom: 1px solid var(--line);
      padding: 8px 10px;
      text-align: left;
      vertical-align: top;
      white-space: normal;
      word-break: break-word;
    }}
    .datatable th {{
      background: #f7eddc;
      position: sticky;
      top: 0;
      z-index: 1;
    }}
    .muted {{
      color: #6d6255;
    }}
    .gallery {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 16px;
      margin-top: 16px;
    }}
    figure {{
      margin: 0;
    }}
    figcaption {{
      margin-top: 8px;
      font-size: 0.82rem;
      color: #5c4d3d;
      word-break: break-word;
    }}
    .subsection + .subsection {{
      margin-top: 18px;
    }}
    @media (max-width: 980px) {{
      main {{
        padding: 16px;
      }}
      .overview-layout {{
        grid-template-columns: 1fr;
      }}
      .definitions dl {{
        grid-template-columns: 1fr;
      }}
      .datatable {{
        min-width: 620px;
      }}
    }}
  </style>
</head>
<body>
  <main>
    <header>
      <h1>{html.escape(args.title)}</h1>
      <p class="muted">Generated from the supplied modulator outputs in the current run context.</p>
    </header>
    {''.join(body)}
  </main>
</body>
</html>
"""

    os.makedirs(os.path.dirname(args.out_html) or ".", exist_ok=True)
    with open(args.out_html, "w") as out:
        out.write(html_doc)


if __name__ == "__main__":
    main()
