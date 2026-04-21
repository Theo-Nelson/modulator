#!/usr/bin/env python3

import argparse
import base64
import glob
import html
import os
from pathlib import Path

import pandas as pd


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


def fmt_int(val):
    try:
        return f"{int(val):,}"
    except Exception:
        return "0"


def fmt_float(val, digits=3):
    try:
        return f"{float(val):.{digits}f}"
    except Exception:
        return "0"


def df_to_html(df, max_rows=25):
    if df is None or df.empty:
        return "<p class='muted'>No data available.</p>"
    return df.head(max_rows).to_html(index=False, escape=True, classes="datatable", border=0)


def section(title, body):
    return f"<section><h2>{html.escape(title)}</h2>{body}</section>"


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
            ("Haplotype Blocks", fmt_int(len(hap_blocks_df) if not hap_blocks_df.empty else 0)),
        ])

    pca_img = embed_png(args.pca_png)
    pca_html = f"<img src='{pca_img}' alt='PCA plot' />" if pca_img else "<p class='muted'>PCA plot unavailable.</p>"

    top_tx_cols = [
        c for c in [
            "zt_label", "gene_name", "read_support", "exact_chain_reads",
            "trunc_assigned_reads", "anchor_reads", "anchor_frac",
            "metagene_index", "zn_index", "classification"
        ] if c in class_df.columns
    ]
    top_tx_df = class_df.sort_values(["read_support", "exact_chain_reads"], ascending=False) if not class_df.empty else class_df

    top_gene_sites_html = "<p class='muted'>No ZN site aggregation available.</p>"
    if not zn_long_df.empty and "gene_name" in zn_long_df.columns:
        gene_sites = (
            zn_long_df.assign(site_key=zn_long_df[["chrom", "start0", "end0", "strand", "mod_code"]].astype(str).agg(":".join, axis=1))
            .groupby(["gene_name", "mod_code"], as_index=False)["site_key"].nunique()
            .rename(columns={"site_key": "n_sites"})
            .sort_values("n_sites", ascending=False)
        )
        top_gene_sites_html = df_to_html(gene_sites, max_rows=args.top_genes)

    overlap_html = "<p class='muted'>No overlap-resolution summaries available.</p>"
    if not overlap_df.empty:
        overlap_html = df_to_html(overlap_df.fillna("0"), max_rows=100)

    diff_html = "<p class='muted'>No differential-site results available.</p>"
    if not diff_df.empty:
        keep = [c for c in [
            "gene_name", "mod_code", "chrom", "start0", "end0", "strand",
            "p_value", "p_adj_bh", "effect_max_abs_frac_diff",
            "pvalue", "fdr", "effect_size"
        ] if c in diff_df.columns]
        sort_col = next(
            (c for c in ["effect_max_abs_frac_diff", "effect_size", "p_adj_bh", "fdr", "p_value", "pvalue"] if c in diff_df.columns),
            keep[0] if keep else None,
        )
        diff_html = df_to_html(
            diff_df[keep].sort_values(sort_col, ascending=(sort_col in {"p_adj_bh", "fdr", "p_value", "pvalue"}))
            if keep and sort_col else diff_df,
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
                    pieces.append(
                        "<figure>"
                        f"<img src='{img}' alt='{html.escape(os.path.basename(path))}' />"
                        f"<figcaption>{html.escape(os.path.basename(path))}</figcaption>"
                        "</figure>"
                    )
            diff_fig_html = "<div class='gallery'>" + "".join(pieces) + "</div>"

    candidate_snps_html = "<p class='muted'>No segregating SNP candidates available.</p>"
    if not candidate_snps_df.empty:
        keep = [c for c in ["snp_id", "gene_names", "metagene_indices", "total_cov", "ref_count", "alt_count", "alt_frac", "samples_with_alt"] if c in candidate_snps_df.columns]
        candidate_snps_html = df_to_html(candidate_snps_df[keep].sort_values(["alt_count", "alt_frac"], ascending=False), max_rows=args.top_genes)

    snp_tx_html = "<p class='muted'>No SNP to transcript associations available.</p>"
    if not snp_tx_assoc_df.empty:
        keep = [c for c in ["snp_id", "gene_names", "metagene_indices", "n_alt_reads", "n_transcripts_tested", "p_value", "p_adj_bh", "effect_max_abs_tx_frac_diff"] if c in snp_tx_assoc_df.columns]
        snp_tx_html = df_to_html(snp_tx_assoc_df[keep].sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False]), max_rows=args.top_genes)

    snp_mod_html = "<p class='muted'>No SNP to epitranscriptome associations available.</p>"
    if not snp_mod_assoc_df.empty:
        keep = [c for c in ["snp_id", "mod_site_id", "target_mod_code", "gene_names", "n_alt_reads", "n_modified", "p_value", "p_adj_bh", "effect_abs_delta_mod_frac"] if c in snp_mod_assoc_df.columns]
        snp_mod_html = df_to_html(snp_mod_assoc_df[keep].sort_values(["p_adj_bh", "effect_abs_delta_mod_frac"], ascending=[True, False]), max_rows=args.top_genes)

    joint_html = "<p class='muted'>No joint SNP-transcript-epitranscriptome dependency results available.</p>"
    if not snp_tx_mod_assoc_df.empty:
        keep = [c for c in ["snp_id", "mod_site_id", "target_mod_code", "n_transcripts_tested", "cmh_p_value", "cmh_p_adj_bh", "weighted_within_tx_effect", "classification"] if c in snp_tx_mod_assoc_df.columns]
        joint_html = df_to_html(snp_tx_mod_assoc_df[keep].sort_values(["cmh_p_adj_bh", "weighted_within_tx_effect"], ascending=[True, False]), max_rows=args.top_genes)

    hap_blocks_html = "<p class='muted'>No haplotype blocks available.</p>"
    if not hap_blocks_df.empty:
        keep = [c for c in ["block_id", "context_key", "chrom", "n_snps", "support_reads", "complete_reads", "haplotypes"] if c in hap_blocks_df.columns]
        hap_blocks_html = df_to_html(hap_blocks_df[keep].sort_values(["complete_reads", "n_snps"], ascending=False), max_rows=args.top_genes)

    hap_assoc_html = "<p class='muted'>No haplotype associations available.</p>"
    pieces = []
    if not hap_tx_assoc_df.empty:
        keep = [c for c in ["block_id", "n_transcripts_tested", "p_value", "p_adj_bh", "effect_max_abs_tx_frac_diff"] if c in hap_tx_assoc_df.columns]
        pieces.append("<h3>Haplotype to Transcript</h3>" + df_to_html(hap_tx_assoc_df[keep].sort_values(["p_adj_bh", "effect_max_abs_tx_frac_diff"], ascending=[True, False]), max_rows=args.top_genes))
    if not hap_mod_assoc_df.empty:
        keep = [c for c in ["block_id", "mod_site_id", "target_mod_code", "p_value", "p_adj_bh", "effect_max_abs_mod_rate_diff"] if c in hap_mod_assoc_df.columns]
        pieces.append("<h3>Haplotype to Epitranscriptome</h3>" + df_to_html(hap_mod_assoc_df[keep].sort_values(["p_adj_bh", "effect_max_abs_mod_rate_diff"], ascending=[True, False]), max_rows=args.top_genes))
    if pieces:
        hap_assoc_html = "".join(pieces)

    assignment_cols = [
        c for c in [
            "code", "gene_name", "gene_index", "transcript_index",
            "metagene_index", "zn_index", "read_support", "exact_chain_reads",
            "trunc_assigned_reads", "anchor_reads", "anchor_frac", "assignment_mode"
        ] if c in partition_map_df.columns
    ]
    sample_cols = [c for c in sample_stats_df.columns]
    read_cols = [c for c in read_stats_df.columns]
    length_cols = [c for c in tx_lengths_df.columns]

    cards_html = "".join(
        f"<div class='card'><div class='label'>{html.escape(label)}</div><div class='value'>{html.escape(value)}</div></div>"
        for label, value in overview_cards
    )

    body = []
    body.append(
        section(
            "Overview",
            f"<div class='cards'>{cards_html}</div><div class='hero'>{pca_html}</div>",
        )
    )
    body.append(section("Top Transcript Partitions", df_to_html(top_tx_df[top_tx_cols], max_rows=args.top_transcripts) if top_tx_cols else "<p class='muted'>No transcript summary columns available.</p>"))
    body.append(section("Sample Stats", df_to_html(sample_stats_df[sample_cols], max_rows=100) if sample_cols else "<p class='muted'>No sample stats available.</p>"))
    body.append(section("Read Funnel", df_to_html(read_stats_df[read_cols], max_rows=100) if read_cols else "<p class='muted'>No read funnel stats available.</p>"))
    body.append(section("Partition Map", df_to_html(partition_map_df[assignment_cols], max_rows=args.top_transcripts) if assignment_cols else "<p class='muted'>No partition map available.</p>"))
    body.append(section("Assigned Read Lengths", df_to_html(tx_lengths_df[length_cols], max_rows=args.top_transcripts) if length_cols else "<p class='muted'>No assigned read lengths available.</p>"))
    body.append(section("Top Gene / Mod Site Burden (ZN)", top_gene_sites_html))
    body.append(section("Overlap Resolution", overlap_html))
    body.append(section("Differential Sites", diff_html + diff_fig_html))
    body.append(section("Segregating SNP Candidates", candidate_snps_html))
    body.append(section("SNP to Transcript Associations", snp_tx_html))
    body.append(section("SNP to Epitranscriptome Associations", snp_mod_html))
    body.append(section("SNP Transcript Epitranscriptome Dependency", joint_html))
    body.append(section("Haplotype Blocks", hap_blocks_html))
    body.append(section("Haplotype Associations", hap_assoc_html))

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
    }}
    body {{
      margin: 0;
      font-family: Georgia, "Times New Roman", serif;
      background: linear-gradient(180deg, #efe5d8 0%, var(--bg) 100%);
      color: var(--ink);
    }}
    main {{
      max-width: 1200px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2 {{
      font-weight: 700;
      letter-spacing: 0.01em;
    }}
    section {{
      background: var(--card);
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 20px;
      margin-bottom: 20px;
      box-shadow: 0 10px 30px rgba(70, 44, 20, 0.08);
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
      gap: 12px;
      margin-bottom: 20px;
    }}
    .card {{
      background: var(--accent-soft);
      border-radius: 14px;
      padding: 14px;
    }}
    .label {{
      font-size: 0.85rem;
      color: #5c4d3d;
      margin-bottom: 6px;
    }}
    .value {{
      font-size: 1.5rem;
      color: var(--accent);
      font-weight: 700;
    }}
    img {{
      max-width: 100%;
      border-radius: 12px;
      border: 1px solid var(--line);
      background: white;
    }}
    .datatable {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.92rem;
    }}
    .datatable th, .datatable td {{
      border-bottom: 1px solid var(--line);
      padding: 8px 10px;
      text-align: left;
      vertical-align: top;
    }}
    .datatable th {{
      background: #f7eddc;
      position: sticky;
      top: 0;
    }}
    .muted {{
      color: #6d6255;
    }}
    .gallery {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
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
  </style>
</head>
<body>
  <main>
    <header>
      <h1>{html.escape(args.title)}</h1>
      <p class="muted">Generated from {html.escape(Path(args.out_html).name)} inputs in the current modulator run.</p>
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
