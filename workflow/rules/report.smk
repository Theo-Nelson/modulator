REPORT_CFG = config.get("report", {})
ENABLE_HTML_REPORT = _as_bool(REPORT_CFG.get("enable", True))

OUT_REPORT = f"results/report/{PREFIX}_report.html"


if ENABLE_HTML_REPORT:
    rule generate_html_report:
        input:
            classification = OUT_CLASS,
            metrics = OUT_METRICS,
            tx_counts = OUT_TX,
            pca_png = OUT_PCA,
            sample_stats = OUT_STATS,
            read_stats = OUT_READ_STATS,
            tx_lengths = OUT_TX_READ_LENGTHS,
            partition_map = OUT_PARTITION_MAP,
            extra = (
                ([f"results/aggregate_zn/{PREFIX}_FILTERED_sites_long.tsv"] if ENABLE_ZN and ENABLE_AGG_ZN else []) +
                ([f"results/aggregate_zt/{PREFIX}_FILTERED_long.tsv"] if ENABLE_ZT and ENABLE_AGG_ZT else []) +
                ([f"results/test_diffs/{PREFIX}__ZN_site_diff_results.tsv"] if ENABLE_ZN and ENABLE_AGG_ZN and ENABLE_TEST_DIFFS else []) +
                ([f"results/test_diffs/{PREFIX}__figs"] if ENABLE_ZN and ENABLE_AGG_ZN and ENABLE_TEST_DIFFS else []) +
                (MULTIGENE_SUMMARIES if ENABLE_MULTIGENE_FILTER else []) +
                ([
                    f"results/genotype/{PREFIX}_candidate_snps.tsv",
                    f"results/genotype/{PREFIX}_snp_transcript_assoc.tsv",
                    f"results/genotype/{PREFIX}_snp_mod_assoc.tsv",
                    f"results/genotype/{PREFIX}_snp_tx_mod_dependency.tsv",
                    f"results/genotype/{PREFIX}_haplotype_blocks.tsv",
                    f"results/genotype/{PREFIX}_haplotype_transcript_assoc.tsv",
                    f"results/genotype/{PREFIX}_haplotype_mod_assoc.tsv",
                ] if ENABLE_GENOTYPE else [])
            )
        output:
            OUT_REPORT
        params:
            script = "workflow/scripts/generate_html_report.py",
            zn_long = (f"results/aggregate_zn/{PREFIX}_FILTERED_sites_long.tsv" if ENABLE_ZN and ENABLE_AGG_ZN else ""),
            zt_long = (f"results/aggregate_zt/{PREFIX}_FILTERED_long.tsv" if ENABLE_ZT and ENABLE_AGG_ZT else ""),
            diff_tsv = (f"results/test_diffs/{PREFIX}__ZN_site_diff_results.tsv" if ENABLE_ZN and ENABLE_AGG_ZN and ENABLE_TEST_DIFFS else ""),
            diff_figs_dir = (f"results/test_diffs/{PREFIX}__figs" if ENABLE_ZN and ENABLE_AGG_ZN and ENABLE_TEST_DIFFS else ""),
            multigene_summary_glob = (f"{ZT_SCRAP_DIR}/*.multigene_filter_summary.tsv" if ENABLE_MULTIGENE_FILTER else ""),
            candidate_snps = (f"results/genotype/{PREFIX}_candidate_snps.tsv" if ENABLE_GENOTYPE else ""),
            snp_tx_assoc = (f"results/genotype/{PREFIX}_snp_transcript_assoc.tsv" if ENABLE_GENOTYPE else ""),
            snp_mod_assoc = (f"results/genotype/{PREFIX}_snp_mod_assoc.tsv" if ENABLE_GENOTYPE else ""),
            snp_tx_mod_assoc = (f"results/genotype/{PREFIX}_snp_tx_mod_dependency.tsv" if ENABLE_GENOTYPE else ""),
            hap_blocks = (f"results/genotype/{PREFIX}_haplotype_blocks.tsv" if ENABLE_GENOTYPE else ""),
            hap_tx_assoc = (f"results/genotype/{PREFIX}_haplotype_transcript_assoc.tsv" if ENABLE_GENOTYPE else ""),
            hap_mod_assoc = (f"results/genotype/{PREFIX}_haplotype_mod_assoc.tsv" if ENABLE_GENOTYPE else ""),
            title = REPORT_CFG.get("title", f"modulator report: {PREFIX}"),
            max_diff_figs = int(REPORT_CFG.get("max_diff_figs", 6)),
            top_transcripts = int(REPORT_CFG.get("top_transcripts", 20)),
            top_genes = int(REPORT_CFG.get("top_genes", 20))
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output})
            python {params.script} \
              --classification {input.classification} \
              --metrics {input.metrics} \
              --tx-counts {input.tx_counts} \
              --pca-png {input.pca_png} \
              --sample-stats {input.sample_stats} \
              --read-stats {input.read_stats} \
              --tx-lengths {input.tx_lengths} \
              --partition-map {input.partition_map} \
              --out-html {output} \
              --title "{params.title}" \
              --max-diff-figs {params.max_diff_figs} \
              --top-transcripts {params.top_transcripts} \
              --top-genes {params.top_genes} \
              --zn-long "{params.zn_long}" \
              --zt-long "{params.zt_long}" \
              --diff-results "{params.diff_tsv}" \
              --diff-figs-dir "{params.diff_figs_dir}" \
              --multigene-summary-glob "{params.multigene_summary_glob}" \
              --candidate-snps "{params.candidate_snps}" \
              --snp-tx-assoc "{params.snp_tx_assoc}" \
              --snp-mod-assoc "{params.snp_mod_assoc}" \
              --snp-tx-mod-assoc "{params.snp_tx_mod_assoc}" \
              --hap-blocks "{params.hap_blocks}" \
              --hap-tx-assoc "{params.hap_tx_assoc}" \
              --hap-mod-assoc "{params.hap_mod_assoc}"
            """
