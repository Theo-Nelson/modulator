MULTIGENE_FILTER_CFG = config.get("multigene_filter", {})
ENABLE_MULTIGENE_FILTER = _as_bool(MULTIGENE_FILTER_CFG.get("enable", True))

ZT_FILTERED_DIR = f"{ASSEMBLE_DIR}/zt_filtered"
ZT_SCRAP_DIR = f"{ASSEMBLE_DIR}/zt_scrap"

ZT_FILTERED_BAM_PATTERN = f"{ZT_FILTERED_DIR}" + "/{sample}.zt_tagged.clean.bam"
ZT_SCRAP_BAM_PATTERN = f"{ZT_SCRAP_DIR}" + "/{sample}.zt_tagged.multigene_scrap.bam"

ZT_FILTERED_BAMS = expand(ZT_FILTERED_BAM_PATTERN, sample=SAMPLES)
ZT_SCRAP_BAMS = expand(ZT_SCRAP_BAM_PATTERN, sample=SAMPLES)
MULTIGENE_SUMMARIES = expand(f"{ZT_SCRAP_DIR}" + "/{sample}.multigene_filter_summary.tsv", sample=SAMPLES)
MULTIGENE_REMOVED_READS = expand(f"{ZT_SCRAP_DIR}" + "/{sample}.multigene_removed_reads.tsv", sample=SAMPLES)
MULTIGENE_SCRAP_TX_COUNT_FILES = expand(f"{ZT_SCRAP_DIR}" + "/{sample}.multigene_scrap_tx_counts.tsv", sample=SAMPLES)


if ENABLE_MULTIGENE_FILTER:
    rule filter_multigene_reads_from_zt_bam:
        input:
            bam = "results/assemble/zt_tagged/{sample}.zt_tagged.bam",
            gtf = OUT_GTF
        output:
            clean_bam = ZT_FILTERED_BAM_PATTERN,
            scrap_bam = ZT_SCRAP_BAM_PATTERN,
            summary = f"{ZT_SCRAP_DIR}" + "/{sample}.multigene_filter_summary.tsv",
            removed = f"{ZT_SCRAP_DIR}" + "/{sample}.multigene_removed_reads.tsv",
            scrap_tx_counts = f"{ZT_SCRAP_DIR}" + "/{sample}.multigene_scrap_tx_counts.tsv"
        params:
            zero_gene_action = MULTIGENE_FILTER_CFG.get("zero_gene_action", "keep")
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            set -euo pipefail

            SCRIPT_A="{workflow.basedir}/scripts/filter_multigene_reads_from_zt_bam.py"
            SCRIPT_B="{workflow.basedir}/workflow/scripts/filter_multigene_reads_from_zt_bam.py"
            if [ -f "$SCRIPT_A" ]; then
                SCRIPT="$SCRIPT_A"
            elif [ -f "$SCRIPT_B" ]; then
                SCRIPT="$SCRIPT_B"
            else
                echo "ERROR: filter_multigene_reads_from_zt_bam.py not found at $SCRIPT_A or $SCRIPT_B" >&2
                exit 2
            fi

            python "$SCRIPT" \
              --bam "{input.bam}" \
              --gtf "{input.gtf}" \
              --sample "{wildcards.sample}" \
              --out-clean-bam "{output.clean_bam}" \
              --out-scrap-bam "{output.scrap_bam}" \
              --out-summary-tsv "{output.summary}" \
              --out-removed-tsv "{output.removed}" \
              --out-scrap-tx-counts-tsv "{output.scrap_tx_counts}" \
              --zero-gene-action "{params.zero_gene_action}"
            """


    rule aggregate_multigene_scrap_tx_counts:
        input:
            counts = MULTIGENE_SCRAP_TX_COUNT_FILES
        output:
            f"{ASSEMBLE_DIR}/{PREFIX}_multigene_scrap_tx_counts.tsv"
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            set -euo pipefail

            SCRIPT_A="{workflow.basedir}/scripts/aggregate_scrap_tx_counts.py"
            SCRIPT_B="{workflow.basedir}/workflow/scripts/aggregate_scrap_tx_counts.py"
            if [ -f "$SCRIPT_A" ]; then
                SCRIPT="$SCRIPT_A"
            elif [ -f "$SCRIPT_B" ]; then
                SCRIPT="$SCRIPT_B"
            else
                echo "ERROR: aggregate_scrap_tx_counts.py not found at $SCRIPT_A or $SCRIPT_B" >&2
                exit 2
            fi

            python "$SCRIPT" \
              --counts {input.counts} \
              --out "{output}"
            """
