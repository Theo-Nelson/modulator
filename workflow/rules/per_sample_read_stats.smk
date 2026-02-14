# workflow/rules/per_sample_read_stats.smk

OUT_READ_STATS = f"{ASSEMBLE_DIR}/{PREFIX}_per_sample_read_stats.tsv"

rule per_sample_read_stats:
    input:
        # ensure assembly products exist
        gtf     = OUT_GTF,
        stats   = OUT_STATS,
        zt_dir  = directory(ZT_DIR),
        zt_bams = ZT_TAGGED_BAMS
    output:
        OUT_READ_STATS
    threads: 1
    conda:
        "../envs/modulator.yaml"
    params:
        bams_dir = _BAMS_DIR,
        bam_glob = _BAM_GLOB,

        primary_only_flag  = "--primary-only" if config.get("assembler", {}).get("primary_only", True) else "",
        min_mapq           = config.get("assembler", {}).get("min_mapq", 10),
        min_introns_read   = config.get("assembler", {}).get("min_introns_read", 1),
        require_softclip3p = config.get("assembler", {}).get("require_softclip3p", 0),
    shell:
        r"""
        set -euo pipefail
        SCRIPT_A="{workflow.basedir}/scripts/per_sample_read_stats.py"
        SCRIPT_B="{workflow.basedir}/workflow/scripts/per_sample_read_stats.py"
        if [ -f "$SCRIPT_A" ]; then
            SCRIPT="$SCRIPT_A"
        elif [ -f "$SCRIPT_B" ]; then
            SCRIPT="$SCRIPT_B"
        else
            echo "ERROR: per_sample_read_stats.py not found at $SCRIPT_A or $SCRIPT_B" >&2
            exit 2
        fi

        python "$SCRIPT" \
          --bams-dir "{params.bams_dir}" \
          --bam-glob "{params.bam_glob}" \
          --zt-tagged-dir "{input.zt_dir}" \
          --out "{output}" \
          {params.primary_only_flag} \
          --min-mapq "{params.min_mapq}" \
          --min-introns-read "{params.min_introns_read}" \
          --require-softclip3p "{params.require_softclip3p}"
        """
