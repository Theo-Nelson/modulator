rule aggregate_modkit_zt:
    input:
        modkit_dirs = expand("results/modkit_zt/{sample}", sample=SAMPLES),
        summary_tsv = f"results/assemble/{config['prefix']}_classification_summary.tsv"
    output:
        raw_long    = f"results/aggregate_zt/{config['prefix']}_RAW_long.tsv",
        raw_frac    = f"results/aggregate_zt/{config['prefix']}_RAW_frac_pivot.tsv",
        raw_cov     = f"results/aggregate_zt/{config['prefix']}_RAW_cov_pivot.tsv",
        raw_nmod    = f"results/aggregate_zt/{config['prefix']}_RAW_Nmod_pivot.tsv",
        filt_long   = f"results/aggregate_zt/{config['prefix']}_FILTERED_long.tsv",
        filt_frac   = f"results/aggregate_zt/{config['prefix']}_FILTERED_frac_pivot.tsv",
        filt_cov    = f"results/aggregate_zt/{config['prefix']}_FILTERED_cov_pivot.tsv",
        filt_nmod   = f"results/aggregate_zt/{config['prefix']}_FILTERED_Nmod_pivot.tsv"
    params:
        out_prefix = f"results/aggregate_zt/{config['prefix']}",
        script     = "scripts/aggregate_by_transcript.py",
        min_cov    = config.get("min_cov", 5),
        enabled    = config.get("toggles", {}).get("enable_zt_aggregate", True),

        emit_raw_flag = (
            "--emit-raw"
            if config.get("aggregation", {}).get("zt", {}).get(
                "emit_raw", config.get("aggregate_outputs", {}).get("emit_raw", True)
            )
            else "--no-emit-raw"
        ),
        emit_filtered_flag = (
            "--emit-filtered"
            if config.get("aggregation", {}).get("zt", {}).get(
                "emit_filtered", config.get("aggregate_outputs", {}).get("emit_filtered", True)
            )
            else "--no-emit-filtered"
        ),
        write_long_flag = (
            "--write-long"
            if config.get("aggregation", {}).get("zt", {}).get(
                "write_long", config.get("aggregate_outputs", {}).get("write_long", True)
            )
            else "--no-write-long"
        ),
        write_pivots_flag = (
            "--write-pivots"
            if config.get("aggregation", {}).get("zt", {}).get(
                "write_pivots", config.get("aggregate_outputs", {}).get("write_pivots", True)
            )
            else "--no-write-pivots"
        ),

        # filtering knobs (same defaults you’re already using)
        filt_enable = config.get("aggregation", {}).get("zt", {}).get(
            "filter_enable", config.get("filters", {}).get("enable_site_filter", True)
        ),
        count_diff_factor = config.get("aggregation", {}).get("zt", {}).get(
            "count_diff_factor", config.get("filters", {}).get("count_diff_factor", 3)
        ),
        mod_fail_margin = config.get("aggregation", {}).get("zt", {}).get(
            "mod_fail_margin", config.get("filters", {}).get("mod_fail_margin", 1)
        ),
        filt_flag = (
            "--filter-enable"
            if config.get("aggregation", {}).get("zt", {}).get(
                "filter_enable", config.get("filters", {}).get("enable_site_filter", True)
            )
            else ""
        )
    conda:
        "../envs/modulator.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.raw_long})
        if [[ "{params.enabled}" =~ ^[Tt]rue$ ]]; then
            python {params.script} \
                --modkit-dir results/modkit_zt \
                --summary-tsv {input.summary_tsv} \
                --out-prefix {params.out_prefix} \
                --min-cov {params.min_cov} \
                {params.filt_flag} \
                --count-diff-factor {params.count_diff_factor} \
                --mod-fail-margin {params.mod_fail_margin} \
                {params.emit_raw_flag} \
                {params.emit_filtered_flag} \
                {params.write_long_flag} \
                {params.write_pivots_flag} \
                --debug-summary \
                --verbose
        else
            # disabled: touch all outputs (empty) so the DAG can finish
            touch {output.raw_long} {output.raw_frac} {output.raw_cov} {output.raw_nmod} \
                  {output.filt_long} {output.filt_frac} {output.filt_cov} {output.filt_nmod}
        fi
        """

