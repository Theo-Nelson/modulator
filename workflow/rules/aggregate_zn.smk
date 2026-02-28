rule aggregate_modkit_zn:
    input:
        modkit_dirs = expand("results/modkit_zn/{sample}", sample=SAMPLES),
        gtf = f"results/assemble/{config['prefix']}.gtf"
    output:
        long_tsv = f"results/aggregate_zn/{config['prefix']}_FILTERED_sites_long.tsv",
        pivots   = directory(f"results/aggregate_zn/{config['prefix']}_FILTERED__per_gene_mod")
    params:
        out_prefix = f"results/aggregate_zn/{config['prefix']}",
        script     = "workflow/scripts/aggregate_by_gene.py",
        min_cov    = config.get("min_cov", 5),
        enabled    = config.get("toggles", {}).get("enable_zn_aggregate", True),

        # NEW: external-sort tuning (optional but recommended)
        tmpdir = config.get("aggregation", {}).get("tmpdir", os.environ.get("TMPDIR", "/tmp")),
        chunk_lines = config.get("aggregation", {}).get("chunk_lines", 2000000),

        emit_raw_flag = (
            "--emit-raw"
            if config.get("aggregation", {}).get("zn", {}).get(
                "emit_raw", config.get("aggregate_outputs", {}).get("emit_raw", True)
            )
            else "--no-emit-raw"
        ),
        emit_filtered_flag = (
            "--emit-filtered"
            if config.get("aggregation", {}).get("zn", {}).get(
                "emit_filtered", config.get("aggregate_outputs", {}).get("emit_filtered", True)
            )
            else "--no-emit-filtered"
        ),
        write_long_flag = (
            "--write-long"
            if config.get("aggregation", {}).get("zn", {}).get(
                "write_long", config.get("aggregate_outputs", {}).get("write_long", True)
            )
            else "--no-write-long"
        ),
        write_pivots_flag = (
            "--write-pivots"
            if config.get("aggregation", {}).get("zn", {}).get(
                "write_pivots", config.get("aggregate_outputs", {}).get("write_pivots", True)
            )
            else "--no-write-pivots"
        ),
        write_raw_per_gene_flag = (
            "--write-raw-per-gene"
            if config.get("aggregation", {}).get("zn", {}).get(
                "write_raw_per_gene",
                config.get("aggregate_outputs", {}).get("write_raw_per_gene_tables",
                    config.get("aggregate_outputs", {}).get("write_raw_per_gene", False)
                )
            )
            else "--no-write-raw-per-gene"
        ),
        write_filtered_per_gene_flag = (
            "--write-filtered-per-gene"
            if config.get("aggregation", {}).get("zn", {}).get(
                "write_filtered_per_gene",
                config.get("aggregate_outputs", {}).get("write_filtered_per_gene_tables",
                    config.get("aggregate_outputs", {}).get("write_filtered_per_gene", True)
                )
            )
            else "--no-write-filtered-per-gene"
        ),

        filt_enable = config.get("aggregation", {}).get("zn", {}).get(
            "filter_enable", config.get("filters", {}).get("enable_site_filter", True)
        ),
        count_diff_factor = config.get("aggregation", {}).get("zn", {}).get(
            "count_diff_factor", config.get("filters", {}).get("count_diff_factor", 3)
        ),
        mod_fail_margin = config.get("aggregation", {}).get("zn", {}).get(
            "mod_fail_margin", config.get("filters", {}).get("mod_fail_margin", 1)
        ),
        filt_flag = "--filter-enable" if config.get("aggregation", {}).get("zn", {}).get(
            "filter_enable", config.get("filters", {}).get("enable_site_filter", True)
        ) else ""
    conda:
        "../envs/modulator.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.long_tsv}) {output.pivots}
        if [[ "{params.enabled}" =~ ^[Tt]rue$ ]]; then
            python {params.script} \
            --modkit-dir results/modkit_zn \
            --gtf {input.gtf} \
            --out-prefix {params.out_prefix} \
            --min-cov {params.min_cov} \
            --tmpdir {params.tmpdir} \
            --chunk-lines {params.chunk_lines} \
            {params.filt_flag} \
            --count-diff-factor {params.count_diff_factor} \
            --mod-fail-margin {params.mod_fail_margin} \
            {params.emit_raw_flag} \
            {params.emit_filtered_flag} \
            {params.write_long_flag} \
            {params.write_pivots_flag} \
            {params.write_raw_per_gene_flag} \
            {params.write_filtered_per_gene_flag} \
            --verbose
        else
            echo -e "sample\tZN_transcript_index\tchrom\tstart0\tend0\tstrand\tmod_code\tNvalid_cov\tNmod\tfrac_modified\tgene_id\tgene_name\tNcanonical\tNother_mod\tNdelete\tNfail\tNdiff\tNnocall" > {output.long_tsv}
            : > {output.pivots}/.placeholder
        fi
        """
