rule aggregate_modkit_zt:
    input:
        modkit_dirs = expand("results/modkit_zt/{sample}", sample=config["samples"]),
        summary_tsv = f"results/assemble/{config['prefix']}_classification_summary.tsv"
    output:
        long_tsv = f"results/aggregate_zt/{config['prefix']}_long.tsv",
        frac_pivot = f"results/aggregate_zt/{config['prefix']}_frac_pivot.tsv",
        cov_pivot = f"results/aggregate_zt/{config['prefix']}_cov_pivot.tsv",
        nmod_pivot = f"results/aggregate_zt/{config['prefix']}_Nmod_pivot.tsv"
    params:
        out_prefix = f"results/aggregate_zt/{config['prefix']}",
        script = "scripts/aggregate_by_transcript.py",
        min_cov = config.get("min_cov", 0)
    conda:
        "../envs/modulator.yaml"
    shell:
        """
        python {params.script} \
            --modkit-dir results/modkit_zt \
            --summary-tsv {input.summary_tsv} \
            --out-prefix {params.out_prefix} \
            --min-cov {params.min_cov} \
            --verbose
        """