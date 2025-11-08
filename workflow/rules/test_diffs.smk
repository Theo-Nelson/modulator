rule test_transcript_diffs:
    input:
        in_tsv = f"results/aggregate_zn/{config['prefix']}_FILTERED_sites_long.tsv"
    output:
        results_tsv = f"results/test_diffs/{config['prefix']}__ZN_site_diff_results.tsv",
        figs_dir = directory(f"results/test_diffs/{config['prefix']}__figs")
    params:
        out_prefix = f"results/test_diffs/{config['prefix']}",
        script = "scripts/test_stoichiometry_diffs.py",
        min_cov = config.get("min_cov_test", 20),
        topk = config.get("topk", 10)
    conda:
        "../envs/modulator.yaml"
    shell:
        """
        python {params.script} \
            --in-tsv {input.in_tsv} \
            --out-prefix {params.out_prefix} \
            --min-cov {params.min_cov} \
            --topk {params.topk} \
            --verbose
        """

