rule test_transcript_diffs:
    input:
        in_tsv = f"results/aggregate_zn/{config['prefix']}_FILTERED_sites_long.tsv"
    output:
        results_tsv = f"results/test_diffs/{config['prefix']}__ZN_site_diff_results.tsv",
        figs_dir = directory(f"results/test_diffs/{config['prefix']}__figs")
    params:
        out_prefix = f"results/test_diffs/{config['prefix']}",
        script = "workflow/scripts/test_stoichiometry_diffs.py",
        # existing knobs
        min_cov = config.get("min_cov_test", 20),
        topk = config.get("topk", 10),
        # NEW: read test_diffs block safely
        test_diffs = config.get("test_diffs", {}),
        test_flag = (lambda td: f"--test {td['test']}" if isinstance(td.get("test"), str) and td["test"] else "")(config.get("test_diffs", {})),
        alt_flag  = (lambda td: f"--alternative {td['alternative']}" if isinstance(td.get("alternative"), str) and td["alternative"] else "")(config.get("test_diffs", {})),
        pc_flag   = (lambda td: f"--pseudocount {td['pseudocount']}" if td.get("pseudocount") is not None else "")(config.get("test_diffs", {})),
        gene_flags = (lambda td: " ".join(f"--gene-filter '{g}'" for g in td["gene_filter"])
                                 if isinstance(td.get("gene_filter"), (list, tuple)) and td["gene_filter"] else "")(config.get("test_diffs", {})),
        mod_flags  = (lambda td: " ".join(f"--mod-filter '{m}'" for m in td["mod_filter"])
                                 if isinstance(td.get("mod_filter"), (list, tuple)) and td["mod_filter"] else "")(config.get("test_diffs", {})),
    conda:
        "../envs/modulator.yaml"
    shell:
        r"""
        python {params.script} \
            --in-tsv {input.in_tsv} \
            --out-prefix {params.out_prefix} \
            --min-cov {params.min_cov} \
            --topk {params.topk} \
            {params.test_flag} \
            {params.pc_flag} \
            {params.alt_flag} \
            {params.gene_flags} \
            {params.mod_flags} \
            --verbose
        """

