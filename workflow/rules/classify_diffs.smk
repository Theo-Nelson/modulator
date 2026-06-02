rule classify_diff_sites:
    input:
        diff_tsv = f"results/test_diffs/{config['prefix']}__ZN_site_diff_results.tsv",
        gtf      = f"results/assemble/{config['prefix']}.gtf",
        zn_long  = f"results/aggregate_zn/{config['prefix']}_FILTERED_sites_long.tsv"
    output:
        classified_tsv = f"results/test_diffs/{config['prefix']}__ZN_site_classified.tsv"
    params:
        script = "workflow/scripts/classify_diff_sites.py",
        cd = config.get("classify_diffs", {}),
        fig_flags = (
                        (
                            f"--arch-figs-dir results/test_diffs/{config['prefix']}__figs_by_category_arch "
                            f"--zn-long results/aggregate_zn/{config['prefix']}_FILTERED_sites_long.tsv "
                            f"--figs-dir results/test_diffs/{config['prefix']}__figs_by_category "
                            f"--figs-per-category {int(config.get('classify_diffs', {}).get('figs_per_category', 10))}"
                        )
                        if _as_bool(config.get("classify_diffs", {}).get("figures", True)) else ""
                     ),
        min_effect = config.get("classify_diffs", {}).get("min_effect", 0.10),
        fdr        = config.get("classify_diffs", {}).get("fdr", 0.05),
        min_cov    = config.get("classify_diffs", {}).get("min_cov", 0),
        tes_tol    = config.get("classify_diffs", {}).get("tes_tol", 200),
        inside_tol = config.get("classify_diffs", {}).get("inside_tol", 50),
        ejc_nt     = config.get("classify_diffs", {}).get("ejc_nt", 150),
        intergenic_gap = config.get("classify_diffs", {}).get("intergenic_gap", 1000),
        # mod_filter precedence: classify_diffs.mod_filter, then test_diffs.mod_filter.
        # When neither is set, emit NO --mod-filter so ALL detected modifications are
        # classified (the diff table already carries every mod_code emitted upstream).
        mod_flags  = (lambda cd, td: " ".join(
                        f"--mod-filter '{m}'" for m in (
                            cd["mod_filter"] if isinstance(cd.get("mod_filter"), (list, tuple)) and cd.get("mod_filter")
                            else (td["mod_filter"] if isinstance(td.get("mod_filter"), (list, tuple)) and td.get("mod_filter") else [])
                        )
                     ))(config.get("classify_diffs", {}), config.get("test_diffs", {})),
    conda:
        "../envs/modulator.yaml"
    shell:
        r"""
        python {params.script} \
            --diff-tsv {input.diff_tsv} \
            --gtf {input.gtf} \
            --out-tsv {output.classified_tsv} \
            --min-effect {params.min_effect} \
            --fdr {params.fdr} \
            --min-cov {params.min_cov} \
            --tes-tol {params.tes_tol} \
            --inside-tol {params.inside_tol} \
            --ejc-nt {params.ejc_nt} \
            --intergenic-gap {params.intergenic_gap} \
            {params.mod_flags} \
            {params.fig_flags} \
            --verbose
        """
