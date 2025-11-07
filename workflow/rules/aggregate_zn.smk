rule aggregate_modkit_zn:
    input:
        modkit_dirs = expand("results/modkit_zn/{sample}", sample=config["samples"]),
        gtf = f"results/assemble/{config['prefix']}.gtf"
    output:
        long_tsv = f"results/aggregate_zn/{config['prefix']}_sites_long.tsv",
        pivots = directory(f"results/aggregate_zn/{config['prefix']}__per_gene_mod")
    params:
        out_prefix = f"results/aggregate_zn/{config['prefix']}",
        script = "scripts/aggregate_by_gene.py",
        min_cov = config.get("min_cov", 5)
    conda:
        "../envs/modulator.yaml"
    shell:
        """
        python {params.script} \
            --modkit-dir results/modkit_zn \
            --gtf {input.gtf} \
            --out-prefix {params.out_prefix} \
            --min-cov {params.min_cov} \
            --verbose
        """