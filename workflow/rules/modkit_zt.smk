rule modkit_pileup_zt:
    input:
        bam = "results/assemble/zt_tagged/{sample}.zt_tagged.bam",
        ref = REF_FA
    output:
        bed_dir = directory("results/modkit_zt/{sample}")
    params:
        mods = config["mods"],
        ref_bases = config["ref_bases"],
        threads = config["threads"],
        modkit_bin = "../bin/modkit"
    conda:
        "../envs/modulator.yaml"
    shell:
        """
        {params.modkit_bin} pileup {input.bam} {output.bed_dir} \
            --ref {input.ref} \
            --filter-threshold A:0.8 \
            --filter-threshold C:0.8 \
            --filter-threshold G:0.8 \
            --filter-threshold T:0.8 \
            --mod-thresholds 17596:0.99 \
            --mod-thresholds a:0.99 \
            --mod-thresholds m:0.99 \
            --mod-thresholds 17802:0.99 \
            --mod-thresholds 69426:0.99 \
            --mod-thresholds 19228:0.99 \
            --mod-thresholds 19229:0.99 \
            --mod-thresholds 19227:0.99 \
            --log-filepath results/modkit_zt/{wildcards.sample}.log \
            --max-depth 1000 \
            --interval-size 100000 \
            --prefix {wildcards.sample} \
            --partition-tag ZT \
            -t {params.threads}
        """
