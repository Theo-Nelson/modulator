rule modkit_pileup_zn:
    input:
        bam = "results/assemble/zt_tagged/{sample}.zt_tagged.bam",
        ref = f"../{config['reference_genome']}"
    output:
        beds = directory("results/modkit_zn/{sample}")
    params:
        mods = config["mods"],
        ref_bases = config["ref_bases"],
        threads = config["threads"],
        modkit_bin = "../bin/modkit"
    conda:
        "../envs/modulator.yaml"
    run:
        import os
        os.makedirs(output.beds, exist_ok=True)
        for mod, ref_base in zip(params.mods, params.ref_bases):
            shell(
                f"{params.modkit_bin} pileup {input.bam} {output.beds}/{mod}_filtered_mod.bed "
                f"--ref {input.ref} "
                f"--filter-threshold {ref_base}:0.8 "
                f"--mod-thresholds {mod}:0.99 "
                f"--partition-tag ZN "
                f"--max-depth 1000 "
                f"--interval-size 100000 "
                f"-t {params.threads} "
                f"--log-filepath {output.beds}/{mod}_bed_log"
            )