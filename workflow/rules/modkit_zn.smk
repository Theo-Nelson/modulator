rule modkit_pileup_zn:
    input:
        bam = "results/assemble/zt_tagged/{sample}.zt_tagged.bam",
        ref = REF_FA
    output:
        beds = directory("results/modkit_zn/{sample}")
    params:
        mods = config["mods"],
        ref_bases = config["ref_bases"],
        common = MODKIT_CFG.get("common", {}),
        zn_cfg = MODKIT_CFG.get("zn", {})
    conda:
        "../envs/modulator.yaml"
    run:
        import os, shlex
        os.makedirs(output.beds, exist_ok=True)

        base_flags = _as_flags_from_common(params.common, sample=wildcards.sample, which="modkit_zn")

        # partition tag (ZN)
        ptag = params.zn_cfg.get("partition_tag", "ZN")
        base_flags_ptag = base_flags + ["--partition-tag", ptag]

        # One output bed per mod (your original behavior)
        for mod, ref_base in zip(params.mods, params.ref_bases):
            out_bed = os.path.join(output.beds, f"{mod}_filtered_mod.bed")

            cmd = (
                ["modkit", "pileup",
                 shlex.quote(input.bam),
                 shlex.quote(out_bed),
                 "--ref", shlex.quote(input.ref)]
                + base_flags_ptag
            )

            shell(" ".join(cmd))

