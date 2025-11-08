rule modkit_pileup_zt:
    input:
        bam = "results/assemble/zt_tagged/{sample}.zt_tagged.bam",
        ref = REF_FA
    output:
        bed_dir = directory("results/modkit_zt/{sample}")
    params:
        common = MODKIT_CFG.get("common", {}),
        zt_cfg = MODKIT_CFG.get("zt", {})
    conda:
        "../envs/modulator.yaml"
    run:
        import os, shlex
        os.makedirs(output.bed_dir, exist_ok=True)

        flags = _as_flags_from_common(params.common, sample=wildcards.sample, which="modkit_zt")

        # partition tag (ZT)
        ptag = params.zt_cfg.get("partition_tag", "ZT")
        flags = flags + ["--partition-tag", ptag]

        cmd = (
            ["modkit", "pileup",
             shlex.quote(input.bam),
             shlex.quote(output.bed_dir),
             "--ref", shlex.quote(input.ref)]
            + flags
        )

        shell(" ".join(cmd))

