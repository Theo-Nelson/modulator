rule modkit_pileup_zn:
    input:
        bam = "results/assemble/zt_tagged/{sample}.zt_tagged.bam",
        ref = REF_FA
    output:
        bed_dir = directory("results/modkit_zn/{sample}")
    params:
        common = MODKIT_CFG.get("common", {}),
        zn_cfg = MODKIT_CFG.get("zn", {}),
        enabled = config.get("toggles", {}).get("enable_zn_pileup", True)
    conda:
        "../envs/modulator.yaml"
    run:
        import os, shlex
        os.makedirs(output.bed_dir, exist_ok=True)

        if not params.enabled:
            open(os.path.join(output.bed_dir, ".SKIPPED"), "w").write("ZN pileup disabled\n")
            return

        flags = _as_flags_from_common(params.common, sample=wildcards.sample, which="modkit_zn")
        ptag = params.zn_cfg.get("partition_tag", "ZN")
        flags = flags + ["--partition-tag", ptag]

        cmd = (
            ["modkit", "pileup",
             shlex.quote(input.bam),
             shlex.quote(output.bed_dir),
             "--ref", shlex.quote(input.ref)]
            + flags
        )
        shell(" ".join(cmd))

