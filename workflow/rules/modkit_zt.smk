rule modkit_pileup_zt:
    input:
        bam = "results/assemble/zt_tagged/{sample}.zt_tagged.bam",
        ref = REF_FA
    output:
        bed_dir = directory("results/modkit_zt/{sample}")
    params:
        common  = MODKIT_CFG.get("common", {}),
        zt_cfg  = MODKIT_CFG.get("zt", {}),
        enabled = config.get("toggles", {}).get("enable_zt_pileup", True)
    conda:
        "../envs/modulator.yaml"
    run:
        import os, shlex

        # ensure output directory exists
        os.makedirs(output.bed_dir, exist_ok=True)

        # allow toggle to disable this step cleanly
        if not _as_bool(params.enabled, True):
            with open(os.path.join(output.bed_dir, ".SKIPPED"), "w") as fh:
                fh.write("ZT pileup disabled\n")
            return

        # render null-safe flags from common config
        flags = _as_flags_from_common(params.common, sample=wildcards.sample, which="modkit_zt")

        # add partition tag if provided
        ptag = params.zt_cfg.get("partition_tag", "ZT")
        if _is_set(ptag):
            flags += ["--partition-tag", str(ptag)]

        # build and run the command (quote paths; leave flags as tokens)
        cmd = (
            ["modkit", "pileup",
             shlex.quote(input.bam),
             shlex.quote(output.bed_dir),
             "--ref", shlex.quote(input.ref)]
            + flags
        )
        shell(" ".join(cmd))

