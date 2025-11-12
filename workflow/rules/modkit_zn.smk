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
    threads: 8
    conda:
        "../envs/modulator.yaml"
    run:
        import os, shlex

        os.makedirs(output.bed_dir, exist_ok=True)

        if not _as_bool(params.enabled, True):
            with open(os.path.join(output.bed_dir, ".SKIPPED"), "w") as fh:
                fh.write("ZN pileup disabled\n")
            return

        # flags from config, but ensure no thread flags sneak in
        flags = _as_flags_from_common(params.common, sample=wildcards.sample, which="modkit_zn")

        def _strip_thread_flags(fl):
            out, skip = [], False
            for tok in fl:
                if skip:
                    skip = False
                    continue
                if tok in ("-t", "--threads"):
                    skip = True
                    continue
                out.append(tok)
            return out

        flags = _strip_thread_flags(flags)

        # partition tag
        ptag = params.zn_cfg.get("partition_tag", "ZN")
        if _is_set(ptag):
            flags += ["--partition-tag", str(ptag)]

        # enforce rule threads
        flags = ["-t", str(threads)] + flags

        cmd = (["modkit", "pileup",
                shlex.quote(input.bam),
                shlex.quote(output.bed_dir),
                "--ref", shlex.quote(input.ref)] + flags)
        shell(" ".join(cmd))

