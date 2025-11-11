# convenience handle for nested assembler config
ASSEMBLER = config.get("assembler", {})

# precompute all CLI bits in params (no Python in the shell block)
def _tes_window_opt():
    v = ASSEMBLER.get("tes_window", None)
    return f"--tes-window {v}" if v not in (None, "null", "") else ""

def _flag(name, default=False, cli=""):
    val = ASSEMBLER.get(name, default)
    return cli if val else ""

rule assemble_transcripts:
    input:
        bamdir = lambda wc: BAMS_DIR,   # marks the directory as an input dep
        gtf    = REF_GTF
    output:
        gtf     = f"results/assemble/{config['prefix']}.gtf",
        summary = f"results/assemble/{config['prefix']}_classification_summary.tsv",
        counts  = f"results/assemble/{config['prefix']}_tx_counts.tsv",
        pcapng  = f"results/assemble/{config['prefix']}_tx_counts.pca.png",
        stats   = f"results/assemble/{config['prefix']}_per_sample_stats.tsv",
        tagged_bams = expand(f"results/assemble/zt_tagged/{{sample}}.zt_tagged.bam", sample=SAMPLES)
    params:
        out_gtf      = f"results/assemble/{config['prefix']}.gtf",
        out_prefix   = f"results/assemble/{config['prefix']}",
        script       = "workflow/scripts/assemble_transcripts.py",
        bam_glob     = BAM_GLOB,
        threads      = config.get("threads", 8),

        # numeric/string options
        min_mapq           = ASSEMBLER.get("min_mapq", 1),
        min_introns_read   = ASSEMBLER.get("min_introns_read", 0),
        require_softclip3p = ASSEMBLER.get("require_softclip3p", 0),
        apa_window         = ASSEMBLER.get("apa_window", 20),
        tes_window_opt     = _tes_window_opt(),
        min_reads          = ASSEMBLER.get("min_reads", 1),
        min_frac           = ASSEMBLER.get("min_frac", 0.0001),
        min_introns        = ASSEMBLER.get("min_introns", 0),
        min_polya_length   = ASSEMBLER.get("min_polya_length", 12),
        min_polya_purity   = ASSEMBLER.get("min_polya_purity", 0.7),
        polya_support_frac = ASSEMBLER.get("polya_support_frac", 0.6),
        tes_match_tol      = ASSEMBLER.get("tes_match_tol", 25),
        exact_tes_tol      = ASSEMBLER.get("exact_tes_tol", 10),
        min_reads_per_sample_for_mod = ASSEMBLER.get("min_reads_per_sample_for_mod", 5),
        min_total_reads_for_mod      = ASSEMBLER.get("min_total_reads_for_mod", 20),

        # boolean flags (emit only when true)
        primary_only_flag        = _flag("primary_only", False, "--primary-only"),
        write_zt_bams_flag       = _flag("write_zt_bams", False, "--write-zt-bams"),
        write_zt_tagged_flag     = _flag("write_zt_tagged_sample_bams", True, "--write-zt-tagged-sample-bams"),
        emit_modkit_manifest_flag= _flag("emit_modkit_manifest", False, "--emit-modkit-manifest")
    conda:
        "../envs/modulator.yaml"
    shell:
        r"""
        python {params.script} \
            --dir {input.bamdir} \
            --glob "{params.bam_glob}" \
            --gtf {input.gtf} \
            --out-gtf {params.out_gtf} \
            --out-prefix {params.out_prefix} \
            --threads {params.threads} \
            {params.primary_only_flag} \
            --min-mapq {params.min_mapq} \
            --min-introns-read {params.min_introns_read} \
            --require-softclip3p {params.require_softclip3p} \
            --apa-window {params.apa_window} \
            {params.tes_window_opt} \
            --min-reads {params.min_reads} \
            --min-frac {params.min_frac} \
            --min-introns {params.min_introns} \
            --min-polya-length {params.min_polya_length} \
            --min-polya-purity {params.min_polya_purity} \
            --polya-support-frac {params.polya_support_frac} \
            --tes-match-tol {params.tes_match_tol} \
            --exact-tes-tol {params.exact_tes_tol} \
            {params.write_zt_bams_flag} \
            {params.write_zt_tagged_flag} \
            {params.emit_modkit_manifest_flag} \
            --min-reads-per-sample-for-mod {params.min_reads_per_sample_for_mod} \
            --min-total-reads-for-mod {params.min_total_reads_for_mod}
        """

