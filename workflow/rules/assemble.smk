rule assemble_transcripts:
    input:
        bams = lambda wildcards: [get_bam_path(s) for s in config["samples"]],
        gtf = REF_GTF
    output:
        gtf = f"results/assemble/{config['prefix']}.gtf",
        summary = f"results/assemble/{config['prefix']}_classification_summary.tsv",
        tagged_bams = expand(f"results/assemble/zt_tagged/{{sample}}.zt_tagged.bam", sample=config["samples"])
    params:
        out_gtf = f"results/assemble/{config['prefix']}.gtf",
        script = "scripts/assemble_transcripts.py",
        min_reads = config.get("min_reads", 10),
        min_frac = config.get("min_frac", 0.05)
    conda:
        "../envs/modulator.yaml"
    shell:
        """
        python {params.script} \
            --bams {input.bams} \
            --gtf {input.gtf} \
            --out-gtf {params.out_gtf} \
            --min-reads {params.min_reads} \
            --min-frac {params.min_frac} \
            --write-zt-tagged-sample-bams
        """
