GENO_CFG = config.get("genotype", {})
ENABLE_GENOTYPE = _as_bool(GENO_CFG.get("enable", False))

GENO_DIR = "results/genotype"
GENO_INPUT_BAMS = expand(MODKIT_INPUT_BAM, sample=SAMPLES)

GENO_READ_ASSIGN = f"{GENO_DIR}/{PREFIX}_read_assignments.tsv"
GENO_SNPS = f"{GENO_DIR}/{PREFIX}_candidate_snps.tsv"
GENO_MOL_SNPS = f"{GENO_DIR}/{PREFIX}_molecule_snps.tsv"
GENO_MOD_SITES = f"{GENO_DIR}/{PREFIX}_candidate_mod_sites.tsv"
GENO_MOD_BED = f"{GENO_DIR}/{PREFIX}_candidate_mod_sites.bed"
GENO_MOL_MODS = f"{GENO_DIR}/{PREFIX}_molecule_mod_calls.tsv"
GENO_SNP_TX = f"{GENO_DIR}/{PREFIX}_snp_transcript_assoc.tsv"
GENO_SNP_MOD = f"{GENO_DIR}/{PREFIX}_snp_mod_assoc.tsv"
GENO_JOINT = f"{GENO_DIR}/{PREFIX}_snp_tx_mod_dependency.tsv"
GENO_HAP_BLOCKS = f"{GENO_DIR}/{PREFIX}_haplotype_blocks.tsv"
GENO_HAP_MOLS = f"{GENO_DIR}/{PREFIX}_molecule_haplotypes.tsv"
GENO_HAP_TX = f"{GENO_DIR}/{PREFIX}_haplotype_transcript_assoc.tsv"
GENO_HAP_MOD = f"{GENO_DIR}/{PREFIX}_haplotype_mod_assoc.tsv"

GENO_MOD_SOURCE_ZN = (f"results/aggregate_zn/{PREFIX}_FILTERED_sites_long.tsv" if ENABLE_ZN and ENABLE_AGG_ZN else "")
GENO_MOD_SOURCE_ZT = (f"results/aggregate_zt/{PREFIX}_FILTERED_long.tsv" if ENABLE_ZT and ENABLE_AGG_ZT else "")


if ENABLE_GENOTYPE:
    rule build_read_assignments:
        input:
            bams = GENO_INPUT_BAMS,
            summary = OUT_CLASS
        output:
            GENO_READ_ASSIGN
        params:
            script = "workflow/scripts/build_read_assignment_table.py"
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output})
            python {params.script} \
              --bams {input.bams} \
              --summary-tsv {input.summary} \
              --out-tsv {output} \
              --primary-only
            """

    rule discover_candidate_snps:
        input:
            bams = GENO_INPUT_BAMS,
            gtf = OUT_GTF,
            ref = REF_FA
        output:
            GENO_SNPS
        params:
            script = "workflow/scripts/discover_candidate_snps.py",
            min_alt_reads = int(GENO_CFG.get("min_alt_reads", 4)),
            min_total_cov = int(GENO_CFG.get("min_total_cov", 8)),
            min_alt_frac = float(GENO_CFG.get("min_alt_frac", 0.10)),
            max_alt_frac = float(GENO_CFG.get("max_alt_frac", 0.90)),
            min_baseq = int(GENO_CFG.get("min_baseq", 20)),
            min_mapq = int(GENO_CFG.get("min_mapq", config.get("assembler", {}).get("min_mapq", 10))),
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output})
            python {params.script} \
              --bams {input.bams} \
              --reference-fa {input.ref} \
              --gtf {input.gtf} \
              --out-tsv {output} \
              --min-alt-reads {params.min_alt_reads} \
              --min-total-cov {params.min_total_cov} \
              --min-alt-frac {params.min_alt_frac} \
              --max-alt-frac {params.max_alt_frac} \
              --min-baseq {params.min_baseq} \
              --min-mapq {params.min_mapq} \
              --primary-only
            """

    rule build_molecule_snp_table:
        input:
            bams = GENO_INPUT_BAMS,
            snps = GENO_SNPS
        output:
            GENO_MOL_SNPS
        params:
            script = "workflow/scripts/build_molecule_snp_table.py",
            min_baseq = int(GENO_CFG.get("min_baseq", 20)),
            min_mapq = int(GENO_CFG.get("min_mapq", config.get("assembler", {}).get("min_mapq", 10))),
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output})
            python {params.script} \
              --bams {input.bams} \
              --candidate-snps {input.snps} \
              --out-tsv {output} \
              --min-baseq {params.min_baseq} \
              --min-mapq {params.min_mapq} \
              --primary-only
            """

    rule build_candidate_mod_sites:
        input:
            extra = (
                ([GENO_MOD_SOURCE_ZN] if GENO_MOD_SOURCE_ZN else []) +
                ([GENO_MOD_SOURCE_ZT] if GENO_MOD_SOURCE_ZT else [])
            )
        output:
            tsv = GENO_MOD_SITES,
            bed = GENO_MOD_BED
        params:
            script = "workflow/scripts/build_candidate_mod_sites.py",
            zn_long = GENO_MOD_SOURCE_ZN,
            zt_long = GENO_MOD_SOURCE_ZT,
            min_total_cov = int(GENO_CFG.get("min_mod_site_cov", 1)),
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output.tsv})
            python {params.script} \
              --zn-long "{params.zn_long}" \
              --zt-long "{params.zt_long}" \
              --out-tsv {output.tsv} \
              --out-bed {output.bed} \
              --min-total-cov {params.min_total_cov}
            """

    rule build_molecule_mod_table:
        input:
            bams = GENO_INPUT_BAMS,
            sites = GENO_MOD_SITES,
            bed = GENO_MOD_BED,
            assignments = GENO_READ_ASSIGN,
            ref = REF_FA
        output:
            GENO_MOL_MODS
        params:
            script = "workflow/scripts/build_molecule_mod_table.py",
            threads = int(min(8, config.get("threads", 4))),
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output})
            python {params.script} \
              --bams {input.bams} \
              --candidate-sites-tsv {input.sites} \
              --candidate-bed {input.bed} \
              --read-assignments {input.assignments} \
              --reference-fa {input.ref} \
              --out-tsv {output} \
              --threads {params.threads}
            """

    rule test_snp_transcript_assoc:
        input:
            GENO_MOL_SNPS
        output:
            GENO_SNP_TX
        params:
            script = "workflow/scripts/test_snp_transcript_assoc.py",
            min_group = int(GENO_CFG.get("min_group_reads", 4)),
            test = GENO_CFG.get("test", "auto"),
            pseudocount = float(GENO_CFG.get("pseudocount", 0.5)),
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output})
            python {params.script} \
              --molecule-snps {input} \
              --out-tsv {output} \
              --min-allele-reads {params.min_group} \
              --min-transcript-reads {params.min_group} \
              --test {params.test} \
              --pseudocount {params.pseudocount}
            """

    rule test_snp_mod_assoc:
        input:
            snps = GENO_MOL_SNPS,
            mods = GENO_MOL_MODS
        output:
            GENO_SNP_MOD
        params:
            script = "workflow/scripts/test_snp_mod_assoc.py",
            min_group = int(GENO_CFG.get("min_group_reads", 4)),
            test = GENO_CFG.get("test", "auto"),
            pseudocount = float(GENO_CFG.get("pseudocount", 0.5)),
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output})
            python {params.script} \
              --molecule-snps {input.snps} \
              --molecule-mods {input.mods} \
              --out-tsv {output} \
              --min-allele-reads {params.min_group} \
              --min-total-reads {params.min_group} \
              --test {params.test} \
              --pseudocount {params.pseudocount}
            """

    rule test_snp_tx_mod_dependency:
        input:
            snps = GENO_MOL_SNPS,
            mods = GENO_MOL_MODS,
            snp_tx = GENO_SNP_TX,
            snp_mod = GENO_SNP_MOD
        output:
            GENO_JOINT
        params:
            script = "workflow/scripts/test_snp_tx_mod_dependency.py",
            min_group = int(GENO_CFG.get("min_group_reads", 4)),
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output})
            python {params.script} \
              --molecule-snps {input.snps} \
              --molecule-mods {input.mods} \
              --snp-transcript-assoc {input.snp_tx} \
              --snp-mod-assoc {input.snp_mod} \
              --out-tsv {output} \
              --min-stratum-reads {params.min_group}
            """

    rule build_haplotype_blocks:
        input:
            GENO_MOL_SNPS
        output:
            blocks = GENO_HAP_BLOCKS,
            mols = GENO_HAP_MOLS
        params:
            script = "workflow/scripts/build_haplotype_blocks.py",
            min_alt_reads = int(GENO_CFG.get("min_alt_reads", 4)),
            min_cocover = int(GENO_CFG.get("min_haplotype_reads", 4)),
            max_snps = int(GENO_CFG.get("max_haplotype_snps", 4)),
            min_hap_reads = int(GENO_CFG.get("min_haplotype_reads", 4)),
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output.blocks})
            python {params.script} \
              --molecule-snps {input} \
              --out-blocks-tsv {output.blocks} \
              --out-molecules-tsv {output.mols} \
              --min-alt-reads {params.min_alt_reads} \
              --min-cocover-reads {params.min_cocover} \
              --max-block-snps {params.max_snps} \
              --min-haplotype-reads {params.min_hap_reads}
            """

    rule test_haplotype_associations:
        input:
            haps = GENO_HAP_MOLS,
            mods = GENO_MOL_MODS
        output:
            hap_tx = GENO_HAP_TX,
            hap_mod = GENO_HAP_MOD
        params:
            script = "workflow/scripts/test_haplotype_associations.py",
            min_group = int(GENO_CFG.get("min_group_reads", 4)),
            min_hap = int(GENO_CFG.get("min_haplotype_reads", 4)),
            test = GENO_CFG.get("test", "auto"),
            pseudocount = float(GENO_CFG.get("pseudocount", 0.5)),
        conda:
            "../envs/modulator.yaml"
        shell:
            r"""
            mkdir -p $(dirname {output.hap_tx})
            python {params.script} \
              --molecule-haplotypes {input.haps} \
              --molecule-mods {input.mods} \
              --out-haplotype-transcript {output.hap_tx} \
              --out-haplotype-mod {output.hap_mod} \
              --min-haplotype-reads {params.min_hap} \
              --min-transcript-reads {params.min_group} \
              --min-total-reads {params.min_group} \
              --test {params.test} \
              --pseudocount {params.pseudocount}
            """
