################################################################################
# assemble_transcripts — constant outputs, conditional per-sample BAMs, robust path
################################################################################
import os
from snakemake.io import glob_wildcards

# ---- Config-derived constants (safe for 'output:') ----
PREFIX        = config["prefix"]
ASSEMBLE_DIR  = "results/assemble"
OUT_GTF       = f"{ASSEMBLE_DIR}/{PREFIX}.gtf"
OUT_CLASS     = f"{ASSEMBLE_DIR}/{PREFIX}_classification_summary.tsv"
OUT_TX        = f"{ASSEMBLE_DIR}/{PREFIX}_tx_counts.tsv"
OUT_PCA       = f"{ASSEMBLE_DIR}/{PREFIX}_tx_counts.pca.png"
OUT_STATS     = f"{ASSEMBLE_DIR}/{PREFIX}_per_sample_stats.tsv"
OUT_TX_READ_LENGTHS = f"{ASSEMBLE_DIR}/{PREFIX}_tx_assigned_read_lengths.tsv"
ZT_DIR        = f"{ASSEMBLE_DIR}/zt_tagged"

# ---- Discover sample names from BAMs at parse time ----
_BAMS_DIR = config["bams_dir"]
_BAM_GLOB = config.get("bam_glob", "*.bam")

_bam_paths = sorted(glob.glob(os.path.join(_BAMS_DIR, _BAM_GLOB)))
_SAMPLE_NAMES = [os.path.basename(p).replace(".bam", "") for p in _bam_paths]

# ---- Optional per-sample outputs if requested ----
_REQUIRE_ZT_SAMPLES = bool(config.get("assembler", {}).get("write_zt_tagged_sample_bams", False))
ZT_TAGGED_BAMS = [f"{ZT_DIR}/{s}.zt_tagged.bam" for s in _SAMPLE_NAMES] if _REQUIRE_ZT_SAMPLES else []

# ---- Final outputs (no callables) ----
ASSEMBLE_OUTPUTS = [
    OUT_GTF,
    OUT_CLASS,
    OUT_TX,
    OUT_PCA,
    OUT_STATS,
    OUT_TX_READ_LENGTHS,
    directory(ZT_DIR),   # always ensure the directory exists
] + ZT_TAGGED_BAMS

rule assemble_transcripts:
    """
    Assemble and annotate transcripts (ZN/ZT options).
    - Per-sample zt_tagged BAMs are required outputs only when
      assembler.write_zt_tagged_sample_bams = true.
    - Script path resolves whether invoked from repo root or workflow/.
    """
    input:
        dir = _BAMS_DIR,
        gtf = config["reference_gtf"]
    output:
        ASSEMBLE_OUTPUTS
    threads:
        int(config.get("threads", 1))
    params:
        # Basics
        glob                 = _BAM_GLOB,
        out_prefix           = f"{ASSEMBLE_DIR}/{PREFIX}",

        # Assembler params
        min_mapq             = config.get("assembler", {}).get("min_mapq", 10),
        min_introns_read     = config.get("assembler", {}).get("min_introns_read", 1),
        require_softclip3p   = config.get("assembler", {}).get("require_softclip3p", 0),
        apa_window           = config.get("assembler", {}).get("apa_window", 20),
        min_reads            = config.get("assembler", {}).get("min_reads", 40),
        min_frac             = config.get("assembler", {}).get("min_frac", 0.00),
        min_introns          = config.get("assembler", {}).get("min_introns", 1),
        min_polya_length     = config.get("assembler", {}).get("min_polya_length", 12),
        min_polya_purity     = config.get("assembler", {}).get("min_polya_purity", 0.5),
        polya_support_frac   = config.get("assembler", {}).get("polya_support_frac", 0.5),
        tes_match_tol        = config.get("assembler", {}).get("tes_match_tol", 25),
        exact_tes_tol        = config.get("assembler", {}).get("exact_tes_tol", 10),
        min_reads_per_sample_for_mod = config.get("assembler", {}).get("min_reads_per_sample_for_mod", 5),
        min_total_reads_for_mod      = config.get("assembler", {}).get("min_total_reads_for_mod", 20),

        # Conditional flags (rendered now)
        primary_only_flag    = "--primary-only" if config.get("assembler", {}).get("primary_only", True) else "",
        # write_zt_bams_flag   = "--write-zt-bams" if config.get("assembler", {}).get("write_zt_bams", False) else "",
        write_zt_tagged_sample_bams_flag = "--write-zt-tagged-sample-bams" if _REQUIRE_ZT_SAMPLES else "",
        # emit_modkit_manifest_flag      = "--emit-modkit-manifest" if config.get("assembler", {}).get("emit_modkit_manifest", False) else "",

        # Optional tes_window flag
        tes_window_flag = "" if str(config.get("assembler", {}).get("tes_window", "null")).lower() in ("", "none", "null")
                         else f"--tes-window {config.get('assembler', {}).get('tes_window')}"
    shell:
        r"""
        set -euo pipefail

        # Resolve script path robustly (works from repo root or workflow/)
        SCRIPT_A="{workflow.basedir}/scripts/assemble_transcripts.py"
        SCRIPT_B="{workflow.basedir}/workflow/scripts/assemble_transcripts.py"
        if [ -f "$SCRIPT_A" ]; then
            SCRIPT="$SCRIPT_A"
        elif [ -f "$SCRIPT_B" ]; then
            SCRIPT="$SCRIPT_B"
        else
            echo "ERROR: assemble_transcripts.py not found at $SCRIPT_A or $SCRIPT_B" >&2
            exit 2
        fi

        python "$SCRIPT" \
            --dir "{input.dir}" \
            --glob "{params.glob}" \
            --gtf "{input.gtf}" \
            --out-gtf "{output[0]}" \
            --out-prefix "{params.out_prefix}" \
            --threads "{threads}" \
            {params.primary_only_flag} \
            --min-mapq "{params.min_mapq}" \
            --min-introns-read "{params.min_introns_read}" \
            --require-softclip3p "{params.require_softclip3p}" \
            --apa-window "{params.apa_window}" \
            {params.tes_window_flag} \
            --min-reads "{params.min_reads}" \
            --min-frac "{params.min_frac}" \
            --min-introns "{params.min_introns}" \
            --min-polya-length "{params.min_polya_length}" \
            --min-polya-purity "{params.min_polya_purity}" \
            --polya-support-frac "{params.polya_support_frac}" \
            --tes-match-tol "{params.tes_match_tol}" \
            --exact-tes-tol "{params.exact_tes_tol}" \
            {params.write_zt_tagged_sample_bams_flag} \
            --min-reads-per-sample-for-mod "{params.min_reads_per_sample_for_mod}" \
            --min-total-reads-for-mod "{params.min_total_reads_for_mod}"

        # Ensure the directory() target always exists
        mkdir -p "{ZT_DIR}"
        """
