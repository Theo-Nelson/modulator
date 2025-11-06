#!/usr/bin/env bash
# run_modkit_by_ZN.sh
# Run modkit pileup partitioned by ZN and per-mod summary, using assembler v8 outputs.

set -euo pipefail

### ====== USER VARIABLES (edit as needed) ======
# Tools & refs
MODKIT_BIN="${MODKIT_BIN:-modkit}"
REFERENCE_GENOME="${REFERENCE_GENOME:-hg38.fa}"
THREADS="${THREADS:-64}"

# Assembler outputs (from assemble_tes_variants_from_reads_v8.py)
GTF="${GTF:-fivegenes_readbacked_annot.gtf}"
SUMMARY="${SUMMARY:-fivegenes_readbacked_annot_classification_summary.tsv}"

# Inputs: ZT/ZN-tagged BAMs produced by the assembler
TAGGED_DIR="${TAGGED_DIR:-zt_tagged}"

# Outputs
MODKIT_OUT="${MODKIT_OUT:-modkit_out_ZN}"

# Mods & reference bases (order matters)
MODS=("17596" "a" "m" "17802" "69426" "19228" "19229" "19227")
REF_BASES=("A"     "A" "C" "T"     "A"     "C"     "G"     "T")
### ============================================

# Sanity checks
if ! command -v "$MODKIT_BIN" >/dev/null 2>&1; then
  echo "[ERROR] modkit not found in PATH (or MODKIT_BIN is wrong): $MODKIT_BIN" >&2
  exit 1
fi
if [[ ! -f "$REFERENCE_GENOME" ]]; then
  echo "[ERROR] Reference genome not found: $REFERENCE_GENOME" >&2
  exit 1
fi
if (( ${#MODS[@]} != ${#REF_BASES[@]} )); then
  echo "[ERROR] MODS and REF_BASES arrays must be the same length." >&2
  exit 1
fi
if [[ ! -d "$TAGGED_DIR" ]]; then
  echo "[ERROR] Tagged BAM directory not found: $TAGGED_DIR" >&2
  exit 1
fi

mkdir -p "$MODKIT_OUT"

# Make globs that don't match expand to nothing
shopt -s nullglob

# Loop over each tagged BAM
for ZN_BAM in "$TAGGED_DIR"/*.zt_tagged.bam; do
  SAMPLE_BASENAME="$(basename "$ZN_BAM" .zt_tagged.bam)"
  SAMPLE_OUT_DIR="${MODKIT_OUT}/${SAMPLE_BASENAME}"
  mkdir -p "$SAMPLE_OUT_DIR"

  echo "[INFO] Processing sample: $SAMPLE_BASENAME"
  echo "[INFO]   BAM: $ZN_BAM"
  echo "[INFO]   Output dir: $SAMPLE_OUT_DIR"

  # Per-mod pileup + summary
  for ((i=0; i<${#MODS[@]}; i++)); do
    MOD="${MODS[$i]}"
    REF="${REF_BASES[$i]}"

    OUT_PREFIX="${SAMPLE_OUT_DIR}/${SAMPLE_BASENAME}_${MOD}_filtered_mod"
    LOG_FILE="${OUT_PREFIX}_bed_log"

    echo "[INFO]   MOD=${MOD} (REF=${REF}) → ${OUT_PREFIX}.bed[._ZN].*"

    # PILEUP (partition by ZN → emits one file per ZN as suffix _<ZN>.bed)
    "$MODKIT_BIN" pileup "$ZN_BAM" \
      "${OUT_PREFIX}.bed" \
      --ref "$REFERENCE_GENOME" \
      --filter-threshold "${REF}":.8 \
      --mod-thresholds "${MOD}":0.99 \
      --partition-tag ZN \
      --max-depth 1000 \
      --interval-size 100000 \
      -t "$THREADS" \
      --log-filepath "$LOG_FILE"
  done
done

echo "[DONE] modkit pileup (ZN partitions) and summaries written under: $MODKIT_OUT"
echo "[HINT] Next: aggregate per-gene:"
echo "  python aggregate_modkit_by_ZN_per_gene_v1.py \\"
echo "    --modkit-dir \"$MODKIT_OUT\" \\"
echo "    --gtf \"$GTF\" \\"
echo "    --out-prefix modkit_by_transcript_ZN \\"
echo "    --min-cov 5 --write-gene-pivots --verbose"

