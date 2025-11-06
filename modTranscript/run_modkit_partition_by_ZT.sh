#!/usr/bin/env bash
set -euo pipefail

# Config (override via env or edit here)
MODKIT_BIN="${MODKIT_BIN:-modkit}"
REFERENCE_GENOME="${REFERENCE_GENOME:-hg38.fa}"   # must have .fai index
THREADS="${THREADS:-8}"

# Inputs/outputs
IN_DIR="${1:-zt_tagged}"      # where assemble_v6 wrote <sample>.zt_tagged.bam
OUT_DIR="${2:-modkit_out}"    # output directory (modkit creates one .bed per ZT code)
mkdir -p "${OUT_DIR}"

# Per-base filter thresholds (canonical base >= 0.8), exactly like your snippet
REF_THR=( --filter-threshold A:.8 --filter-threshold C:.8 --filter-threshold G:.8 --filter-threshold T:.8 )

# Per-mod pass thresholds (0.99) for the 8 mods you listed
MOD_THR=(
  --mod-thresholds 17596:0.99
  --mod-thresholds a:0.99
  --mod-thresholds m:0.99
  --mod-thresholds 17802:0.99
  --mod-thresholds 69426:0.99
  --mod-thresholds 19228:0.99
  --mod-thresholds 19229:0.99
  --mod-thresholds 19227:0.99
)

for BAM in "${IN_DIR}"/*.zt_tagged.bam; do
  [ -e "$BAM" ] || { echo "No BAMs found in ${IN_DIR}"; exit 1; }
  SAMPLE=$(basename "$BAM" .zt_tagged.bam)
  LOG_FILE="${OUT_DIR}/${SAMPLE}.log"

  "${MODKIT_BIN}" pileup "${BAM}" "${OUT_DIR}" \
    --ref "${REFERENCE_GENOME}" \
    "${REF_THR[@]}" \
    "${MOD_THR[@]}" \
    --log-filepath "${LOG_FILE}" \
    --max-depth 1000 \
    --interval-size 100000 \
    --prefix "${SAMPLE}" \
    --partition-tag ZT \
    -t "${THREADS}"
done

echo "[OK] Finished modkit pileup; outputs in ${OUT_DIR}"

