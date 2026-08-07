#!/usr/bin/env bash
# Run the full modulator pipeline on the synthetic 3-exon dataset.
# Usage: bash resources/synthetic_3exon/run_pipeline.sh   (from the repo root)
set -uo pipefail
ENV=/home/fs01/thn4005/.local/share/mamba/envs/modulator
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
export TMPDIR="$ROOT/results/tmp"; mkdir -p "$TMPDIR"

PYTHONPATH=src PATH="$ENV/bin:$PATH" "$ENV/bin/python" -m modulator run \
  --config config/config.yaml --workdir "$ROOT" --jobs 4 ${RESUME:+--resume} \
  --set reference_fa=resources/synthetic_3exon/reference/synthetic_ref.fa \
        reference_gtf=resources/synthetic_3exon/reference/synthetic_ref.gtf \
        samplesheet=resources/synthetic_3exon/samples.tsv \
        bams_dir=resources/synthetic_3exon/bams \
        prefix=syn3exon \
        threads=8 \
        genotype.enable=true \
        genotype.jobs=4 \
        genotype.mod_sites_require_snp_link=false \
        hierarchical_stoich.enable=true
