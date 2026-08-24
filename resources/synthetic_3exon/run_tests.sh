#!/usr/bin/env bash
# End-to-end regression suite for modulator, driven by the synthetic fixture.
#   1. repo smoke checks (assembler + genotype unit logic)
#   2. classify_diffs taxonomy unit test (all 13 structural categories)
#   3. full 15-stage pipeline on the synthetic dataset
#   4. validate_outputs.py -- 30 ground-truth assertions
# Exits non-zero if anything fails. Run from the repo root.
set -uo pipefail
ENV=${ENV:-/home/fs01/thn4005/.local/share/mamba/envs/modulator}   # override: ENV=/path ./run_tests.sh
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
PY="$ENV/bin/python"
HERE=resources/synthetic_3exon
rc=0

echo "### 1/4 repo smoke checks"
"$PY" workflow/scripts/regression_smoke_checks.py || rc=1
"$PY" workflow/scripts/genotype_regression_smoke_checks.py || rc=1

echo "### 2/4 classify_diffs taxonomy + snp-at-mod-base unit tests"
"$PY" "$HERE/test_classify_categories.py" || rc=1
"$PY" "$HERE/test_snp_at_mod_base.py" || rc=1
"$PY" "$HERE/test_dispersion_scaling.py" || rc=1
PYTHONPATH=src "$PY" "$HERE/test_scaling_features.py" || rc=1

echo "### 3/4 full pipeline on the synthetic dataset"
bash "$HERE/run_pipeline.sh" || rc=1

# stream-vs-sort aggregation parity: the pipeline step above just wrote results/modkit_zn +
# results/assemble/syn3exon.gtf, so this runs for real (not a vacuous skip) on the synthetic fixture.
echo "### 3b/4 aggregation engine parity (stream vs sort)"
"$PY" "$HERE/test_aggregate_engine_parity.py" || rc=1

echo "### 4/4 validate outputs against ground truth"
"$PY" "$HERE/validate_outputs.py" --results results --prefix syn3exon || rc=1

echo
if [ "$rc" -eq 0 ]; then echo "ALL REGRESSION CHECKS PASSED"; else echo "REGRESSION FAILURES (rc=$rc)"; fi
exit "$rc"
