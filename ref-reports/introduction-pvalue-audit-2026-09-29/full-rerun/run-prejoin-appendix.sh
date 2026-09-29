#!/usr/bin/env bash
set -euo pipefail
project=/home/ed/projects/takeup
export OPENBLAS_NUM_THREADS=1 TAKEUP_THREADS=1 TAKEUP_DISTANCE_SPEC=realized
export TAKEUP_ANALYSIS_CONTEXT="$project/build/introduction-pvalue-audit-20260929/corrected-context.rds"
export TAKEUP_AUDIT_EXPECTED="$project/data/cluster_expected_dist.csv"
cd "$project/build/expected-distance-prejoin-audit-20260929/work/realized"
for script in paper-main-analysis-sample-attrition-table paper-knowledge-table-attrition-tables prior-deworming-robustness-table; do
 export TAKEUP_AUDIT_SCRIPT="scripts/appendix/$script.R"
 printf '%s started\n' "$script"
 Rscript --vanilla "$project/ref-reports/introduction-pvalue-audit-2026-09-29/full-rerun/launch-historical.R" > "$script.log" 2>&1
 printf '%s completed\n' "$script"
done
