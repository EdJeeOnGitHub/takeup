#!/usr/bin/env bash
set -euo pipefail
project=/home/ed/projects/takeup
export OPENBLAS_NUM_THREADS=1 TAKEUP_THREADS=1 TAKEUP_DISTANCE_SPEC=realized
export TAKEUP_ANALYSIS_CONTEXT="$project/build/introduction-pvalue-audit-20260929/corrected-context.rds"
for version in full historical; do
  cd "$project/build/expected-distance-${version}-audit-20260929/work/realized"
  launcher="$project/ref-reports/introduction-pvalue-audit-2026-09-29/full-rerun/launch.R"
  if [[ "$version" == historical ]]; then
    launcher="$project/ref-reports/introduction-pvalue-audit-2026-09-29/full-rerun/launch-historical.R"
  fi
  for script in paper-main-analysis-sample-attrition-table paper-knowledge-table-attrition-tables prior-deworming-robustness-table; do
    if [[ -f "$script.completed" ]]; then continue; fi
    export TAKEUP_AUDIT_SCRIPT="scripts/appendix/$script.R"
    printf '%s %s started\n' "$version" "$script"
    Rscript --vanilla "$launcher" > "$script.log" 2>&1
    touch "$script.completed"
    printf '%s %s completed\n' "$version" "$script"
  done
done
