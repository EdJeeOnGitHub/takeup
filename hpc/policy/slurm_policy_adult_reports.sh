#!/usr/bin/env bash
set -euo pipefail
: "${REPO_ROOT:?}" "${RUN_ROOT:?}" "${DISTANCE_DATA:?}" "${POLICY_CENSUS:?}" "${POLICY_GUROBI_ROOT:?}" "${STAGE:?}"
module load R/4.2.0
cd "$REPO_ROOT"
export PATH="$POLICY_GUROBI_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$POLICY_GUROBI_ROOT/lib:${LD_LIBRARY_PATH:-}"
export GRB_LICENSE_FILE="$POLICY_GUROBI_ROOT/gurobi.lic"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
REPORT_ROOT=${REPORT_ROOT:-$RUN_ROOT/reports}
mkdir -p "$REPORT_ROOT/tables" "$REPORT_ROOT/figures"
benchmark="$RUN_ROOT/cap-3500/benchmark"
case "$STAGE" in
  cost-baseline|cost-bootstrap)
    input=$benchmark
    analysis=baseline-posterior
    pooling=0
    if [[ "$STAGE" == cost-bootstrap ]]; then
      input="$RUN_ROOT/cap-3500/cluster-weighted"
      analysis=exponential-cluster-weights
      pooling=0,0.5,1
    fi
    test -s "$input/policy-collection-audit.csv"
    Rscript scripts/policy/run-population-cost.R \
      "--parameter-csv=$input/policy-model-parameters.csv" --parameter-type=canonical \
      "--analysis-id=$analysis" "--distance-data=$DISTANCE_DATA" \
      "--output-path=$REPORT_ROOT" "--work-path=${SLURM_TMPDIR:-/tmp}" \
      "--cores=${SLURM_CPUS_PER_TASK:-8}" --solver=gurobi \
      "--allocation-path=$input" "--pooling-rhos=$pooling" --include-legacy=true
    ;;
  median)
    Rscript scripts/policy/generate-distance-cap-table.R \
      "--parameter-csv=$benchmark/policy-model-parameters.csv" --parameter-type=canonical \
      "--distance-data=$DISTANCE_DATA" "--csv-path=$REPORT_ROOT/policy-distance-cap-diagnostics.csv" \
      "--policy-path=$benchmark" \
      "--table-path=$REPORT_ROOT/tables/policy-distance-cap-feasibility.tex"
    Rscript scripts/policy/render-paper-figures.R \
      "--parameter-csv=$benchmark/policy-model-parameters.csv" --parameter-draw=median \
      "--policy-path=$benchmark" "--median-allocation-path=$REPORT_ROOT/median-allocations" \
      "--distance-data=$DISTANCE_DATA" "--output-path=$REPORT_ROOT/figures"
    ;;
  assemble)
    Rscript scripts/policy/assemble-adult-policy-overview.R "--run-root=$RUN_ROOT" "--output-path=$REPORT_ROOT"
    Rscript scripts/policy/write-policy-amplification-contrasts.R "--report-path=$REPORT_ROOT"
    Rscript scripts/policy/render-policy-cap-table.R "--report-path=$REPORT_ROOT"
    Rscript scripts/policy/validate-population-cost.R "--input-path=$REPORT_ROOT" "--distance-data=$DISTANCE_DATA"
    Rscript scripts/policy/assemble-population-cost.R "--input-path=$REPORT_ROOT" \
      "--table-path=$REPORT_ROOT/tables" "--figure-path=$REPORT_ROOT/figures"
    Rscript scripts/policy/audit-complete-model-robustness.R "--package-path=$RUN_ROOT/cap-3500"
    Rscript scripts/policy/render-model-scenario-table.R "--input-root=$RUN_ROOT/cap-3500" \
      "--output=$REPORT_ROOT/tables/optim-policy-model-scenarios.tex"
    ;;
  *) echo "Unknown reporting stage: $STAGE" >&2; exit 2 ;;
esac
touch "$REPORT_ROOT/_SUCCESS-$STAGE"
