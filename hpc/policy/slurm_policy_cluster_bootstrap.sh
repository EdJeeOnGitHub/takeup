#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=policy-cluster-bs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=28G
#SBATCH --time=04:00:00
#SBATCH --output=temp/log/policy-cluster-bs-%A-%a.log
#SBATCH --error=temp/log/policy-cluster-bs-%A-%a.log

set -euo pipefail

STAGE=${STAGE:?Set STAGE to prepare, predict, optimize, summarize, or population}
PROJECT_ROOT=${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}
NUM_REPLICATES=${NUM_REPLICATES:-999}
WEIGHT_METHOD=${WEIGHT_METHOD:-exponential}
DISTANCE_DEFINITION=${DISTANCE_DEFINITION:-assigned}
DISTANCE_CAP=${DISTANCE_CAP:-3500}
POPULATION_WEIGHTING=${POPULATION_WEIGHTING:-equal-community}
TARGET_MODE=${TARGET_MODE:-legacy-fixed}
CENSUS_DATA=${CENSUS_DATA:-${POLICY_CENSUS:-${PROJECT_ROOT}/data/takeup_census.RData}}
export POLICY_CENSUS="$CENSUS_DATA"
NUM_CORES=${NUM_CORES:-12}
optimizer_cpu_budget=${SLURM_CPUS_PER_TASK:-${NUM_CORES}}
OPTIMIZE_CORES=${OPTIMIZE_CORES:-$(( optimizer_cpu_budget < 8 ? optimizer_cpu_budget : 8 ))}
DRAW_BATCH_SIZE=${DRAW_BATCH_SIZE:-25}
SOLVER_THREADS=${SOLVER_THREADS:-1}
SOLVER_SEED=${SOLVER_SEED:-0}
POLICY_SOLVER=${POLICY_SOLVER:-gurobi}
POLICY_SCRATCH=${POLICY_SCRATCH:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

MODEL=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
OUTPUT_PATH=${OUTPUT_PATH:-optim/data/${MODEL}/agg-full-many-pots-exponential-cluster-weights}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
if [[ -z "${WEIGHTED_PATH:-}" ]]; then
  if [[ -d ${ANALYSIS_ROOT}/main-core-exponential-cluster-weight-999 ]]; then
    WEIGHTED_PATH=${ANALYSIS_ROOT}/main-core-exponential-cluster-weight-999
  else
    WEIGHTED_PATH=${ANALYSIS_ROOT}/main-core-weighted-modes
  fi
fi
DISTANCE_DATA=${DISTANCE_DATA:-optim/data/full-many-pots-experiment.rds}
TARGET_CSV=${TARGET_CSV:-optim/data/${MODEL}/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv}
if [[ "$TARGET_MODE" == draw-specific-experimental-control ]]; then
  TARGET_CSV="$OUTPUT_PATH/policy-experimental-targets.csv"
fi
TABLE_PATH=${TABLE_PATH:-presentations/tables/fit105/optim-summ-exponential-cluster-weights.tex}
POLICY_REVIEW_OUTPUT=${POLICY_REVIEW_OUTPUT:-ref-reports/policy-cost-sensitivity}
POLICY_WORK_PATH=${POLICY_WORK_PATH:-temp-data/policy-cost-sensitivity}

cd "${PROJECT_ROOT}"
module load -f R/4.2.0
if [[ -n "${POLICY_GUROBI_ROOT:-}" ]]; then
  export PATH="$POLICY_GUROBI_ROOT/bin:$PATH"
  export LD_LIBRARY_PATH="$POLICY_GUROBI_ROOT/lib:${LD_LIBRARY_PATH:-}"
  export GRB_LICENSE_FILE="$POLICY_GUROBI_ROOT/gurobi.lic"
elif [[ -n "${GUROBI_MODULE:-}" ]]; then
  module load -f "${GUROBI_MODULE}"
elif [[ "${SLURM_JOB_PARTITION:-}" == "caslake" ]]; then
  module load -f gurobi/11.0
else
  module load -f gurobi/9.2
fi
if [[ -z "${R_LIBS_USER:-}" ]]; then
  if [[ -d /home/edjee/R/x86_64-pc-linux-gnu-library/4.2 ]]; then
    R_LIBS_USER=/home/edjee/R/x86_64-pc-linux-gnu-library/4.2
  else
    R_LIBS_USER=/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu
  fi
fi
export R_LIBS_USER
mkdir -p temp/log "${OUTPUT_PATH}"

case "${STAGE}" in
  prepare)
    Rscript --no-save --no-restore scripts/policy/prepare-cluster-bootstrap.R \
      "--weighted-path=${WEIGHTED_PATH}" "--output-path=${OUTPUT_PATH}" \
      "--num-replicates=${NUM_REPLICATES}" "--method=${WEIGHT_METHOD}" \
      "--distance-definition=${DISTANCE_DEFINITION}"
    ;;
  predict)
    Rscript --no-save --no-restore scripts/policy/predict-cluster-bootstrap.R \
      "--parameter-csv=${OUTPUT_PATH}/policy-bootstrap-parameters.csv" \
      "--distance-data=${DISTANCE_DATA}" "--output-path=${OUTPUT_PATH}" \
      "--distance-cap=${DISTANCE_CAP}" "--num-cores=${NUM_CORES}" \
      "--population-weighting=${POPULATION_WEIGHTING}" "--census-data=${CENSUS_DATA}" \
      "--num-replicates=${NUM_REPLICATES}"
    ;;
  optimize)
    : "${SLURM_ARRAY_TASK_ID:?Optimize requires scenario array 1-5}"
    Rscript --no-save --no-restore scripts/policy/optimize-cluster-bootstrap.R \
      "--input-path=${OUTPUT_PATH}" "--target-csv=${TARGET_CSV}" \
      "--population-weighting=${POPULATION_WEIGHTING}" "--target-mode=${TARGET_MODE}" \
      "--distance-data=${DISTANCE_DATA}" "--census-data=${CENSUS_DATA}" \
      "--num-cores=${OPTIMIZE_CORES}" "--draw-batch-size=${DRAW_BATCH_SIZE}" \
      "--solver=${POLICY_SOLVER}" "--solver-threads=${SOLVER_THREADS}" \
      "--solver-seed=${SOLVER_SEED}" "--scratch-path=${POLICY_SCRATCH}" \
      "--scenario-id=${SLURM_ARRAY_TASK_ID}" \
      "--num-replicates=${NUM_REPLICATES}"
    ;;
  summarize)
    Rscript --no-save --no-restore scripts/policy/summarize-cluster-bootstrap.R \
      "--input-path=${OUTPUT_PATH}" "--table-path=${TABLE_PATH}" \
      "--num-replicates=${NUM_REPLICATES}" "--method=${WEIGHT_METHOD}"
    ;;
  population)
    REUSE_OPTIONS=()
    if [[ "$POPULATION_WEIGHTING" == adult-census ]]; then REUSE_OPTIONS+=("--allocation-path=$OUTPUT_PATH"); fi
    Rscript --no-save --no-restore scripts/policy/run-population-cost.R \
      "--parameter-csv=${OUTPUT_PATH}/policy-bootstrap-parameters.csv" \
      --parameter-type=canonical --analysis-id=exponential-cluster-weights \
      "--distance-data=${DISTANCE_DATA}" \
      "--output-path=${POLICY_REVIEW_OUTPUT}" \
      "--work-path=${POLICY_WORK_PATH}" "--cores=${NUM_CORES}" \
      "--max-draws=${NUM_REPLICATES}" "--solver=${POLICY_SOLVER}" \
      "--solver-threads=${SOLVER_THREADS}" "--solver-seed=${SOLVER_SEED}" \
      ${REUSE_OPTIONS[@]+"${REUSE_OPTIONS[@]}"}
    ;;
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
