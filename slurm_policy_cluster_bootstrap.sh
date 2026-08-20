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
NUM_CORES=${NUM_CORES:-12}
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
TABLE_PATH=${TABLE_PATH:-presentations/tables/fit105/optim-summ-exponential-cluster-weights.tex}
POLICY_REVIEW_OUTPUT=${POLICY_REVIEW_OUTPUT:-ref-reports/policy-cost-sensitivity}
POLICY_WORK_PATH=${POLICY_WORK_PATH:-temp-data/policy-cost-sensitivity}

cd "${PROJECT_ROOT}"
module load -f R/4.2.0
if [[ -n "${GUROBI_MODULE:-}" ]]; then
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
    Rscript --no-save --no-restore optim/prepare-policy-cluster-bootstrap.R \
      "--weighted-path=${WEIGHTED_PATH}" "--output-path=${OUTPUT_PATH}" \
      "--num-replicates=${NUM_REPLICATES}" "--method=${WEIGHT_METHOD}"
    ;;
  predict)
    Rscript --no-save --no-restore optim/predict-policy-cluster-bootstrap.R \
      "--parameter-csv=${OUTPUT_PATH}/policy-bootstrap-parameters.csv" \
      "--distance-data=${DISTANCE_DATA}" "--output-path=${OUTPUT_PATH}" \
      --distance-cap=3500 "--num-cores=${NUM_CORES}" \
      "--num-replicates=${NUM_REPLICATES}"
    ;;
  optimize)
    : "${SLURM_ARRAY_TASK_ID:?Optimize requires scenario array 1-5}"
    Rscript --no-save --no-restore optim/optimize-policy-cluster-bootstrap.R \
      "--input-path=${OUTPUT_PATH}" "--target-csv=${TARGET_CSV}" \
      "--scenario-id=${SLURM_ARRAY_TASK_ID}" \
      "--num-replicates=${NUM_REPLICATES}"
    ;;
  summarize)
    Rscript --no-save --no-restore optim/summarize-policy-cluster-bootstrap.R \
      "--input-path=${OUTPUT_PATH}" "--table-path=${TABLE_PATH}" \
      "--num-replicates=${NUM_REPLICATES}" "--method=${WEIGHT_METHOD}"
    ;;
  population)
    Rscript --no-save --no-restore optim/run-policy-population-cost.R \
      "--parameter-csv=${OUTPUT_PATH}/policy-bootstrap-parameters.csv" \
      --parameter-type=canonical --analysis-id=exponential-cluster-weights \
      "--distance-data=${DISTANCE_DATA}" \
      "--output-path=${POLICY_REVIEW_OUTPUT}" \
      "--work-path=${POLICY_WORK_PATH}" "--cores=${NUM_CORES}" \
      "--max-draws=${NUM_REPLICATES}" --solver=auto
    ;;
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
