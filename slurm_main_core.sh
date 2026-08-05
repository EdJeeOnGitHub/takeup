#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=main-core
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=1-12:00:00
#SBATCH --output=temp/log/main-core-%A_%a.log
#SBATCH --error=temp/log/main-core-%A_%a.log

set -euo pipefail

CHAIN_ID=${CHAIN_ID:-${SLURM_ARRAY_TASK_ID:-}}
if [[ -z "${CHAIN_ID}" ]]; then
  echo "Set CHAIN_ID or submit this script with --array=1-4." >&2
  exit 2
fi
if [[ ! "${CHAIN_ID}" =~ ^[1-4]$ ]]; then
  echo "CHAIN_ID must be 1, 2, 3, or 4; got ${CHAIN_ID}." >&2
  exit 2
fi

MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP}
INPUT_PATH=${INPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data}
OUTPUT_PATH=${OUTPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-production}
STAN_PATH=${STAN_PATH:-stan_models}
STAN_FILE=${STAN_FILE:-takeup_struct_main_core.stan}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.33.1}
OUTPUT_BASENAME=${OUTPUT_BASENAME:-dist_fit106_MAIN_CORE_chain${CHAIN_ID}}
THREADS_PER_CHAIN=${THREADS_PER_CHAIN:-8}
ITER_WARMUP=${ITER_WARMUP:-1000}
ITER_SAMPLING=${ITER_SAMPLING:-1000}
ADAPT_DELTA=${ADAPT_DELTA:-0.999}
MAX_TREEDEPTH=${MAX_TREEDEPTH:-12}
METRIC=${METRIC:-diag_e}
SEED=${SEED:-20260801}
REFRESH=${REFRESH:-25}
INIT_FILE=${INIT_FILE:-}
CLUSTER_WEIGHT_FILE=${CLUSTER_WEIGHT_FILE:-}
USE_CORE_CLUSTER_SHOCK=${USE_CORE_CLUSTER_SHOCK:-0}
CORE_CLUSTER_SHOCK_SD_PRIOR=${CORE_CLUSTER_SHOCK_SD_PRIOR:-0.1}
CORE_LAMBDA_STRUCTURE=${CORE_LAMBDA_STRUCTURE:-0}
CORE_LAMBDA_LOG_RATIO_SD_PRIOR=${CORE_LAMBDA_LOG_RATIO_SD_PRIOR:-0.25}
CORE_OBSERVATION_MODEL=${CORE_OBSERVATION_MODEL:-0}
CORE_RECOGNITION_STRUCTURE=${CORE_RECOGNITION_STRUCTURE:-0}
CORE_REPORT_STRUCTURE=${CORE_REPORT_STRUCTURE:-0}
CORE_REPORT_ARM_DIST_HIERARCHICAL=${CORE_REPORT_ARM_DIST_HIERARCHICAL:-0}
CORE_REPORT_ARM_DIST_PRIOR_SCALE=${CORE_REPORT_ARM_DIST_PRIOR_SCALE:-0.25}

mkdir -p temp/log "${OUTPUT_PATH}"
module load gcc/10.2.0
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

args=(
  "--model=${MODEL}"
  "--input-path=${INPUT_PATH}"
  "--output-path=${OUTPUT_PATH}"
  "--stan-path=${STAN_PATH}"
  "--stan-file=${STAN_FILE}"
  "--cmdstan-path=${CMDSTAN_PATH}"
  "--output-basename=${OUTPUT_BASENAME}"
  --chains=1
  --parallel-chains=1
  "--threads-per-chain=${THREADS_PER_CHAIN}"
  "--iter-warmup=${ITER_WARMUP}"
  "--iter-sampling=${ITER_SAMPLING}"
  "--adapt-delta=${ADAPT_DELTA}"
  "--max-treedepth=${MAX_TREEDEPTH}"
  "--metric=${METRIC}"
  "--seed=$((SEED + CHAIN_ID))"
  "--refresh=${REFRESH}"
  "--use-core-cluster-shock=${USE_CORE_CLUSTER_SHOCK}"
  "--core-cluster-shock-sd-prior=${CORE_CLUSTER_SHOCK_SD_PRIOR}"
  "--core-lambda-structure=${CORE_LAMBDA_STRUCTURE}"
  "--core-lambda-log-ratio-sd-prior=${CORE_LAMBDA_LOG_RATIO_SD_PRIOR}"
  "--core-observation-model=${CORE_OBSERVATION_MODEL}"
  "--core-recognition-structure=${CORE_RECOGNITION_STRUCTURE}"
  "--core-report-structure=${CORE_REPORT_STRUCTURE}"
  "--core-report-arm-dist-hierarchical=${CORE_REPORT_ARM_DIST_HIERARCHICAL}"
  "--core-report-arm-dist-prior-scale=${CORE_REPORT_ARM_DIST_PRIOR_SCALE}"
)
if [[ -n "${INIT_FILE}" ]]; then
  if [[ ! -f "${INIT_FILE}" ]]; then
    echo "Initial-value file not found: ${INIT_FILE}" >&2
    exit 2
  fi
  args+=("--init-files=${INIT_FILE}")
fi
if [[ -n "${CLUSTER_WEIGHT_FILE}" ]]; then
  if [[ ! -f "${CLUSTER_WEIGHT_FILE}" ]]; then
    echo "Cluster-weight file not found: ${CLUSTER_WEIGHT_FILE}" >&2
    exit 2
  fi
  args+=("--cluster-weight-file=${CLUSTER_WEIGHT_FILE}")
fi

Rscript --no-save --no-restore --no-init-file \
  scratch/sample-slim-individual-fp.R "${args[@]}"
