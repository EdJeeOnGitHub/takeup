#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=indiv-fp-slim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-12:00:00
#SBATCH --output=temp/log/indiv-fp-slim-%j.log
#SBATCH --error=temp/log/indiv-fp-slim-%j.log

set -euo pipefail

INPUT_PATH=${INPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data}
OUTPUT_PATH=${OUTPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data}
STAN_PATH=${STAN_PATH:-stan_models_fit105}
STAN_FILE=${STAN_FILE:-takeup_struct_indiv_fp_slim.stan}
OUTPUT_BASENAME=${OUTPUT_BASENAME:-dist_fit106_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP_SLIM}
CHAINS=${CHAINS:-4}
PARALLEL_CHAINS=${PARALLEL_CHAINS:-4}
THREADS_PER_CHAIN=${THREADS_PER_CHAIN:-2}
ITER_WARMUP=${ITER_WARMUP:-1000}
ITER_SAMPLING=${ITER_SAMPLING:-1000}
ADAPT_DELTA=${ADAPT_DELTA:-0.999}
MAX_TREEDEPTH=${MAX_TREEDEPTH:-12}
METRIC=${METRIC:-dense_e}
SEED=${SEED:-20260724}
REFRESH=${REFRESH:-25}
INIT_FILES=${INIT_FILES:-}
FIT_MODEL=${FIT_MODEL:-}
FIT_DIST_MODEL=${FIT_DIST_MODEL:-}
FIT_BELIEFS_MODEL=${FIT_BELIEFS_MODEL:-}
FIT_WTP_MODEL=${FIT_WTP_MODEL:-}

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}

args=(
  "--input-path=${INPUT_PATH}" \
  "--output-path=${OUTPUT_PATH}" \
  "--stan-path=${STAN_PATH}" \
  "--stan-file=${STAN_FILE}" \
  "--output-basename=${OUTPUT_BASENAME}" \
  "--chains=${CHAINS}" \
  "--parallel-chains=${PARALLEL_CHAINS}" \
  "--threads-per-chain=${THREADS_PER_CHAIN}" \
  "--iter-warmup=${ITER_WARMUP}" \
  "--iter-sampling=${ITER_SAMPLING}" \
  "--adapt-delta=${ADAPT_DELTA}" \
  "--max-treedepth=${MAX_TREEDEPTH}" \
  "--metric=${METRIC}" \
  "--seed=${SEED}" \
  "--refresh=${REFRESH}"
)

if [[ -n "${INIT_FILES}" ]]; then
  args+=("--init-files=${INIT_FILES}")
fi
if [[ -n "${FIT_MODEL}" ]]; then
  args+=("--fit-model=${FIT_MODEL}")
fi
if [[ -n "${FIT_DIST_MODEL}" ]]; then
  args+=("--fit-dist-model=${FIT_DIST_MODEL}")
fi
if [[ -n "${FIT_BELIEFS_MODEL}" ]]; then
  args+=("--fit-beliefs-model=${FIT_BELIEFS_MODEL}")
fi
if [[ -n "${FIT_WTP_MODEL}" ]]; then
  args+=("--fit-wtp-model=${FIT_WTP_MODEL}")
fi

Rscript --no-save --no-restore --no-init-file \
  scratch/sample-slim-individual-fp.R "${args[@]}"
