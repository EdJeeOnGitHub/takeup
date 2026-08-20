#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=indiv-fp-rs-prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=1-12:00:00
#SBATCH --output=temp/log/indiv-fp-rs-prod-%j.log
#SBATCH --error=temp/log/indiv-fp-rs-prod-%j.log

set -euo pipefail

CHAIN_ID=${CHAIN_ID:?Set CHAIN_ID to 1, 2, 3, or 4}
if [[ ! "${CHAIN_ID}" =~ ^[1-4]$ ]]; then
  echo "CHAIN_ID must be 1, 2, 3, or 4; got ${CHAIN_ID}." >&2
  exit 2
fi

INPUT_PATH=${INPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data}
OUTPUT_PATH=${OUTPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data/indiv-fp-reduce-sum-production}
STAN_PATH=${STAN_PATH:-stan_models_fit105}
STAN_FILE=${STAN_FILE:-takeup_struct_indiv_fp_reduce_sum.stan}
OUTPUT_BASENAME=${OUTPUT_BASENAME:-dist_fit106_INDIV_FP_REDUCE_SUM_PROD_chain${CHAIN_ID}}
THREADS_PER_CHAIN=${THREADS_PER_CHAIN:-8}
ITER_WARMUP=${ITER_WARMUP:-750}
ITER_SAMPLING=${ITER_SAMPLING:-750}
ADAPT_DELTA=${ADAPT_DELTA:-0.999}
MAX_TREEDEPTH=${MAX_TREEDEPTH:-12}
METRIC=${METRIC:-diag_e}
SEED=${SEED:-20260728}
REFRESH=${REFRESH:-25}
INIT_FILE=${INIT_FILE:-/project/akaring/takeup-data/temp-data/indiv-fp-warm-inits/init-chain-${CHAIN_ID}.json}

if [[ ! -f "${INIT_FILE}" ]]; then
  echo "Warm-start file not found: ${INIT_FILE}" >&2
  exit 2
fi

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}

Rscript --no-save --no-restore --no-init-file \
  scratch/sample-slim-individual-fp.R \
  "--input-path=${INPUT_PATH}" \
  "--output-path=${OUTPUT_PATH}" \
  "--stan-path=${STAN_PATH}" \
  "--stan-file=${STAN_FILE}" \
  "--output-basename=${OUTPUT_BASENAME}" \
  --chains=1 \
  --parallel-chains=1 \
  "--threads-per-chain=${THREADS_PER_CHAIN}" \
  "--iter-warmup=${ITER_WARMUP}" \
  "--iter-sampling=${ITER_SAMPLING}" \
  "--adapt-delta=${ADAPT_DELTA}" \
  "--max-treedepth=${MAX_TREEDEPTH}" \
  "--metric=${METRIC}" \
  "--seed=$((SEED + CHAIN_ID))" \
  "--refresh=${REFRESH}" \
  "--init-files=${INIT_FILE}"
