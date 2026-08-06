#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=core-finite-mix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=temp/log/core-finite-mix-%A-%a.log
#SBATCH --error=temp/log/core-finite-mix-%A-%a.log

set -euo pipefail
STAGE=${STAGE:?Set STAGE to mode, sample, gq, or summarize}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
OUTPUT_PATH=${OUTPUT_PATH:-${ANALYSIS_ROOT}/main-core-finite-mixture-robustness}
MODE_PATH=${MODE_PATH:-${OUTPUT_PATH}/mode/finite-mixture}
BASELINE_GQ_DIR=${BASELINE_GQ_DIR:-${ANALYSIS_ROOT}/main-core-lambda-identification/gq/common}
STAN_PATH=${STAN_PATH:-stan_models}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.33.1}
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP}
ITER_WARMUP=${ITER_WARMUP:-400}
ITER_SAMPLING=${ITER_SAMPLING:-400}

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

case "${STAGE}" in
  mode)
    Rscript --no-save --no-restore --no-init-file \
      scratch/fit-main-core-weighted-mode.R \
      "--workspace=${WORKSPACE}" "--model=${MODEL}" \
      "--output-path=${OUTPUT_PATH}/mode" --label=finite-mixture \
      --core-type-distribution=2 \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --threads=8 --seed=20260921
    ;;
  sample)
    : "${SLURM_ARRAY_TASK_ID:?sample requires --array=1-4}"
    chain=${SLURM_ARRAY_TASK_ID}
    mkdir -p "${OUTPUT_PATH}/fits/finite-mixture"
    Rscript --no-save --no-restore --no-init-file \
      scratch/sample-slim-individual-fp.R \
      "--model=${MODEL}" "--input-path=$(dirname "${WORKSPACE}")" \
      "--output-path=${OUTPUT_PATH}/fits/finite-mixture" \
      "--stan-path=${STAN_PATH}" --stan-file=takeup_struct_main_core.stan \
      "--cmdstan-path=${CMDSTAN_PATH}" \
      "--output-basename=finite-mixture-chain${chain}" \
      --chains=1 --parallel-chains=1 --threads-per-chain=8 \
      "--iter-warmup=${ITER_WARMUP}" "--iter-sampling=${ITER_SAMPLING}" \
      --adapt-delta=0.999 --max-treedepth=12 --metric=diag_e \
      "--seed=$((20260940 + chain))" \
      "--init-files=${MODE_PATH}/mode-init.json" \
      --core-type-distribution=2
    ;;
  gq)
    fit_csvs=$(find "${OUTPUT_PATH}/fits/finite-mixture" -maxdepth 1 -type f \
      -name 'finite-mixture-chain*-*.csv' | sort | paste -sd, -)
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-compact-gq.R \
      "--workspace=${WORKSPACE}" "--model=${MODEL}" \
      "--fit-csvs=${fit_csvs}" \
      "--output-path=${OUTPUT_PATH}/gq/finite-mixture" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --core-type-distribution=2 \
      "--force-recompile=${FORCE_RECOMPILE:-1}" \
      --threads=2 --parallel-chains=4
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scratch/summarize-main-core-finite-mixture.R \
      "--baseline-gq-dir=${BASELINE_GQ_DIR}" \
      "--mixture-fit-dir=${OUTPUT_PATH}/fits/finite-mixture" \
      "--mixture-gq-dir=${OUTPUT_PATH}/gq/finite-mixture" \
      "--output-path=${OUTPUT_PATH}"
    ;;
  *) echo "Unknown STAGE=${STAGE}" >&2; exit 2 ;;
esac
