#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=core-student-t
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1-12:00:00
#SBATCH --output=temp/log/core-student-t-%A-%a.log
#SBATCH --error=temp/log/core-student-t-%A-%a.log

set -euo pipefail
STAGE=${STAGE:?Set STAGE to mode, sample, gq, or summarize}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
OUTPUT_PATH=${OUTPUT_PATH:-${ANALYSIS_ROOT}/main-core-student-t-robustness}
MODE_PATH=${MODE_PATH:-${OUTPUT_PATH}/mode/student-t5}
BASELINE_GQ_DIR=${BASELINE_GQ_DIR:-${ANALYSIS_ROOT}/main-core-lambda-identification/gq/common}
STAN_PATH=${STAN_PATH:-stan_models_fit105}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.33.1}
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP}
DF=${DF:-5}
COMPONENTS=${COMPONENTS:-12}
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
      "--output-path=${OUTPUT_PATH}/mode" --label=student-t5 \
      --core-type-distribution=1 "--core-student-t-df=${DF}" \
      "--core-student-t-components=${COMPONENTS}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --threads=8 --seed=20260831
    ;;
  sample)
    : "${SLURM_ARRAY_TASK_ID:?sample requires --array=1-4}"
    chain=${SLURM_ARRAY_TASK_ID}
    mkdir -p "${OUTPUT_PATH}/fits/student-t5"
    Rscript --no-save --no-restore --no-init-file \
      scratch/sample-slim-individual-fp.R \
      "--model=${MODEL}" "--input-path=$(dirname "${WORKSPACE}")" \
      "--output-path=${OUTPUT_PATH}/fits/student-t5" \
      "--stan-path=${STAN_PATH}" --stan-file=takeup_struct_main_core.stan \
      "--cmdstan-path=${CMDSTAN_PATH}" \
      "--output-basename=student-t5-chain${chain}" \
      --chains=1 --parallel-chains=1 --threads-per-chain=8 \
      "--iter-warmup=${ITER_WARMUP}" "--iter-sampling=${ITER_SAMPLING}" \
      --adapt-delta=0.999 --max-treedepth=12 --metric=diag_e \
      "--seed=$((20260900 + chain))" \
      "--init-files=${MODE_PATH}/mode-init.json" \
      --core-type-distribution=1 "--core-student-t-df=${DF}" \
      "--core-student-t-components=${COMPONENTS}"
    ;;
  gq)
    fit_csvs=$(find "${OUTPUT_PATH}/fits/student-t5" -maxdepth 1 -type f \
      -name 'student-t5-chain*-*.csv' | sort | paste -sd, -)
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-compact-gq.R \
      "--workspace=${WORKSPACE}" "--model=${MODEL}" \
      "--fit-csvs=${fit_csvs}" "--output-path=${OUTPUT_PATH}/gq/student-t5" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --core-type-distribution=1 "--core-student-t-df=${DF}" \
      "--core-student-t-components=${COMPONENTS}" --threads=2 --parallel-chains=4
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scratch/summarize-main-core-student-t.R \
      "--baseline-gq-dir=${BASELINE_GQ_DIR}" \
      "--student-fit-dir=${OUTPUT_PATH}/fits/student-t5" \
      "--student-gq-dir=${OUTPUT_PATH}/gq/student-t5" \
      "--output-path=${OUTPUT_PATH}"
    ;;
  *) echo "Unknown STAGE=${STAGE}" >&2; exit 2 ;;
esac
