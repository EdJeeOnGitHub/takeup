#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=lambda-profile
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --time=06:00:00
#SBATCH --output=temp/log/lambda-profile-%A-%a.log
#SBATCH --error=temp/log/lambda-profile-%A-%a.log

set -euo pipefail
STAGE=${STAGE:?Set STAGE to prepare, profile, or summarize}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
ROOT_OUTPUT=${ROOT_OUTPUT:-${ANALYSIS_ROOT}/main-core-lambda-identification}
OUTPUT_PATH=${OUTPUT_PATH:-${ROOT_OUTPUT}/profile}
MANIFEST=${MANIFEST:-${OUTPUT_PATH}/lambda-profile-manifest.csv}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP}
STAN_PATH=${STAN_PATH:-stan_models_fit105}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.33.1}
INIT_JSON=${INIT_JSON:-${ROOT_OUTPUT}/modes/grouped-sd0p25/mode-init.json}

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

case "${STAGE}" in
  prepare)
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-lambda-profile-manifest.R \
      "--output-path=${OUTPUT_PATH}" --grid-min=-3 --grid-max=3 --grid-step=0.1
    ;;
  profile)
    : "${SLURM_ARRAY_TASK_ID:?profile requires a Slurm array}"
    Rscript --no-save --no-restore --no-init-file \
      scratch/profile-main-core-lambda.R \
      "--manifest=${MANIFEST}" "--task-id=${SLURM_ARRAY_TASK_ID}" \
      "--workspace=${WORKSPACE}" "--model=${MODEL}" \
      "--output-path=${OUTPUT_PATH}" "--init-json=${INIT_JSON}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" --threads=8
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scratch/summarize-main-core-lambda-profile.R \
      "--manifest=${MANIFEST}" "--output-path=${OUTPUT_PATH}"
    ;;
  *) echo "Unknown STAGE=${STAGE}" >&2; exit 2 ;;
esac
