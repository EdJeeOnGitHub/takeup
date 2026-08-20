#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=lambda-recovery
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=temp/log/lambda-recovery-%A-%a.log
#SBATCH --error=temp/log/lambda-recovery-%A-%a.log

set -euo pipefail
STAGE=${STAGE:?Set STAGE to generate, fit, hmc, hmc_gq, summarize, or assemble}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
ROOT_OUTPUT=${ROOT_OUTPUT:-${ANALYSIS_ROOT}/main-core-lambda-identification}
OUTPUT_PATH=${OUTPUT_PATH:-${ROOT_OUTPUT}/recovery}
MANIFEST=${MANIFEST:-${OUTPUT_PATH}/lambda-recovery-manifest.csv}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP}
STAN_PATH=${STAN_PATH:-stan_models_fit105}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.33.1}
BASELINE_CSV_DIR=${BASELINE_CSV_DIR:-${ANALYSIS_ROOT}/main-core-baseline-production}
BASELINE_INIT=${BASELINE_INIT:-${ROOT_OUTPUT}/modes/common/mode-init.json}

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

manifest_value() {
  local row=$1
  local name=$2
  Rscript --no-save --no-restore --no-init-file -e \
    "x <- read.csv('${MANIFEST}', stringsAsFactors=FALSE); cat(x[['${name}']][${row}])"
}
audit_value() {
  local row=$1
  local name=$2
  Rscript --no-save --no-restore --no-init-file -e \
    "x <- subset(read.csv('${MANIFEST}', stringsAsFactors=FALSE), replicate <= 2); cat(x[['${name}']][${row}])"
}

case "${STAGE}" in
  generate)
    baseline_csvs=$(find "${BASELINE_CSV_DIR}" -maxdepth 1 -type f \
      -name 'dist_fit106_MAIN_CORE_chain*-1.csv' | sort | paste -sd, -)
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-lambda-recovery.R \
      "--workspace=${WORKSPACE}" "--model=${MODEL}" \
      "--fit-csvs=${baseline_csvs}" "--output-path=${OUTPUT_PATH}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --replicates=50 --seed=20263000
    ;;
  fit)
    : "${SLURM_ARRAY_TASK_ID:?fit requires --array=1-750}"
    task=${SLURM_ARRAY_TASK_ID}
    data_json=$(manifest_value "${task}" data_json)
    structure=$(manifest_value "${task}" lambda_structure)
    prior=$(manifest_value "${task}" prior_sd)
    label=$(manifest_value "${task}" label)
    seed=$(manifest_value "${task}" seed)
    Rscript --no-save --no-restore --no-init-file \
      scratch/fit-main-core-weighted-mode.R \
      "--data-json=${data_json}" "--label=${label}" \
      "--output-path=${OUTPUT_PATH}/fits" "--init-json=${BASELINE_INIT}" \
      "--core-lambda-structure=${structure}" \
      "--core-lambda-log-ratio-sd-prior=${prior}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      "--seed=${seed}" --threads=8
    ;;
  hmc)
    : "${SLURM_ARRAY_TASK_ID:?hmc requires --array=1-120}"
    task=$((SLURM_ARRAY_TASK_ID - 1))
    audit_row=$((task / 4 + 1))
    chain=$((task % 4 + 1))
    data_json=$(audit_value "${audit_row}" data_json)
    structure=$(audit_value "${audit_row}" lambda_structure)
    prior=$(audit_value "${audit_row}" prior_sd)
    label=$(audit_value "${audit_row}" label)
    seed=$(audit_value "${audit_row}" seed)
    init_file="${OUTPUT_PATH}/fits/${label}/mode-init.json"
    hmc_path="${OUTPUT_PATH}/hmc/${label}"
    mkdir -p "${hmc_path}"
    Rscript --no-save --no-restore --no-init-file \
      scratch/sample-slim-individual-fp.R \
      "--model=${MODEL}" "--input-path=$(dirname "${WORKSPACE}")" \
      "--data-json=${data_json}" "--output-path=${hmc_path}" \
      "--stan-path=${STAN_PATH}" --stan-file=takeup_struct_main_core.stan \
      "--cmdstan-path=${CMDSTAN_PATH}" \
      "--output-basename=${label}-chain${chain}" \
      --chains=1 --parallel-chains=1 --threads-per-chain=8 \
      --iter-warmup=200 --iter-sampling=200 --adapt-delta=0.999 \
      --max-treedepth=12 --metric=diag_e \
      "--seed=$((seed + chain))" "--init-files=${init_file}" \
      "--core-lambda-structure=${structure}" \
      "--core-lambda-log-ratio-sd-prior=${prior}"
    ;;
  hmc_gq)
    : "${SLURM_ARRAY_TASK_ID:?hmc_gq requires --array=1-30}"
    audit_row=${SLURM_ARRAY_TASK_ID}
    data_json=$(audit_value "${audit_row}" data_json)
    structure=$(audit_value "${audit_row}" lambda_structure)
    prior=$(audit_value "${audit_row}" prior_sd)
    label=$(audit_value "${audit_row}" label)
    fit_csvs=$(find "${OUTPUT_PATH}/hmc/${label}" -maxdepth 1 -type f \
      -name "${label}-chain*-*.csv" | sort | paste -sd, -)
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-compact-gq.R \
      "--data-json=${data_json}" "--fit-csvs=${fit_csvs}" \
      "--output-path=${OUTPUT_PATH}/hmc-gq/${label}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      "--core-lambda-structure=${structure}" \
      "--core-lambda-log-ratio-sd-prior=${prior}"
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scratch/summarize-main-core-lambda-recovery.R \
      "--manifest=${MANIFEST}" "--output-path=${OUTPUT_PATH}"
    ;;
  assemble)
    Rscript --no-save --no-restore --no-init-file \
      scratch/assemble-main-core-lambda-appendix.R \
      "--source-path=${ROOT_OUTPUT}" \
      "--output-path=${ROOT_OUTPUT}/paper-package"
    ;;
  *) echo "Unknown STAGE=${STAGE}" >&2; exit 2 ;;
esac
