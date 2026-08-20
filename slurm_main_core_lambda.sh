#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=core-lambda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --time=1-12:00:00
#SBATCH --output=temp/log/core-lambda-%A-%a.log
#SBATCH --error=temp/log/core-lambda-%A-%a.log

set -euo pipefail

STAGE=${STAGE:?Set STAGE to prepare, reuse_common, mode, sample, gq, or summarize}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
OUTPUT_PATH=${OUTPUT_PATH:-${ANALYSIS_ROOT}/main-core-lambda-identification}
MANIFEST=${MANIFEST:-${OUTPUT_PATH}/lambda-specification-manifest.csv}
STAN_PATH=${STAN_PATH:-stan_models_fit105}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.33.1}
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP}
ITER_WARMUP=${ITER_WARMUP:-400}
ITER_SAMPLING=${ITER_SAMPLING:-400}
BASELINE_CSV_DIR=${BASELINE_CSV_DIR:-${ANALYSIS_ROOT}/main-core-baseline-production}

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

manifest_value() {
  local row=$1
  local column=$2
  Rscript --no-save --no-restore --no-init-file -e \
    "x <- read.csv('${MANIFEST}', stringsAsFactors=FALSE); cat(x[[${column}]][${row}])"
}

case "${STAGE}" in
  prepare)
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-lambda-manifest.R \
      "--output-path=${OUTPUT_PATH}"
    ;;
  reuse_common)
    common_path="${OUTPUT_PATH}/fits/common"
    mkdir -p "${common_path}"
    index=0
    while IFS= read -r csv; do
      index=$((index + 1))
      cp "${csv}" "${common_path}/common-chain${index}-1.csv"
    done < <(find "${BASELINE_CSV_DIR}" -maxdepth 1 -type f \
      -name 'dist_fit106_MAIN_CORE_chain*-1.csv' | sort)
    if [[ "${index}" -ne 4 ]]; then
      echo "Expected four baseline chains; copied ${index}." >&2
      exit 2
    fi
    ;;
  mode)
    : "${SLURM_ARRAY_TASK_ID:?mode requires --array=1-9}"
    spec=${SLURM_ARRAY_TASK_ID}
    structure=$(manifest_value "${spec}" 2)
    prior=$(manifest_value "${spec}" 3)
    label=$(manifest_value "${spec}" 5)
    seed=$(manifest_value "${spec}" 6)
    Rscript --no-save --no-restore --no-init-file \
      scratch/fit-main-core-weighted-mode.R \
      "--workspace=${WORKSPACE}" "--model=${MODEL}" \
      "--output-path=${OUTPUT_PATH}/modes" "--label=${label}" \
      "--core-lambda-structure=${structure}" \
      "--core-lambda-log-ratio-sd-prior=${prior}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      "--seed=${seed}" --threads=8
    ;;
  sample)
    : "${SLURM_ARRAY_TASK_ID:?sample requires --array=1-36}"
    task=$((SLURM_ARRAY_TASK_ID - 1))
    spec=$((task / 4 + 1))
    chain=$((task % 4 + 1))
    structure=$(manifest_value "${spec}" 2)
    prior=$(manifest_value "${spec}" 3)
    label=$(manifest_value "${spec}" 5)
    seed=$(manifest_value "${spec}" 6)
    init_file="${OUTPUT_PATH}/modes/${label}/mode-init.json"
    output_dir="${OUTPUT_PATH}/fits/${label}"
    mkdir -p "${output_dir}"
    Rscript --no-save --no-restore --no-init-file \
      scripts/structural/sample-main-core.R \
      "--model=${MODEL}" "--input-path=$(dirname "${WORKSPACE}")" \
      "--output-path=${output_dir}" "--stan-path=${STAN_PATH}" \
      --stan-file=takeup_struct_main_core.stan \
      "--cmdstan-path=${CMDSTAN_PATH}" \
      "--output-basename=${label}-chain${chain}" \
      --chains=1 --parallel-chains=1 --threads-per-chain=8 \
      "--iter-warmup=${ITER_WARMUP}" \
      "--iter-sampling=${ITER_SAMPLING}" \
      --adapt-delta=0.999 --max-treedepth=12 --metric=diag_e \
      "--seed=$((seed + chain))" "--init-files=${init_file}" \
      "--core-lambda-structure=${structure}" \
      "--core-lambda-log-ratio-sd-prior=${prior}"
    ;;
  gq)
    : "${SLURM_ARRAY_TASK_ID:?gq requires --array=1-9}"
    spec=${SLURM_ARRAY_TASK_ID}
    structure=$(manifest_value "${spec}" 2)
    prior=$(manifest_value "${spec}" 3)
    label=$(manifest_value "${spec}" 5)
    fit_csvs=$(find "${OUTPUT_PATH}/fits/${label}" -maxdepth 1 \
      -type f -name "${label}-chain*-*.csv" | sort | paste -sd, -)
    Rscript --no-save --no-restore --no-init-file \
      scripts/structural/generate-compact-gq.R \
      "--workspace=${WORKSPACE}" "--model=${MODEL}" \
      "--fit-csvs=${fit_csvs}" \
      "--output-path=${OUTPUT_PATH}/gq/${label}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      "--core-lambda-structure=${structure}" \
      "--core-lambda-log-ratio-sd-prior=${prior}"
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scratch/summarize-main-core-lambda-grid.R \
      "--manifest=${MANIFEST}" "--output-path=${OUTPUT_PATH}"
    ;;
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
