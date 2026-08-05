#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=core-boot999
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=temp/log/core-boot999-%A-%a.log
#SBATCH --error=temp/log/core-boot999-%A-%a.log

set -euo pipefail
STAGE=${STAGE:?Set STAGE to weights, baseline, bootstrap, or summarize}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
WEIGHTS_PATH=${WEIGHTS_PATH:-${ANALYSIS_ROOT}/main-core-cluster-bootstrap-999-weights}
OUTPUT_PATH=${OUTPUT_PATH:-${ANALYSIS_ROOT}/main-core-cluster-bootstrap-999}
PAPER_OUTPUT=${PAPER_OUTPUT:-${ANALYSIS_ROOT}/main-core-cluster-bootstrap-999-results}
MANIFEST=${MANIFEST:-${WEIGHTS_PATH}/weight-manifest.csv}
BASELINE_INIT=${BASELINE_INIT:-${OUTPUT_PATH}/unweighted/mode-init.json}
STAN_PATH=${STAN_PATH:-stan_models}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.33.1}

mkdir -p temp/log "${WEIGHTS_PATH}" "${OUTPUT_PATH}" "${PAPER_OUTPUT}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

case "${STAGE}" in
  weights)
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-cluster-weights.R \
      "--workspace=${WORKSPACE}" "--output-path=${WEIGHTS_PATH}" \
      --methods=multinomial --replicates=999 --seed=20260802
    ;;
  baseline)
    Rscript --no-save --no-restore --no-init-file \
      scratch/fit-main-core-weighted-mode.R \
      "--workspace=${WORKSPACE}" "--output-path=${OUTPUT_PATH}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --threads=8
    ;;
  bootstrap)
    : "${SLURM_ARRAY_TASK_ID:?bootstrap stage requires an array task}"
    Rscript --no-save --no-restore --no-init-file \
      scratch/fit-main-core-weighted-mode.R \
      "--workspace=${WORKSPACE}" "--manifest=${MANIFEST}" \
      "--task-id=${SLURM_ARRAY_TASK_ID}" "--init-json=${BASELINE_INIT}" \
      "--output-path=${OUTPUT_PATH}" "--stan-path=${STAN_PATH}" \
      "--cmdstan-path=${CMDSTAN_PATH}" --threads=8
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scratch/summarize-main-core-cluster-robustness.R \
      "--weighted-path=${OUTPUT_PATH}" "--output-path=${PAPER_OUTPUT}" \
      "--table-path=${PAPER_OUTPUT}/main-core-cluster-bootstrap-multipliers.tex" \
      "--figure-path=${PAPER_OUTPUT}/main-core-cluster-bootstrap-multipliers.pdf"
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-cluster-bootstrap-ate-table.R \
      "--draws=${PAPER_OUTPUT}/estimand-draws.csv" \
      "--table-path=${PAPER_OUTPUT}/main-core-cluster-bootstrap-overall-te-table.tex" \
      "--summary-path=${PAPER_OUTPUT}/cluster-bootstrap-structural-results.csv"
    ;;
  *) echo "Unknown STAGE=${STAGE}" >&2; exit 2 ;;
esac
