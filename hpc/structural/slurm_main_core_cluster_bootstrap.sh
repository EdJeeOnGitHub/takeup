#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=core-bootstrap
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=temp/log/core-bootstrap-%A-%a.log
#SBATCH --error=temp/log/core-bootstrap-%A-%a.log

set -euo pipefail

STAGE=${STAGE:?Set STAGE to prepare, baseline, bootstrap, or summarize}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
WORKSPACE_DIR=${WORKSPACE_DIR:-${ANALYSIS_ROOT}/main-core-bootstrap-input}
OUTPUT_PATH=${OUTPUT_PATH:-${ANALYSIS_ROOT}/main-core-exponential-cluster-weights}
STAN_PATH=${STAN_PATH:-stan_models}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.33.1}
MANIFEST=${MANIFEST:-${WORKSPACE_DIR}/cluster-weights/weight-manifest.csv}
BASELINE_INIT=${BASELINE_INIT:-${OUTPUT_PATH}/unweighted/mode-init.json}
PAPER_OUTPUT=${PAPER_OUTPUT:-${ANALYSIS_ROOT}/main-core-exponential-cluster-weight-results}

mkdir -p temp/log "${WORKSPACE_DIR}" "${OUTPUT_PATH}" "${PAPER_OUTPUT}"
if [[ "${STAGE}" == "prepare" ]]; then
  module load -f gdal/2.4.1 udunits/2.2 proj/6.1 cmake R/4.2.0
else
  module load -f R/4.2.0
fi
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

case "${STAGE}" in
  prepare)
    Rscript --no-save --no-restore --no-init-file scripts/structural/run-model.R takeup fit \
      --data-only --cmdstanr \
      --models=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP \
      --outputname=dist_fit105 "--output-path=${WORKSPACE_DIR}" \
      --num-mix-groups=1 --threads=4 --chains=1 --iter=2
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/generate-main-core-cluster-weights.R \
      "--workspace=${WORKSPACE_DIR}/dist_fit105.RData" \
      "--output-path=${WORKSPACE_DIR}/cluster-weights" \
      --methods=exponential --replicates=200 --seed=20260802
    ;;
  baseline)
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/fit-main-core-weighted-mode.R \
      "--workspace=${WORKSPACE_DIR}/dist_fit105.RData" \
      "--output-path=${OUTPUT_PATH}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --threads=8
    ;;
  bootstrap)
    : "${SLURM_ARRAY_TASK_ID:?bootstrap stage requires a Slurm array}"
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/fit-main-core-weighted-mode.R \
      "--workspace=${WORKSPACE_DIR}/dist_fit105.RData" \
      "--manifest=${MANIFEST}" "--task-id=${SLURM_ARRAY_TASK_ID}" \
      "--init-json=${BASELINE_INIT}" "--output-path=${OUTPUT_PATH}" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --threads=8
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/summarize-main-core-cluster-robustness.R \
      "--weighted-path=${OUTPUT_PATH}" "--output-path=${PAPER_OUTPUT}" \
      "--table-path=${PAPER_OUTPUT}/main-core-exponential-cluster-weight-multipliers.tex" \
      "--figure-path=${PAPER_OUTPUT}/main-core-exponential-cluster-weight-multipliers.pdf"
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/generate-main-core-cluster-bootstrap-ate-table.R \
      "--draws=${PAPER_OUTPUT}/estimand-draws.csv" \
      "--table-path=${PAPER_OUTPUT}/main-core-exponential-cluster-weight-overall-te-table.tex" \
      "--summary-path=${PAPER_OUTPUT}/exponential-cluster-weight-structural-results.csv" \
      --method=exponential
    ;;
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
