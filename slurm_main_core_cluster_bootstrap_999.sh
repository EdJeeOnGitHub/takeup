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
STAGE=${STAGE:?Set STAGE to weights, baseline, bootstrap, launch, or summarize}
PROJECT_ROOT=${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}
WORKER_SCRIPT=${WORKER_SCRIPT:-${PROJECT_ROOT}/slurm_main_core_cluster_bootstrap_999.sh}
TOTAL_REPLICATES=${TOTAL_REPLICATES:-999}
SLURM_PARTITION_NAME=${SLURM_PARTITION_NAME:-${SLURM_JOB_PARTITION:-broadwl}}
SLURM_ACCOUNT_NAME=${SLURM_ACCOUNT_NAME:-}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
WEIGHT_METHOD=${WEIGHT_METHOD:-exponential}
WEIGHTS_PATH=${WEIGHTS_PATH:-${ANALYSIS_ROOT}/main-core-exponential-cluster-weight-999-weights}
OUTPUT_PATH=${OUTPUT_PATH:-${ANALYSIS_ROOT}/main-core-exponential-cluster-weight-999}
PAPER_OUTPUT=${PAPER_OUTPUT:-${ANALYSIS_ROOT}/main-core-exponential-cluster-weight-999-results}
MANIFEST=${MANIFEST:-${WEIGHTS_PATH}/weight-manifest.csv}
BASELINE_INIT=${BASELINE_INIT:-${OUTPUT_PATH}/unweighted/mode-init.json}
STAN_PATH=${STAN_PATH:-stan_models}
if [[ -z "${CMDSTAN_PATH:-}" ]]; then
  if [[ -d /home/edjee/.cmdstan/cmdstan-2.35.0 ]]; then
    CMDSTAN_PATH=/home/edjee/.cmdstan/cmdstan-2.35.0
  else
    CMDSTAN_PATH=/home/edjee/.cmdstan/cmdstan-2.33.1
  fi
fi

cd "${PROJECT_ROOT}"

mkdir -p temp/log "${WEIGHTS_PATH}" "${OUTPUT_PATH}" "${PAPER_OUTPUT}"
module load -f R/4.2.0
if [[ -z "${R_LIBS_USER:-}" ]]; then
  if [[ -d /home/edjee/R/x86_64-pc-linux-gnu-library/4.2 ]]; then
    R_LIBS_USER=/home/edjee/R/x86_64-pc-linux-gnu-library/4.2
  else
    R_LIBS_USER=/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu
  fi
fi
export R_LIBS_USER
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

sbatch_self() {
  local site_args=(--partition="${SLURM_PARTITION_NAME}")
  if [[ -n "${SLURM_ACCOUNT_NAME}" ]]; then
    site_args+=(--account="${SLURM_ACCOUNT_NAME}")
  fi
  sbatch --parsable "${site_args[@]}" "$@" "${WORKER_SCRIPT}"
}

case "${STAGE}" in
  weights)
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-cluster-weights.R \
      "--workspace=${WORKSPACE}" "--output-path=${WEIGHTS_PATH}" \
      "--methods=${WEIGHT_METHOD}" "--replicates=${TOTAL_REPLICATES}" \
      --seed=20260802
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
  launch)
    : "${NEXT_START:?launch stage requires NEXT_START}"
    batch_size=${BATCH_SIZE:-333}
    max_concurrent=${MAX_CONCURRENT:-32}
    total_replicates=${TOTAL_REPLICATES}
    if (( NEXT_START > total_replicates )); then
      summary_job=$(sbatch_self --export=ALL,STAGE=summarize)
      printf 'summary %s\n' "${summary_job}"
      exit 0
    fi
    batch_end=$((NEXT_START + batch_size - 1))
    if (( batch_end > total_replicates )); then batch_end=${total_replicates}; fi
    bootstrap_job=$(sbatch_self \
      --array="${NEXT_START}-${batch_end}%${max_concurrent}" \
      --export=ALL,STAGE=bootstrap)
    if (( batch_end < total_replicates )); then
      next_job=$(sbatch_self --dependency="afterok:${bootstrap_job}" \
        --export=ALL,STAGE=launch,NEXT_START=$((batch_end + 1)),BATCH_SIZE="${batch_size}",MAX_CONCURRENT="${max_concurrent}" \
      )
      printf 'bootstrap_%s_%s %s\nlauncher %s\n' \
        "${NEXT_START}" "${batch_end}" "${bootstrap_job}" "${next_job}"
    else
      summary_job=$(sbatch_self --dependency="afterok:${bootstrap_job}" \
        --export=ALL,STAGE=summarize)
      printf 'bootstrap_%s_%s %s\nsummary %s\n' \
        "${NEXT_START}" "${batch_end}" "${bootstrap_job}" "${summary_job}"
    fi
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scratch/summarize-main-core-cluster-robustness.R \
      "--weighted-path=${OUTPUT_PATH}" "--output-path=${PAPER_OUTPUT}" \
      "--table-path=${PAPER_OUTPUT}/main-core-exponential-cluster-weight-multipliers.tex" \
      "--figure-path=${PAPER_OUTPUT}/main-core-exponential-cluster-weight-multipliers.pdf"
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-cluster-bootstrap-ate-table.R \
      "--draws=${PAPER_OUTPUT}/estimand-draws.csv" \
      "--table-path=${PAPER_OUTPUT}/main-core-exponential-cluster-weight-overall-te-table.tex" \
      "--summary-path=${PAPER_OUTPUT}/exponential-cluster-weight-structural-results.csv" \
      "--method=${WEIGHT_METHOD}"
    ;;
  *) echo "Unknown STAGE=${STAGE}" >&2; exit 2 ;;
esac
