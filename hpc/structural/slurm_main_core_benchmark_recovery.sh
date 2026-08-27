#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=benchmark-recovery
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=temp/log/benchmark-recovery-%A-%a.log
#SBATCH --error=temp/log/benchmark-recovery-%A-%a.log

set -euo pipefail

STAGE=${STAGE:?Set STAGE to generate, compile, fit, or summarize}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
OUTPUT_PATH=${OUTPUT_PATH:-${ANALYSIS_ROOT}/main-core-benchmark-recovery}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
BASELINE_CSV_DIR=${BASELINE_CSV_DIR:-}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.35.0}
MANIFEST=${MANIFEST:-${OUTPUT_PATH}/benchmark-recovery-manifest.csv}

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

case "${STAGE}" in
  generate)
    : "${BASELINE_CSV_DIR:?Set BASELINE_CSV_DIR to the latest assigned-distance benchmark chain directory}"
    baseline_csvs=$(find "${BASELINE_CSV_DIR}" -maxdepth 1 -type f \
      -name '*.csv' ! -name '*diagnostics*' ! -name '*manifest*' | sort | paste -sd, -)
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/generate-main-core-benchmark-recovery.R \
      "--workspace=${WORKSPACE}" "--fit-csvs=${baseline_csvs}" \
      "--output-path=${OUTPUT_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --replicates=50 --grid-points=21 --seed=20260827
    ;;
  compile)
    Rscript --no-save --no-restore --no-init-file -e \
      "cmdstanr::set_cmdstan_path('${CMDSTAN_PATH}'); cmdstanr::cmdstan_model('stan_models/takeup_struct_main_core.stan', include_paths='stan_models', cpp_options=list(stan_threads=TRUE)); cmdstanr::cmdstan_model('stan_models/takeup_struct_main_core_compact_gq.stan', include_paths='stan_models', cpp_options=list(stan_threads=TRUE))"
    ;;
  fit)
    : "${SLURM_ARRAY_TASK_ID:?fit requires --array=1-50}"
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/fit-main-core-benchmark-recovery.R \
      "--manifest=${MANIFEST}" "--task-id=${SLURM_ARRAY_TASK_ID}" \
      "--output-path=${OUTPUT_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --threads=8 --chains=1 --iter-warmup=1000 --iter-sampling=1000 \
      --rerun-warmup=2000 --rerun-sampling=2000
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/summarize-main-core-benchmark-recovery.R \
      "--input-path=${OUTPUT_PATH}" --expected-replicates=50 \
      "--table-path=${OUTPUT_PATH}/paper/main-core-benchmark-recovery.tex" \
      "--figure-path=${OUTPUT_PATH}/paper/main-core-benchmark-likelihood.pdf"
    ;;
  *) echo "Unknown STAGE=${STAGE}" >&2; exit 2 ;;
esac
