#!/usr/bin/env bash
#SBATCH --account=pi-akaring
#SBATCH --partition=caslake
#SBATCH --job-name=core-priors
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1-12:00:00
#SBATCH --output=temp/log/core-priors-%A-%a.log
#SBATCH --error=temp/log/core-priors-%A-%a.log

set -euo pipefail

STAGE=${STAGE:?Set STAGE to prepare, compile, sample, gq, or summarize}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
OUTPUT_PATH=${OUTPUT_PATH:-${ANALYSIS_ROOT}/main-core-prior-grid}
MANIFEST=${MANIFEST:-${OUTPUT_PATH}/prior-grid-manifest.csv}
RERUN_AUDIT=${RERUN_AUDIT:-${OUTPUT_PATH}/prior-grid-needs-rerun.csv}
STAN_PATH=${STAN_PATH:-stan_models}
CMDSTAN_PATH=${CMDSTAN_PATH:-${CMDSTAN:-$HOME/.cmdstan/cmdstan-2.35.0}}
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP}
MODE_INIT=${MODE_INIT:-${ANALYSIS_ROOT}/main-core-lambda-identification/modes/common/mode-init.json}
ITER_WARMUP=${ITER_WARMUP:-400}
ITER_SAMPLING=${ITER_SAMPLING:-400}
MAX_TREEDEPTH=${MAX_TREEDEPTH:-12}
RERUN_FIT_ROOT=${RERUN_FIT_ROOT:-fits-rerun}
RERUN_GQ_ROOT=${RERUN_GQ_ROOT:-gq-rerun}

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/R/x86_64-pc-linux-gnu-library/4.2}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE
export MAKEFLAGS=${MAKEFLAGS:--j8}

manifest_value() {
  local row=$1
  local column=$2
  Rscript --no-save --no-restore --no-init-file -e \
    "x <- read.csv('${MANIFEST}', stringsAsFactors=FALSE, check.names=FALSE); cat(x[['${column}']][${row}])"
}

prior_args() {
  local spec=$1
  local columns=(
    mu_rep_sd dist_beta_v_sd raw_u_sd_alpha raw_u_sd_beta
    visibility_prior_multiplier beta_intercept_sd beta_ink_effect_sd
    beta_calendar_effect_sd beta_bracelet_effect_sd
    wtp_value_utility_mean wtp_value_utility_sd
    lnorm_wtp_value_utility_prior wtp_mu_prior_sd
  )
  local options=(
    mu-rep-sd dist-beta-v-sd raw-u-sd-alpha raw-u-sd-beta
    core-visibility-prior-multiplier beta-intercept-sd beta-ink-effect-sd
    beta-calendar-effect-sd beta-bracelet-effect-sd
    wtp-value-utility-mean wtp-value-utility-sd
    lnorm-wtp-value-utility-prior core-wtp-mu-prior-sd
  )
  local index
  for index in "${!columns[@]}"; do
    printf -- '--%s=%s\n' "${options[$index]}" \
      "$(manifest_value "${spec}" "${columns[$index]}")"
  done
}

rerun_spec() {
  local row=$1
  Rscript --no-save --no-restore --no-init-file -e \
    "x <- read.csv('${RERUN_AUDIT}', stringsAsFactors=FALSE); cat(x[['spec_id']][${row}])"
}

case "${STAGE}" in
  prepare)
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/generate-main-core-prior-grid-manifest.R \
      "--output-path=${OUTPUT_PATH}"
    ;;
  compile)
    Rscript --no-save --no-restore --no-init-file -e \
      "options(cmdstanr_no_ver_check=TRUE); Sys.setenv(CMDSTAN='${CMDSTAN_PATH}', CMDSTANR_NO_VER_CHECK='TRUE'); library(cmdstanr); set_cmdstan_path('${CMDSTAN_PATH}'); cmdstan_model('${STAN_PATH}/takeup_struct_main_core.stan', include_paths='${STAN_PATH}', cpp_options=list(stan_threads=TRUE)); cmdstan_model('${STAN_PATH}/takeup_struct_main_core_compact_gq.stan', include_paths='${STAN_PATH}', cpp_options=list(stan_threads=TRUE))"
    ;;
  sample|sample-rerun)
    : "${SLURM_ARRAY_TASK_ID:?sample stage requires an array task}"
    task=$((SLURM_ARRAY_TASK_ID - 1))
    if [[ "${STAGE}" == sample-rerun ]]; then
      rerun_row=$((task / 4 + 1))
      spec=$(rerun_spec "${rerun_row}")
      fit_root=${RERUN_FIT_ROOT}
    else
      spec=$((task / 4 + 1))
      fit_root=fits
    fi
    chain=$((task % 4 + 1))
    label=$(manifest_value "${spec}" label)
    seed=$(manifest_value "${spec}" seed)
    output_dir="${OUTPUT_PATH}/${fit_root}/${label}"
    mkdir -p "${output_dir}"
    mapfile -t priors < <(prior_args "${spec}")
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
      --adapt-delta=0.999 "--max-treedepth=${MAX_TREEDEPTH}" --metric=diag_e \
      "--seed=$((seed + chain))" "--init-files=${MODE_INIT}" \
      "${priors[@]}"
    ;;
  gq|gq-rerun)
    : "${SLURM_ARRAY_TASK_ID:?gq stage requires an array task}"
    if [[ "${STAGE}" == gq-rerun ]]; then
      spec=$(rerun_spec "${SLURM_ARRAY_TASK_ID}")
      fit_root=${RERUN_FIT_ROOT}
      gq_root=${RERUN_GQ_ROOT}
    else
      spec=${SLURM_ARRAY_TASK_ID}
      fit_root=fits
      gq_root=gq
    fi
    label=$(manifest_value "${spec}" label)
    fit_csvs=$(find "${OUTPUT_PATH}/${fit_root}/${label}" -maxdepth 1 -type f \
      -name "${label}-chain*-*.csv" | sort | paste -sd, -)
    if [[ $(tr ',' '\n' <<< "${fit_csvs}" | grep -c .) -ne 4 ]]; then
      echo "Expected four fit CSVs for ${label}." >&2
      exit 2
    fi
    mapfile -t priors < <(prior_args "${spec}")
    Rscript --no-save --no-restore --no-init-file \
      scripts/structural/generate-compact-gq.R \
      "--workspace=${WORKSPACE}" "--model=${MODEL}" \
      "--fit-csvs=${fit_csvs}" \
      "--output-path=${OUTPUT_PATH}/${gq_root}/${label}" \
      "--output-basename=${label}-compact" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      "${priors[@]}"
    ;;
  summarize)
    Rscript --no-save --no-restore --no-init-file \
      scripts/appendix/summarize-main-core-prior-grid.R \
      "--manifest=${MANIFEST}" "--output-path=${OUTPUT_PATH}"
    ;;
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
