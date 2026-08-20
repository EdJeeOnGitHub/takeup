#!/usr/bin/env bash
#SBATCH --account=pi-akaring
#SBATCH --partition=caslake
#SBATCH --job-name=core-fmix-800
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=temp/log/core-fmix-800-%A-%a.log
#SBATCH --error=temp/log/core-fmix-800-%A-%a.log

set -euo pipefail

STAGE=${STAGE:?Set STAGE to compile, sample, gq, or summarize}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
WORKSPACE=${WORKSPACE:-${ANALYSIS_ROOT}/main-core-bootstrap-input/dist_fit105.RData}
OUTPUT_PATH=${OUTPUT_PATH:-${ANALYSIS_ROOT}/main-core-finite-mixture-robustness-800}
MODE_INIT=${MODE_INIT:-${ANALYSIS_ROOT}/main-core-finite-mixture-robustness/mode/finite-mixture/mode-init.json}
BASELINE_GQ_DIR=${BASELINE_GQ_DIR:-${ANALYSIS_ROOT}/main-core-lambda-identification/gq/common}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.35.0}
STAN_PATH=${STAN_PATH:-stan_models}
MODEL=${MODEL:-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP}
ITER_WARMUP=${ITER_WARMUP:-800}
ITER_SAMPLING=${ITER_SAMPLING:-800}

mkdir -p temp/log "${OUTPUT_PATH}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/R/x86_64-pc-linux-gnu-library/4.2}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE
export MAKEFLAGS=${MAKEFLAGS:--j8}
Rscript --no-save --no-restore --no-init-file -e '
  if (packageVersion("cmdstanr") < package_version("0.9.0")) {
    stop("The Midway3/CmdStan 2.35 workflow requires cmdstanr >= 0.9.0.")
  }
'

case "${STAGE}" in
  compile)
    Rscript --no-save --no-restore --no-init-file -e '
      suppressPackageStartupMessages(library(cmdstanr))
      set_cmdstan_path(Sys.getenv("CMDSTAN"))
      for (stan_file in c(
        "stan_models/takeup_struct_main_core.stan",
        "stan_models/takeup_struct_main_core_compact_gq.stan"
      )) {
        message("Compiling ", stan_file)
        executable <- sub("[.]stan$", "", stan_file)
        tryCatch(
          cmdstan_model(
            stan_file,
            include_paths = "stan_models",
            cpp_options = list(stan_threads = TRUE),
            force_recompile = !file.exists(executable)
          ),
          error = function(error) {
            # CmdStanR 0.6.1 can fail while reopening an intermediate file
            # after CmdStan 2.35 has already produced a valid executable.
            # Accept that narrow wrapper failure only when the executable is
            # present; the shell validates both binaries below.
            if (!file.exists(executable)) stop(error)
            message(
              "CmdStanR post-compile wrapper error; retaining executable: ",
              conditionMessage(error)
            )
          }
        )
      }
    '
    for executable in \
      stan_models/takeup_struct_main_core \
      stan_models/takeup_struct_main_core_compact_gq; do
      test -x "${executable}"
      "${executable}" info >/dev/null
    done
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
      "--seed=$((20261020 + chain))" \
      "--init-files=${MODE_INIT}" \
      --core-type-distribution=2
    ;;
  gq)
    fit_csvs=$(find "${OUTPUT_PATH}/fits/finite-mixture" -maxdepth 1 -type f \
      -name 'finite-mixture-chain*-*.csv' ! -name '*profile*' | sort | paste -sd, -)
    if [[ $(tr ',' '\n' <<<"${fit_csvs}" | grep -c .) -ne 4 ]]; then
      echo "Expected exactly four fitted chain CSVs; got: ${fit_csvs}" >&2
      exit 1
    fi
    Rscript --no-save --no-restore --no-init-file \
      scratch/generate-main-core-compact-gq.R \
      "--workspace=${WORKSPACE}" "--model=${MODEL}" \
      "--fit-csvs=${fit_csvs}" \
      "--output-path=${OUTPUT_PATH}/gq/finite-mixture" \
      "--stan-path=${STAN_PATH}" "--cmdstan-path=${CMDSTAN_PATH}" \
      --core-type-distribution=2 --force-recompile=0 \
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
  *)
    echo "Unknown STAGE=${STAGE}" >&2
    exit 2
    ;;
esac
