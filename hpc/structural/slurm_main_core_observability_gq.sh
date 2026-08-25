#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=obs-gq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=temp/log/obs-gq-%A_%a.log
#SBATCH --error=temp/log/obs-gq-%A_%a.log

set -euo pipefail

CHAIN_ID=${CHAIN_ID:-${SLURM_ARRAY_TASK_ID:-}}
SPECIFICATION=${SPECIFICATION:?Set an observability specification}
FIT_ROOT=${FIT_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-observability-ladder}
INPUT_PATH=${INPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data/main-core-asym-input}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.33.1}
DISTANCE_DEFINITION=${DISTANCE_DEFINITION:-assigned}
SM_EVALUATION_DISTANCE_M=${SM_EVALUATION_DISTANCE_M:-}

case "${SPECIFICATION}" in
  f1) observation=1; recognition=0; report=1; hierarchical=0; prior_scale=0.25 ;;
  f2) observation=1; recognition=2; report=1; hierarchical=0; prior_scale=0.25 ;;
  f3) observation=1; recognition=2; report=2; hierarchical=0; prior_scale=0.25 ;;
  u3) observation=2; recognition=1; report=2; hierarchical=0; prior_scale=0.25 ;;
  hierarchical) observation=1; recognition=0; report=0; hierarchical=1; prior_scale=0.25 ;;
  tight) observation=1; recognition=0; report=0; hierarchical=0; prior_scale=0.10 ;;
  design-pooled) observation=1; recognition=0; report=0; hierarchical=2; prior_scale=0.25 ;;
  *) echo "Unknown SPECIFICATION=${SPECIFICATION}" >&2; exit 2 ;;
esac

fit_csv=${FIT_ROOT}/${SPECIFICATION}/dist_fit106_MAIN_CORE_chain${CHAIN_ID}-1.csv
output_path=${GQ_OUTPUT_PATH:-${FIT_ROOT}/${SPECIFICATION}/${GQ_SUBDIR:-gq}}
mkdir -p temp/log "${output_path}"
module load gcc/10.2.0
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/takeup/renv/library/R-4.2/x86_64-pc-linux-gnu}
export CMDSTAN=${CMDSTAN_PATH}
export CMDSTANR_NO_VER_CHECK=TRUE

extra=()
if [[ -n ${SM_EVALUATION_DISTANCE_M} ]]; then
  extra+=("--sm-evaluation-distance-m=${SM_EVALUATION_DISTANCE_M}")
fi

Rscript --no-save --no-restore --no-init-file \
  scripts/structural/generate-compact-gq.R \
  "--workspace=${INPUT_PATH}/dist_fit105.RData" \
  --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS \
  --stan-path=stan_models \
  "--fit-csvs=${fit_csv}" \
  "--output-path=${output_path}" \
  "--output-basename=${SPECIFICATION}-compact-chain${CHAIN_ID}" \
  "--core-observation-model=${observation}" \
  "--core-recognition-structure=${recognition}" \
  "--core-report-structure=${report}" \
  "--core-report-arm-dist-hierarchical=${hierarchical}" \
  "--core-report-arm-dist-prior-scale=${prior_scale}" \
  "--distance-definition=${DISTANCE_DEFINITION}" \
  "--force-recompile=${FORCE_RECOMPILE:-0}" \
  --threads=4 --parallel-chains=1 \
  "--cmdstan-path=${CMDSTAN_PATH}" \
  "${extra[@]}"
