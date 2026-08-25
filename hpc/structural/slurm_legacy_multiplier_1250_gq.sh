#!/usr/bin/env bash
#SBATCH --partition=caslake
#SBATCH --account=pi-akaring
#SBATCH --job-name=sm-1250-legacy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --array=1-5
#SBATCH --output=temp/log/sm-1250-legacy-%A-%a.log
#SBATCH --error=temp/log/sm-1250-legacy-%A-%a.log

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}
INPUT_ROOT=${INPUT_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
OUTPUT_ROOT=${OUTPUT_ROOT:?Set OUTPUT_ROOT to an isolated candidate directory}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.35.0}
task=${SLURM_ARRAY_TASK_ID:?This launcher requires an array task from 1 to 2}
workspace_root=${WORKSPACE_ROOT:-${OUTPUT_ROOT}/workspaces}

case ${task} in
  1)
    model=indiv_dist_community_fp_indiv_vis
    label=private-distance-community-image
    ;;
  2)
    model=indiv_dist_indiv_fp
    label=full-information
    ;;
  3)
    model=no_outliers
    label=exclude-dispersed
    INPUT_ROOT=${workspace_root}
    fit_stem=dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS
    ;;
  4)
    model=correct_observability
    label=correct-classification
    INPUT_ROOT=${workspace_root}
    fit_stem=dist_fit106_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS
    ;;
  5)
    model=second_order_observability
    label=second-order-observability
    INPUT_ROOT=${workspace_root}
    fit_stem=dist_fit106_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB
    ;;
  *)
    echo "Unsupported array task: ${task}" >&2
    exit 2
    ;;
esac

cd "${PROJECT_ROOT}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/R/x86_64-pc-linux-gnu-library/4.2}
export CMDSTAN_PATH
mkdir -p temp/log "${OUTPUT_ROOT}/${label}"
extra=()
if [[ -n ${fit_stem:-} ]]; then
  fit_csvs=()
  for chain in 1 2 3 4; do
    fit_csvs+=("/project/akaring/takeup-data/data/stan_analysis_data/${fit_stem}-${chain}.csv")
  done
  fit_csv_arg=$(IFS=,; echo "${fit_csvs[*]}")
  extra+=(--fit-csvs="${fit_csv_arg}")
fi
Rscript --no-save --no-restore \
  scripts/appendix/generate-compact-individual-sm-gq.R \
  --model="${model}" \
  --input-path="${INPUT_ROOT}" \
  --output-path="${OUTPUT_ROOT}/${label}" \
  --parallel-chains=4 --threads-per-chain=1 \
  --sm-evaluation-distance-m=1250 \
  "${extra[@]}"
