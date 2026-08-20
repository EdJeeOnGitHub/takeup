#!/usr/bin/env bash
#SBATCH --account=pi-akaring
#SBATCH --partition=amd
#SBATCH --qos=amd
#SBATCH --job-name=struct-sm-robustness
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G
#SBATCH --time=0-10:00:00
#SBATCH --output=temp/log/struct-sm-robustness-%A_%a.log
#SBATCH --error=temp/log/struct-sm-robustness-%A_%a.log
#SBATCH --array=0-5

set -euo pipefail

INPUT_PATH=${INPUT_PATH:-/project/akaring/takeup-data/data/stan_analysis_data}
OUTPUT_PATH=${OUTPUT_PATH:-/project/akaring/takeup-data/temp-data}
DATA_PATH=${DATA_PATH:-data}

fits=(
  105
  105
  105
  105
  106
  106
)

models=(
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB"
)

required_quantities=(
  "cluster_social_multiplier"
  "cluster_mu_rep"
  "cluster_mu_rep_deriv"
  "cluster_delta"
  "cluster_delta_deriv"
  "dist_beta_v"
)

task_id=${SLURM_ARRAY_TASK_ID:-${1:-}}
if [[ -z "${task_id}" || ! "${task_id}" =~ ^[0-5]$ ]]; then
  echo "Usage: sbatch --array=0-5 $0 (or: $0 TASK_ID)" >&2
  exit 2
fi

fit=${fits[$task_id]}
model=${models[$task_id]}

echo "Checking fit ${fit}, model ${model}"
for chain in 1 2 3 4; do
  chain_file="${INPUT_PATH}/dist_fit${fit}_${model}-${chain}.csv"
  if [[ ! -f "${chain_file}" ]]; then
    echo "Missing chain file: ${chain_file}" >&2
    exit 3
  fi

  header_file=$(mktemp)
  awk '!/^#/ { print; exit }' "${chain_file}" |
    tr ',' '\n' |
    sed -E 's/\..*$//' |
    sort -u > "${header_file}"

  for quantity in "${required_quantities[@]}"; do
    if ! grep -Fxq "${quantity}" "${header_file}"; then
      echo "Chain ${chain} lacks generated quantity: ${quantity}" >&2
      rm -f "${header_file}"
      exit 4
    fi
  done
  rm -f "${header_file}"
done

mkdir -p "${OUTPUT_PATH}"
module load -f R/4.3.1
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/projects/aggregate-transfers/renv/library/R-4.3/x86_64-pc-linux-gnu}

Rscript --no-save --no-restore quick_roc_postprocess.R \
  "${fit}" \
  "--input-path=${INPUT_PATH}" \
  "--output-path=${OUTPUT_PATH}" \
  "--data-path=${DATA_PATH}" \
  "--model=${model}" \
  --sm \
  1 2 3 4
