#!/usr/bin/env bash
#SBATCH --partition=caslake
#SBATCH --account=pi-akaring
#SBATCH --job-name=slim-active-gq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --array=1-4
#SBATCH --output=temp/log/slim-active-gq-%A-%a.log
#SBATCH --error=temp/log/slim-active-gq-%A-%a.log

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
PREPARED_ROOT=${PREPARED_ROOT:-/project/akaring/takeup-data/candidate-multiplier-1250-20260825/workspaces}
OUTPUT_ROOT=${OUTPUT_ROOT:-${ANALYSIS_ROOT}/streamlined-active-robustness}
task=${SLURM_ARRAY_TASK_ID:?Submit as an array from 1 to 4}

case ${task} in
  1)
    spec_id=private-distance-community-image
    model_key=indiv_dist_community_fp_indiv_vis
    input_root=${ANALYSIS_ROOT}
    ;;
  2)
    spec_id=full-information
    model_key=indiv_dist_indiv_fp
    input_root=${ANALYSIS_ROOT}
    ;;
  3)
    spec_id=exclude-dispersed
    model_key=no_outliers
    input_root=${PREPARED_ROOT}
    ;;
  4)
    spec_id=second-order-observability
    model_key=second_order_observability
    input_root=${PREPARED_ROOT}
    ;;
  *) exit 2 ;;
esac

fit_csvs=()
for chain in 1 2 3 4; do
  fit_csvs+=("${OUTPUT_ROOT}/${spec_id}/fits/${spec_id}-slim-chain${chain}-1.csv")
done
fit_csv_arg=$(IFS=,; echo "${fit_csvs[*]}")
gq_path=${OUTPUT_ROOT}/${spec_id}/gq-1250

cd "${PROJECT_ROOT}"
mkdir -p temp/log "${gq_path}"
module load gcc/10.2.0
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/R/x86_64-pc-linux-gnu-library/4.2}
export CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.35.0}
export CMDSTANR_NO_VER_CHECK=TRUE

Rscript --no-save --no-restore --no-init-file \
  scripts/appendix/generate-compact-individual-sm-gq.R \
  --model="${model_key}" \
  --input-path="${input_root}" \
  --output-path="${gq_path}" \
  --fit-csvs="${fit_csv_arg}" \
  --parallel-chains=4 --threads-per-chain=1 \
  --sm-evaluation-distance-m=1250
