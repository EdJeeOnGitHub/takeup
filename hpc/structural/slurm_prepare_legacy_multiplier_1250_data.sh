#!/usr/bin/env bash
#SBATCH --partition=caslake
#SBATCH --account=pi-akaring
#SBATCH --job-name=sm-1250-data
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --array=1-3
#SBATCH --output=temp/log/sm-1250-data-%A-%a.log
#SBATCH --error=temp/log/sm-1250-data-%A-%a.log

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}
OUTPUT_ROOT=${OUTPUT_ROOT:?Set OUTPUT_ROOT to an isolated candidate directory}
task=${SLURM_ARRAY_TASK_ID:?This launcher requires an array task from 1 to 3}

case ${task} in
  1)
    model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS
    output_name=dist_fit1250
    belief_args=()
    ;;
  2)
    model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS
    output_name=dist_fit1251
    belief_args=(--beliefs-outcome=correct-observability)
    ;;
  3)
    model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB
    output_name=dist_fit1252
    belief_args=(--beliefs-outcome=sob)
    ;;
  *) exit 2 ;;
esac

cd "${PROJECT_ROOT}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/R/x86_64-pc-linux-gnu-library/4.2}
export DISTANCE_DEFINITION=assigned
mkdir -p temp/log "${OUTPUT_ROOT}/workspaces"
Rscript --no-save --no-restore scripts/structural/run-model.R takeup fit \
  --cmdstanr --data-only --models="${model}" --num-mix-groups=1 \
  --outputname="${output_name}" --output-path="${OUTPUT_ROOT}/workspaces" \
  "${belief_args[@]}"
