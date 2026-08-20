#!/bin/bash
#SBATCH --partition=bigmem2
#SBATCH --job-name=robust-pp-106
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=0-10:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/robust-pp-106-combined-%j.log
#SBATCH --error=temp/log/robust-pp-106-combined-%j.log
#SBATCH --export=IN_SLURM=1

set -euo pipefail

cd /home/edjee/projects/takeup

module load R/4.2.0
source quick_postprocess.sh

DATA_PATH="/project/akaring/takeup-data/data/stan_analysis_data"

models=(
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB"
  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB_MISSING_MARGINALIZED"
)

for model in "${models[@]}"; do
  postprocess_struct_models "${model}" "106" \
    "--input-path=${DATA_PATH} --output-path=${DATA_PATH}"
done
