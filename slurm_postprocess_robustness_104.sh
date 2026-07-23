#!/bin/bash
#SBATCH --partition=bigmem2
#SBATCH --job-name=robust-pp-104
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=0-10:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=edjee96@gmail.com
#SBATCH --output=temp/log/robust-pp-104-%j.log
#SBATCH --error=temp/log/robust-pp-104-%j.log
#SBATCH --export=IN_SLURM=1

set -euo pipefail

VERSION=${1:-104}
INPUT_PATH="/project/akaring/takeup-data/data/stan_analysis_data"
OUTPUT_PATH="/project/akaring/takeup-data/temp-data"

if [[ -v IN_SLURM ]]; then
  module load -f R/4.2.0
fi

models=(
#  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS"
#  "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP"
   "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"
)

for model in "${models[@]}"; do
  Rscript --no-save \
          --no-restore \
          --verbose \
          quick_ate_postprocess.R \
          "${VERSION}" \
          --model="${model}" \
          --input-path="${INPUT_PATH}" \
          --output-path="${OUTPUT_PATH}" \
          1 2 3 4

  Rscript --no-save \
          --no-restore \
          --verbose \
          quick_roc_postprocess.R \
          "${VERSION}" \
          --model="${model}" \
          --sm \
          --input-path="${INPUT_PATH}" \
          --output-path="${OUTPUT_PATH}" \
          1 2 3 4
done
