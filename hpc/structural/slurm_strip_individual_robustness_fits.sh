#!/usr/bin/env bash
#SBATCH --partition=broadwl
#SBATCH --job-name=strip-indiv-fits
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --array=1-2
#SBATCH --output=temp/log/strip-indiv-fits-%A-%a.log
#SBATCH --error=temp/log/strip-indiv-fits-%A-%a.log

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}
ANALYSIS_ROOT=${ANALYSIS_ROOT:-/project/akaring/takeup-data/data/stan_analysis_data}
OUTPUT_ROOT=${OUTPUT_ROOT:-${ANALYSIS_ROOT}/streamlined-active-robustness}
SMOKE_ROOT=${SMOKE_ROOT:-${ANALYSIS_ROOT}/streamlined-active-robustness-smoke}
task=${SLURM_ARRAY_TASK_ID:?Submit as array 1-2}

case ${task} in
  1)
    spec_id=private-distance-community-image
    legacy_stem=dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS
    ;;
  2)
    spec_id=full-information
    legacy_stem=dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP
    ;;
  *) exit 2 ;;
esac

schema=${SMOKE_ROOT}/${spec_id}/fits/${spec_id}-slim-chain1-1.csv
sources=()
for chain in 1 2 3 4; do
  sources+=("${ANALYSIS_ROOT}/${legacy_stem}-${chain}.csv")
done

cd "${PROJECT_ROOT}"
mkdir -p temp/log "${OUTPUT_ROOT}/${spec_id}/fits"
python3 scripts/structural/strip-cmdstan-to-slim-schema.py \
  --schema="${schema}" \
  --output-dir="${OUTPUT_ROOT}/${spec_id}/fits" \
  --output-prefix="${spec_id}-slim" \
  "${sources[@]}"
