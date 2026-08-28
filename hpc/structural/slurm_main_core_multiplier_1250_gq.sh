#!/usr/bin/env bash
#SBATCH --partition=caslake
#SBATCH --account=pi-akaring
#SBATCH --job-name=sm-1250-gq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-17
#SBATCH --output=temp/log/sm-1250-gq-%A-%a.log
#SBATCH --error=temp/log/sm-1250-gq-%A-%a.log

set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}
FIT_ROOT=${FIT_ROOT:-/project/akaring/takeup-data/candidate-hpc-cd5f295-assigned/work}
OUTPUT_ROOT=${OUTPUT_ROOT:?Set OUTPUT_ROOT to an isolated candidate directory}
CMDSTAN_PATH=${CMDSTAN_PATH:-/home/edjee/.cmdstan/cmdstan-2.35.0}
WORKSPACE=${WORKSPACE:-${PROJECT_ROOT}/build/structural-workspace/main-core-input.RData}
task=${SLURM_ARRAY_TASK_ID:?This launcher requires an array task from 1 to 17}

prior_labels=(
  baseline image-tight image-diffuse distance-tight distance-diffuse
  heterogeneity-tight heterogeneity-diffuse visibility-tight visibility-diffuse
  private-tight private-diffuse joint-tight joint-diffuse
)

label=
fit_dir=
extra=()
if (( task <= 13 )); then
  label=${prior_labels[$((task - 1))]}
  fit_dir=${FIT_ROOT}/prior-grid/fits/${label}
  case ${label} in
    image-tight) extra+=(--mu-rep-sd=0.5) ;;
    image-diffuse) extra+=(--mu-rep-sd=2) ;;
    distance-tight) extra+=(--dist-beta-v-sd=0.125) ;;
    distance-diffuse) extra+=(--dist-beta-v-sd=0.5) ;;
    heterogeneity-tight) extra+=(--raw-u-sd-alpha=8 --raw-u-sd-beta=3.34782608695652) ;;
    heterogeneity-diffuse) extra+=(--raw-u-sd-alpha=2.5 --raw-u-sd-beta=0.717391304347826) ;;
    visibility-tight) extra+=(--core-visibility-prior-multiplier=0.5) ;;
    visibility-diffuse) extra+=(--core-visibility-prior-multiplier=2) ;;
    private-tight)
      extra+=(--beta-intercept-sd=0.5 --beta-ink-effect-sd=0.125)
      extra+=(--beta-calendar-effect-sd=0.125 --beta-bracelet-effect-sd=0.125)
      extra+=(--wtp-value-utility-sd=0.00005 --core-wtp-mu-prior-sd=1)
      ;;
    private-diffuse)
      extra+=(--beta-intercept-sd=2 --beta-ink-effect-sd=0.5)
      extra+=(--beta-calendar-effect-sd=0.5 --beta-bracelet-effect-sd=0.5)
      extra+=(--wtp-value-utility-mean=-10 --wtp-value-utility-sd=4)
      extra+=(--lnorm-wtp-value-utility-prior=1 --core-wtp-mu-prior-sd=4)
      ;;
    joint-tight)
      extra+=(--mu-rep-sd=0.5 --dist-beta-v-sd=0.125)
      extra+=(--raw-u-sd-alpha=8 --raw-u-sd-beta=3.34782608695652)
      extra+=(--core-visibility-prior-multiplier=0.5)
      extra+=(--beta-intercept-sd=0.5 --beta-ink-effect-sd=0.125)
      extra+=(--beta-calendar-effect-sd=0.125 --beta-bracelet-effect-sd=0.125)
      extra+=(--wtp-value-utility-sd=0.00005 --core-wtp-mu-prior-sd=1)
      ;;
    joint-diffuse)
      extra+=(--mu-rep-sd=2 --dist-beta-v-sd=0.5)
      extra+=(--raw-u-sd-alpha=2.5 --raw-u-sd-beta=0.717391304347826)
      extra+=(--core-visibility-prior-multiplier=2)
      extra+=(--beta-intercept-sd=2 --beta-ink-effect-sd=0.5)
      extra+=(--beta-calendar-effect-sd=0.5 --beta-bracelet-effect-sd=0.5)
      extra+=(--wtp-value-utility-mean=-10 --wtp-value-utility-sd=4)
      extra+=(--lnorm-wtp-value-utility-prior=1 --core-wtp-mu-prior-sd=4)
      ;;
  esac
elif (( task == 14 )); then
  label=grouped-lambda
  fit_dir=${FIT_ROOT}/lambda/fits/grouped-sd0p25
  extra+=(--core-lambda-structure=1 --core-lambda-log-ratio-sd-prior=0.25)
elif (( task == 15 )); then
  label=arm-specific-lambda
  fit_dir=${FIT_ROOT}/lambda/fits/arm-sd0p25
  extra+=(--core-lambda-structure=2 --core-lambda-log-ratio-sd-prior=0.25)
elif (( task == 16 )); then
  label=student-t5
  fit_dir=${FIT_ROOT}/student-t/fits/student-t5
  extra+=(--core-type-distribution=1 --core-student-t-df=5)
elif (( task == 17 )); then
  label=finite-mixture
  fit_dir=${FIT_ROOT}/finite-mixture/fits/finite-mixture
  extra+=(--core-type-distribution=2)
else
  echo "Unsupported array task: ${task}" >&2
  exit 2
fi

mapfile -t fits < <(find "${fit_dir}" -maxdepth 1 -type f -name '*-chain[1-4]-1.csv' | sort)
if (( ${#fits[@]} != 4 )); then
  echo "Expected four retained fit CSVs for ${label}; found ${#fits[@]}" >&2
  exit 1
fi
fit_csvs=$(IFS=,; echo "${fits[*]}")

cd "${PROJECT_ROOT}"
module load -f R/4.2.0
export R_LIBS_USER=${R_LIBS_USER:-/home/edjee/R/x86_64-pc-linux-gnu-library/4.2}
mkdir -p temp/log "${OUTPUT_ROOT}/${label}"
Rscript --no-save --no-restore scripts/structural/generate-compact-gq.R \
  --workspace="${WORKSPACE}" \
  --fit-csvs="${fit_csvs}" \
  --output-path="${OUTPUT_ROOT}/${label}" \
  --output-basename="${label}-1250" \
  --distance-definition=assigned \
  --sm-evaluation-distance-m=1250 \
  --cmdstan-path="${CMDSTAN_PATH}" \
  --parallel-chains=4 --threads=4 \
  "${extra[@]}"
