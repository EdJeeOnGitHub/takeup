#!/usr/bin/env bash

set -euo pipefail

SOURCE_ROOT=${SOURCE_ROOT:-/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-model-robustness}
TIMESTAMP=${TIMESTAMP:-$(date -u +%Y%m%d-%H%M%S)}
OUTPUT_PATH=${OUTPUT_PATH:-${SOURCE_ROOT}/policy-model-robustness-complete-${TIMESTAMP}}
REPO_ROOT=${REPO_ROOT:-/home/edjee/projects/takeup-ed-refine-todos}

declare -A source=(
  [benchmark]=benchmark
  [private-distance-community-image]=private-distance-community-image
  [full-information]=full-information
  [exclude-dispersed]=exclude-dispersed
  [cluster-shock]=cluster-shock
  [tight-multinomial]=tight-multinomial
  [second-order-observability]=second-order-observability
  [grouped-lambda]=grouped-lambda
  [arm-lambda]=arm-lambda
  [student-t5]=student-t5
  [finite-mixture]=finite-mixture
)
models=(
  benchmark private-distance-community-image full-information exclude-dispersed
  cluster-shock tight-multinomial second-order-observability grouped-lambda
  arm-lambda student-t5 finite-mixture
)

mkdir -p "${OUTPUT_PATH}/audit"
for model in "${models[@]}"; do
  from=${SOURCE_ROOT}/${source[$model]}
  if [[ "${model}" == "cluster-shock" && ! -d "${from}" ]]; then
    from=${SOURCE_ROOT}/cluster-shock-mapped
  fi
  [[ -d "${from}" ]] || { echo "Missing source directory: ${from}" >&2; exit 1; }
  cp -al "${from}" "${OUTPUT_PATH}/${model}"
done

cd "${REPO_ROOT}"
Rscript --no-save --no-restore scripts/policy/audit-complete-model-robustness.R \
  "--package-path=${OUTPUT_PATH}"
: > "${OUTPUT_PATH}/audit/jobs.txt"
(
  cd "${OUTPUT_PATH}"
  find . -type f ! -path './audit/artifact-manifest.sha256' -print0 |
    sort -z | xargs -0 sha256sum > audit/artifact-manifest.sha256
)
echo "${OUTPUT_PATH}"
