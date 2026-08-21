#!/usr/bin/env bash
set -euo pipefail

specification=${1:-assigned}
case "$specification" in
  assigned|realized) ;;
  *) echo "Distance specification must be assigned or realized." >&2; exit 2 ;;
esac

threads=${TAKEUP_THREADS:-8}
threads_per_chain=$((threads / 4))
if ((threads_per_chain < 1)); then threads_per_chain=1; fi

output_dir="build/structural-fit/${specification}"
mkdir -p "$output_dir"

Rscript --no-save --no-restore --no-init-file \
  scripts/structural/sample-main-core.R \
  --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP \
  --input-path=data/stan_analysis_data \
  --workspace="${TAKEUP_STRUCTURAL_WORKSPACE:-data/stan_analysis_data/dist_fit104.RData}" \
  --output-path="$output_dir" \
  --stan-path=stan_models \
  --stan-file=takeup_struct_main_core.stan \
  --output-basename="dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP" \
  --distance-definition="$specification" \
  --chains=4 \
  --parallel-chains=4 \
  --threads-per-chain="$threads_per_chain" \
  --iter-warmup="${ITER_WARMUP:-1000}" \
  --iter-sampling="${ITER_SAMPLING:-1000}" \
  --adapt-delta="${ADAPT_DELTA:-0.999}" \
  --max-treedepth="${MAX_TREEDEPTH:-12}" \
  --metric="${METRIC:-diag_e}" \
  --seed="${SEED:-20260820}"
