# Policy population-weighting rerun audit

Date: 2026-09-06. Read-only audit of estimation and allocation code and existing outputs. No estimation, optimization, or manuscript edits performed.

## Critical finding: existing population weights are wrong

`scripts/policy/run-population-cost.R:116`, `scripts/policy/generate-distance-cap-table.R:59`, and `scripts/policy/run-cost-sensitivity.R:43` aggregate `census$num.individuals` by community. However, `data/takeup_census.RData` contains one row per individual, with household size repeated on each individual row. The cleaning code constructs this variable using household-level `n()` inside `mutate()` (`rct-design-fieldwork/clean_hh_census_data.R:627`). Summing it counts household size repeatedly.

Direct inspection for the 144 policy communities found:

- 39,301 census rows and 39,301 distinct `KEY.individ` values, all with `age.census >= 18`.
- 19,601 distinct household keys.
- Summing `num.individuals` gives 96,052, exactly the population recorded in both existing policy input audits.

Use unique census adults per community, joined by `cluster.id`. This changes relative community weights as well as aggregate participant counts. Consequently the existing population-weighted 25-site saving and $261/$483 break-even thresholds are not verified results under correct adult-population weighting. Do not simply rescale the costs. Rerun allocations and accounting.

## Required work

| Component | Existing evidence | Required action |
|---|---|---|
| Benchmark at 3.5 km | Population runner has 1,600 successful posterior draws, but only endogenous Control and Bracelet, with incorrect weights | Correct population construction and run all five scenarios: Control, Bracelet, fixed-at-0.5-km Control, fixed-at-0.5-km Bracelet, and no social image. Recompute experimental Control target within each draw using correct adult weights. Report paired differences and uncertainty |
| Ten alternative structural models | Current complete table uses equal-weight/common-target allocation pipeline | Rerun five allocation scenarios for each model with correct weights and a consistent target convention. No structural refitting required |
| Cluster weighted-likelihood inference | 999 successful weighted modes have Control/Bracelet allocations with incorrect population weights | Reuse modes and rerun allocations. Include both fixed-return scenarios if reporting bootstrap uncertainty for the multiplier's allocation contribution |
| Larger distance caps | Existing robustness table has 3.5, 4.5, 5.5, and 10 km panels from the legacy workflow | Rerun benchmark at 4.5, 5.5, and 10 km with corrected population weights and the same target/scenarios as the main exercise. No need to cross every alternative model with every cap merely to retain current coverage |
| Short-cap diagnostic | 2.5, 2.75, 3, 3.25, and 3.5 km, evaluated at component-wise posterior median parameters, uses incorrect population construction | Rerun Control/Bracelet allocations with corrected weights. Feasible links, shareable sites, and demand-free geographic minimum are independent of population weights and reusable. Do not assume the old finding that both regimes attain the geographic floor still holds |
| Observability under consolidation | 999-mode population runs include 0%, 50%, and 100% attenuation | Rerun these allocations with corrected weights if retaining this sensitivity. The 0% case duplicates baseline and can be reused within the new run |
| Break-even accounting | Posterior and cluster-mode cost grids exist, based on incorrect populations | Recompute bracelet costs, expected participant travel, thresholds, and probabilities from corrected allocations. This is accounting after allocation, not a new resource-cost minimization |
| Distance-distribution figure and allocation maps | Renderer reads stored scenario allocations | Generate/save corrected population-weighted allocations at the chosen representative parameter vector and rerender. Current population runner saves summaries, not full edge assignments. State whether the plotted distance distribution weights communities or adults |
| Demand curves | Structural take-up as a function of distance | Reuse if parameter draws/specification are unchanged. Population weighting of the planner objective does not alter these curves |

## Alternative-model inventory

The current complete catalog has ten alternatives (plus benchmark):

1. Individual travel costs (`private-distance-community-image`).
2. Individual distance observed by peers (`full-information`).
3. Excluding geographically dispersed communities (`exclude-dispersed`).
4. Unobserved community heterogeneity (`cluster-shock`).
5. Correct classification of take-up (`tight-multinomial`).
6. Perceived community observability (`second-order-observability`).
7. Social-image effects by public-signal status (`grouped-lambda`).
8. Social-image effects by treatment arm (`arm-lambda`).
9. Heavy-tailed intrinsic motivation (`student-t5`).
10. Mixture intrinsic motivation (`finite-mixture`).

See `scripts/policy/audit-complete-model-robustness.R` and `temp-data/policy-model-robustness-complete-streamlined-20260828/`. The latter has summaries/status/provenance for all rows but does not contain the parameter CSVs or edge-demand matrices checked in this audit. Some parameters exist in other local model directories. Recover the exact production inputs/caches where available, otherwise regenerate policy predictions from existing structural fits. Do not confuse older binary correct-observability fits with the promoted tight multinomial model.

## Pipeline changes before running

- Share a validated adult-population vector across runners. Assert 144 communities, 39,301 unique adults, and correct community joins. Existing validation checks internal arithmetic but misses the population error.
- Update the model-aware allocation pipeline rather than routing all models through the current population runner: `run-population-cost.R` forces `model_family = "gaussian"`. The alternative predictor already handles household distance, community shocks, mixture, and other families.
- `optimize-cluster-bootstrap.R` currently uses unit community weights and an across-draw mean target. Change both weighting and target handling explicitly. Preserve a common experimental Control target across counterfactual scenarios within each draw. Specify whether alternative models preserve their own experimental predictions or a benchmark target; do not silently inherit the old fixed CSV.
- Retain model-specific household geography and community-shock mappings when constructing experimental-allocation targets. Candidate-site predictions can be reused at unchanged distances when available; new larger-cap edges need predictions.
- Write new outputs to separate directories. The optimizer skips existing allocation files without checking their weighting/target definition, so reusing old paths could silently preserve stale results.
- Save assignment-level outputs needed for figures and resource accounting, plus solver status, target attainment, parameter provenance, population checksum, and cap.
- Retain equal-weight results as robustness. A clean comparison to new population weighting should hold the target convention fixed, because the current rewrite uses a fixed target while the population runner uses draw-specific targets.

## Scope

This is a policy allocation and reporting rerun, not a rerun of structural estimation, social-multiplier estimation, or reduced-form results. Core coverage is 11 models × 5 policy scenarios at 3.5 km, plus benchmark larger-cap checks, the short-cap diagnostic, and the retained inference/pooling sensitivities. Existing structural draws can be reused. No jobs launched during this audit.

## Server-agent handoff: implementation required before submission

The commands below specify the intended production interface. **They are not ready to submit against the unmodified repository.** In particular, `DISTANCE_CAP`, `POPULATION_WEIGHTING`, and `TARGET_MODE` below are required additions to the stage wrappers, not currently supported controls. The server agent should implement, smoke-test, and then run this handoff. Do not submit structural sampling jobs.

### Existing Bash entry points and required edits

| Existing file | Role and required correction |
|---|---|
| `hpc/policy/slurm_policy_model_robustness.sh` | Main prepare/predict/optimize/summarize stage worker. Retain its model-specific extraction adapters. Add population/target/cap controls, use fresh output roots, and replace the empty `benchmark` prepare branch with preparation of the existing assigned-distance benchmark draws. Pass explicit mapping workspace for community shocks and preserve `dist_fit104.RData` household geography for the two individual-distance models. Its current `ROOT` is hard-coded, and its module/repository defaults are host-specific |
| `hpc/policy/submit_policy_model_robustness.sh` | Current six-model list is stale and includes obsolete `correct-observability`. Replace with the eleven-row list below, or use the explicit submission loop below after patching the stage worker |
| `hpc/policy/slurm_policy_cluster_bootstrap.sh` | Reuse existing 999 exponential weighted modes. Correct target and population handling in optimize. Its population stage calls the buggy population runner and must also be corrected |
| `hpc/policy/submit_policy_cluster_bootstrap.sh` | Existing dependency chain prepare → predict → optimize array → summarize → population. Can be used after the stage worker is patched; explicit output paths are mandatory |
| `hpc/policy/assemble_complete_model_robustness.sh` | Packages all eleven model directories and runs the completeness audit. Use the new population output root. Packaging uses hard links: treat the package as immutable after creation |
| `scripts/policy/create-counterfactuals.sh` | Legacy driver, not the recommended new runner. Its old target/optimization path should not be used for the population-weighted distance-cap jobs |

### Exact model inputs

Server data root in the existing worker: `/project/akaring/takeup-data/data/stan_analysis_data`. Resolve actual mount/repository locations on the server before running. Paths below are relative to that root except the benchmark. Do not replace missing current inputs with older fits without documenting it.

| MODEL_ID | Fit file pattern | Family / extraction options | Existing complete-package draw count |
|---|---|---|---:|
| `benchmark` | Assigned-distance slim chains `dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-{1,2,3,4}.csv`; existing preparer defaults to repository `build/structural-fit/assigned` | Gaussian; use `prepare-baseline-posterior.R` and retain chain/iteration IDs | 200 balanced draws in main rewrite; separate cost run has 1,600 |
| `private-distance-community-image` | `streamlined-active-robustness/private-distance-community-image/fits/private-distance-community-image-slim-chain{1,2,3,4}-1.csv` | `private_distance_community_image`; household workspace `dist_fit104.RData` | 800 |
| `full-information` | `streamlined-active-robustness/full-information/fits/full-information-slim-chain{1,2,3,4}-1.csv` | `full_information`; household workspace `dist_fit104.RData` | 800 |
| `exclude-dispersed` | `streamlined-active-robustness/exclude-dispersed/fits/exclude-dispersed-slim-chain{1,2,3,4}-1.csv` | Gaussian | 4,000 |
| `cluster-shock` | `main-core-cluster-shock-production/dist_fit106_MAIN_CORE_SHOCK_SD0.1_chain{1,2,3,4}-1.csv` | Gaussian; `--include-cluster-shock 144`; cluster map from matching `dist_fit106.RData` | 4,000 |
| `tight-multinomial` | `main-core-report-distance-priors/tight/dist_fit106_MAIN_CORE_chain{1,2,3,4}-1.csv` | `asymmetric_conditional`; `--include-asymmetric` | 1,600 |
| `second-order-observability` | `streamlined-active-robustness/second-order-observability/fits/second-order-observability-slim-chain{1,2,3,4}-1.csv` | Gaussian; `--beliefs-order 2` | 4,000 |
| `grouped-lambda` | `main-core-lambda-identification/fits/grouped-sd0p25/grouped-sd0p25-chain{1,2,3,4}-1.csv` | Gaussian, grouped lambda; `--include-lambda grouped` | 1,600 |
| `arm-lambda` | `main-core-lambda-identification/fits/arm-sd0p25/arm-sd0p25-chain{1,2,3,4}-1.csv` | Gaussian, arm lambda; `--include-lambda arm` | 1,600 |
| `student-t5` | `main-core-student-t-robustness/fits/student-t5/student-t5-chain{1,2,3,4}-1.csv` | `student_t5` | 1,600 |
| `finite-mixture` | `main-core-finite-mixture-robustness-800/fits/finite-mixture/finite-mixture-chain{1,2,3,4}-1.csv` | `finite_mixture`; `--include-finite-mixture` | 3,200 |

Use all retained draws for alternatives. For benchmark, preserve the current 200 balanced draws for a direct comparison first. If the 1,600-draw cost input is verified as the same underlying fit, use a common 1,600-draw set for final benchmark allocations, costs, and distance robustness, retaining the 200-draw comparison. If it is not the same fit, regenerate the benchmark draw set from the intended assigned-distance chains. Record the resolution; do not splice 200-draw allocation estimates and unrelated 1,600-draw cost estimates.

Weighted modes: inspect the existing worker's candidates under the data root, `main-core-exponential-cluster-weight-999` and `main-core-weighted-modes`, as well as the preparer's older default `main-core-exponential-cluster-weights`. Select the directory with the 999 complete exponential, assigned-distance statuses. Existing canonical CSV name is `policy-bootstrap-parameters.csv`.

### Required R implementation and target contract

1. Add a shared validated adult-count loader (for example in `R/policy/cost-sensitivity.R`) and use it in all three buggy runners. Count unique `KEY.individ` within the 144 matched `cluster.id` values. Reject duplicates/inconsistent community assignments rather than silently summing household size.
2. Extend the model-aware preparation/prediction pipeline to save experimental Control take-up by community and draw, including household-specific distance and community shocks where applicable. Produce a target CSV with `draw`, `replicate`, `target_expected_adults`, `target_rate`, `population_total`, and model/input provenance.
3. Adopt this explicit default for the new exercise: each model/draw preserves **its own predicted number of adults treated under the original experimental pairings and Control observability**. Use that same target across all five policy scenarios and all caps for that model/draw. The previous common fixed target can be retained as a separately labeled sensitivity, but must not enter these runs silently.
4. Extend `optimize-cluster-bootstrap.R` to consume the adult vector and draw-matched target instead of `rep(1, num_villages)` and the across-draw target mean. Required CLI additions: `--population-weighting=adult-census`, `--distance-data=...`, and `--target-mode=draw-specific-experimental-control`. The stage worker must create and pass its own `${OUTPUT_PATH}/policy-experimental-targets.csv` as `--target-csv`.
5. Propagate `DISTANCE_CAP` into prediction and provenance. The five array IDs are exactly: 1 `control`, 2 `bracelet`, 3 `static-control`, 4 `static-bracelet`, 5 `suppress-reputation`. Experimental allocation is a reference row, not a sixth optimization.
6. Update `summarize-model-results.R`, `summarize-cluster-bootstrap.R`, `standardize-baseline-posterior.R`, and the renderers as necessary so the displayed take-up/targets are adult-weighted. Preserve paired draw-level contrasts, undefined-equilibrium flags, and target-infeasibility shares. Do not describe an infeasible fallback as meeting the target.
7. Make `run-population-cost.R` reuse the benchmark's selected canonical draws and corrected allocations for cost accounting (or verify exact agreement if solving again). Keep equal-weight comparisons at the same draw-specific target convention. Retain 0/50/100% pooling attenuation in the 999-mode analysis. Do not apply its forced Gaussian adapter to alternative models.
8. Make `generate-distance-cap-table.R` use corrected adult counts. Geographic floors are reusable, but previous equality of the structural allocation and geographic floor is a result to recheck, not an assertion to hard-code in `validate-population-cost.R`.

### Production submission commands (after the preceding changes pass smoke tests)

Run from the server checkout. These commands use existing stage-script filenames with the required new environment controls above. The server agent must replace site-specific paths/partition if different, and record the final literal commands and job IDs in the returned run log.

```bash
set -euo pipefail
export REPO_ROOT="$PWD"
export PROJECT_ROOT="$PWD"
export DISTANCE_DATA=/project/akaring/takeup-data/optim/data/full-many-pots-experiment.rds
export RUN_ROOT=/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-adult-population-20260906
export NUM_CORES=8
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export POPULATION_WEIGHTING=adult-census
export TARGET_MODE=draw-specific-experimental-control
mkdir -p temp/log "$RUN_ROOT"

# These are the current Midway3 settings used by the existing bootstrap launcher.
# Resolve the actual server partition/account before submitting.
SBATCH_SITE=(--partition=caslake --account=pi-akaring --cpus-per-task=8)
STAGE_SCRIPT=hpc/policy/slurm_policy_model_robustness.sh
MODELS=(benchmark private-distance-community-image full-information
  exclude-dispersed cluster-shock tight-multinomial second-order-observability
  grouped-lambda arm-lambda student-t5 finite-mixture)

submit_model_cap() {
  local model_id=$1 cap_m=$2 prep pred opt summ
  export MODEL_ID="$model_id"
  export DISTANCE_CAP="$cap_m"
  export OUTPUT_PATH="$RUN_ROOT/cap-${cap_m}/${model_id}"
  prep=$(sbatch --parsable "${SBATCH_SITE[@]}" --export=ALL,STAGE=prepare "$STAGE_SCRIPT")
  pred=$(sbatch --parsable "${SBATCH_SITE[@]}" --dependency="afterok:${prep}" \
    --export=ALL,STAGE=predict "$STAGE_SCRIPT")
  opt=$(sbatch --parsable "${SBATCH_SITE[@]}" --dependency="afterok:${pred}" \
    --array=1-5 --export=ALL,STAGE=optimize "$STAGE_SCRIPT")
  summ=$(sbatch --parsable "${SBATCH_SITE[@]}" --dependency="afterok:${opt}" \
    --export=ALL,STAGE=summarize "$STAGE_SCRIPT")
  echo "$model_id $cap_m prepare=$prep predict=$pred optimize=$opt summarize=$summ" \
    >> "$RUN_ROOT/jobs.txt"
}

for model_id in "${MODELS[@]}"; do
  submit_model_cap "$model_id" 3500
done
for cap_m in 4500 5500 10000; do
  submit_model_cap benchmark "$cap_m"
done

# Separate inference sensitivity, reusing weighted modes rather than refitting.
unset MODEL_ID
export DISTANCE_CAP=3500
export OUTPUT_PATH="$RUN_ROOT/cluster-weighted"
export POLICY_REVIEW_OUTPUT="$RUN_ROOT/costs"
export POLICY_WORK_PATH="$RUN_ROOT/work"
export TABLE_PATH="$RUN_ROOT/cluster-weighted-summary.tex"
NUM_REPLICATES=999 bash hpc/policy/submit_policy_cluster_bootstrap.sh \
  > "$RUN_ROOT/cluster-jobs.txt"
```

Before production, run a bounded smoke test in a different output root: two draws for each model, all five scenarios, and benchmark caps 3.5 and 10 km. The predictor's existing `MAX_DRAWS` can bound this; ensure the target rows, expected-draw audit, and summary agree with the subset. Include a tight-multinomial infeasibility case if present. Confirm new flags actually reach the R code and are recorded in outputs. These tests do not require structural refits.

### Remaining jobs and reporting commands

After model jobs complete, execute the following in an allocated server session with its R dependencies loaded. `BENCHMARK_CANONICAL_CSV` must be set to the verified common benchmark draw set produced above (normally `cap-3500/benchmark/policy-model-parameters.csv` after standardization). These are existing R options, but the population bug must be fixed first.

```bash
export BENCHMARK_CANONICAL_CSV="$RUN_ROOT/cap-3500/benchmark/policy-model-parameters.csv"
Rscript scripts/policy/run-population-cost.R \
  "--parameter-csv=$BENCHMARK_CANONICAL_CSV" --parameter-type=canonical \
  --analysis-id=baseline-posterior "--distance-data=$DISTANCE_DATA" \
  "--output-path=$RUN_ROOT/costs" "--work-path=$RUN_ROOT/work" \
  --cores=8 --solver=auto --include-legacy=true --pooling-rhos=0

Rscript scripts/policy/generate-distance-cap-table.R \
  "--parameter-csv=$BENCHMARK_CANONICAL_CSV" --parameter-type=canonical \
  "--distance-data=$DISTANCE_DATA" \
  "--csv-path=$RUN_ROOT/costs/policy-distance-cap-diagnostics.csv" \
  "--table-path=$RUN_ROOT/policy-distance-cap-feasibility.tex"

Rscript scripts/policy/assemble-population-cost.R \
  "--input-path=$RUN_ROOT/costs" "--table-path=$RUN_ROOT/tables" \
  "--figure-path=$RUN_ROOT/figures"
Rscript scripts/policy/validate-population-cost.R \
  "--input-path=$RUN_ROOT/costs" "--distance-data=$DISTANCE_DATA"
Rscript scripts/policy/render-model-scenario-table.R \
  "--input-root=$RUN_ROOT/cap-3500" \
  "--output=$RUN_ROOT/tables/optim-policy-model-scenarios.tex"
Rscript scripts/policy/audit-complete-model-robustness.R \
  "--package-path=$RUN_ROOT/cap-3500"
```

For Anne's figure, solve Control and Bracelet once at the common benchmark's component-wise posterior median parameters, using corrected population weights and the corresponding experimental target. Save assignments in the existing renderer's `median-allocations/{control,bracelet}/replicate-0001.rds` layout, retaining its allocation fields. The renderer command is:

```bash
Rscript scripts/policy/render-paper-figures.R \
  "--parameter-csv=$BENCHMARK_CANONICAL_CSV" --parameter-draw=median \
  --allocation-draw=1 "--policy-path=$RUN_ROOT/cap-3500/benchmark" \
  "--distance-data=$DISTANCE_DATA" "--output-path=$RUN_ROOT/figures"
```

Keep the distribution figure community-weighted to preserve its current meaning (distribution of community-to-site distances), but state that allocations preserve adult-weighted coverage. This is different from the adult-weighted mean distance reported in the new table. Label both explicitly. The demand-curve figure itself does not change solely because the planner's weights change.

### Deliverables and acceptance criteria

- Return code changes, the literal Bash commands/job IDs, and an input manifest with source fit paths, hashes, draw selection, population counts, distance cap, target convention, and solver settings. Archive old outputs separately.
- All eleven models × five scenarios at 3.5 km must have complete draw accounting. Benchmark larger caps must use identical draws and targets. Target infeasibility is a reported economic outcome, not an optimization failure to discard.
- Independently verify 39,301 unique adults across the canonical 144 communities. Check coverage and site assignment constraints directly from saved assignments. Record solver optimality/gaps, not merely exit success.
- Verify corrected Control/Bracelet cost-run allocations agree with the main benchmark draws and target definition. Report the change from old incorrect weights. Do not claim the old $261/$483 thresholds or 25-site saving survived until recomputed.
- Supply benchmark table and paired contrasts, full alternative-model table, distance-cap table and diagnostic, bootstrap and pooling summaries, break-even grid/figure, distance-density figure, and allocation maps. Preserve full numerical CSVs/RDS and figure input data alongside TeX/PDF outputs.
- Produce a short change summary with equal-weight versus corrected adult-weighted estimates on matched draws/target conventions. Flag any remaining differences in draw provenance or solver tie-breaking that affect distance/cost summaries.
- Do not update or publish the Overleaf manuscript as part of the server run. Return reviewable artifacts for incorporation here.
