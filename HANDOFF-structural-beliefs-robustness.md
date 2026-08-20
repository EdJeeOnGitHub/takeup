# Structural Beliefs Robustness Handoff

## Context

This implements the structural model robustness reruns requested in `ref-reports/edited-feedback-refine.md`, focused on alternative beliefs/observability inputs:

- Correct observability, matching Panel C of `presentations/rf-tables/main-specs/rf_discrete_belief_decomposition_tbl.tex`.
- Second-order beliefs, using the existing SOB counts.
- Missing first-order beliefs handled in the Bayesian model by marginalizing unobserved belief counts, not by Lee bounds.

The Overleaf paper file at `~/projects/overleaf/overleaf-takeup/ECM ReStud.tex` was not edited.

## Files Changed

- `run_takeup.R`
  - Adds `--beliefs-outcome=<fob|sob|correct-observability>`.
  - Adds `--beliefs-missing=<drop|latent>`.
  - Builds belief-count variants for default FOB, SOB, correct observability, and FOB with latent missing rows.
  - Adds model aliases:
    - `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS`
    - `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB`
    - `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB_LATENT_MISSING`
- `dist_structural_util.R`
  - Adds a model-specific `stan_data_preprocess` hook so each alias can swap the belief data passed to Stan.
- `stan_models/beliefs_data_sec.stan`
  - Adds `belief_observed`.
- `stan_models/beliefs_model_sec.stan`
  - Observed belief rows enter the binomial likelihood.
  - Latent/missing rows contribute no likelihood term, so their counts are marginalized under the fitted belief model.
- `stan_models/beliefs_generated_quantities.stan`
  - Adds posterior-predictive sampled counts for missing belief rows.
- `slurm_beliefs_robustness.sh`
  - New Midway Slurm array script for the three robustness runs.

## Run On Midway

From the repo root on Midway:

```bash
sbatch slurm_beliefs_robustness.sh 106
```

The array runs:

- task 0: `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_CORRECT_OBS`
- task 1: `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB`
- task 2: `STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB_LATENT_MISSING`

Default output path:

```bash
/project/akaring/takeup-data/data/stan_analysis_data
```

The script calls `postprocess_struct_models` after each fit. It uses `ITER=400` by default; override with:

```bash
ITER=800 sbatch slurm_beliefs_robustness.sh 106
```

## Validation Done

These checks passed locally:

```bash
Rscript -e 'parse("run_takeup.R"); parse("dist_structural_util.R")'
git diff --check
bash -n slurm_beliefs_robustness.sh
/home/ed/.cmdstan/cmdstan-2.36.0/bin/stanc --include-paths=stan_models stan_models/takeup_struct.stan --name=takeup_struct_smoke_model --o=/tmp/takeup_struct_smoke.hpp
/home/ed/.cmdstan/cmdstan-2.36.0/bin/stanc --include-paths=stan_models stan_models/takeup_struct_private_info.stan --name=takeup_struct_private_info_smoke_model --o=/tmp/takeup_struct_private_info_smoke.hpp
/home/ed/.cmdstan/cmdstan-2.36.0/bin/stanc --include-paths=stan_models stan_models/takeup_struct_no_generated_quantities.stan --name=takeup_struct_no_gq_smoke_model --o=/tmp/takeup_struct_no_gq_smoke.hpp
```

I did not complete a local R smoke fit because the active local R 4.6 renv library is missing `stringr`/`tidyverse`, while the older renv libraries are not compatible with R 4.6. The Stan syntax itself translates successfully with `stanc`.

## Follow-Up Checks After Midway Runs

- Confirm each job writes `dist_fit106_<MODEL>*.csv` and postprocessed RDS outputs.
- Run the existing HMC diagnostics path used for TODO 11 on each new alias.
- Compare belief probabilities/ATEs against the corresponding reduced-form tables before using the structural robustness results in paper text.
