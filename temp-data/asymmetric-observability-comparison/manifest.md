# Asymmetric observability comparison

- Git commit: `b59675c40305c3fcaab3c0af1dbe36be314362d6`
- Distance definition: `assigned`
- Posterior sampling rerun: no
- 1.5 km GQ rerun: no; existing compact GQs were used
- Pre-emptive exact-1.25 km GQ root: `/project/akaring/takeup-data/candidate-asymmetric-observability-1250-20260825`
- Pre-emptive exact-1.25 km GQ status: completed successfully; jobs 55033251, 55033252, and 55033253; four chains each; generated quantities only

## Specifications and source directories

### Tight multinomial

- Options: `observation_model=1; recognition_structure=0; report_structure=0; report_arm_dist_prior_scale=0.10`
- Fit directory: `/project/akaring/takeup-data/data/stan_analysis_data/main-core-report-distance-priors/tight`
- Saved 500/1500/2500m GQ directory: `/project/akaring/takeup-data/data/stan_analysis_data/main-core-report-distance-priors/tight/gq`

### F3 two-stage conditional

- Options: `observation_model=1; recognition_structure=2; report_structure=2; report_arm_dist_prior_scale=0.25`
- Fit directory: `/project/akaring/takeup-data/data/stan_analysis_data/main-core-observability-ladder/f3`
- Saved 500/1500/2500m GQ directory: `/project/akaring/takeup-data/data/stan_analysis_data/main-core-observability-ladder/f3/gq`

### U3 two-stage unconditional

- Options: `observation_model=2; recognition_structure=1; report_structure=2; report_arm_dist_prior_scale=0.25`
- Fit directory: `/project/akaring/takeup-data/data/stan_analysis_data/main-core-observability-ladder/u3`
- Saved 500/1500/2500m GQ directory: `/project/akaring/takeup-data/data/stan_analysis_data/main-core-observability-ladder/u3/gq`

