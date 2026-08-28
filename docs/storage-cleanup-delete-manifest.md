# Storage cleanup: confirmed deletion manifest

Status: **completed; all nine listed files were deleted and verified absent**
Prepared: 2026-08-28
Deleted: 2026-08-28
Expected recovery: **29,933,619,985 bytes (27.88 GiB / 29.93 GB)**

This first cleanup tranche contains only large derived files that are not read
by the active slim-Stan, compact-GQ, focused structural-render, or sparse-policy
workflows. It deliberately excludes old structural robustness fits and other
historical `temp-data` objects whose compact replacements have not yet been
verified one by one.

## Files to delete

| Bytes | Path | Reason |
|---:|---|---|
| 12,182,088,457 | `temp-data/struct-postprocess/rvar_processed_dist_fit105_sm_draws_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS_1-4.rds` | Full social-multiplier draw cube. The current multiplier renderer constructs the required focused curve summary analytically from balanced slim-chain posterior parameters. `make structural-render` explicitly does not load this object. |
| 4,284,257,021 | `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_2.zip` | Legacy dense fit-86 policy prediction archive; replaced by sparse feasible-edge prediction from compact posterior parameters. |
| 3,909,091,022 | `optim/data/pred-demand-dist-fit86-static-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv` | Legacy dense fit-86 static-demand matrix; not used by the current policy workflow. |
| 3,221,581,038 | `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86.zip` | Legacy dense fit-86 policy prediction archive. |
| 2,209,833,237 | `optim/data/agg-log-full-many-csvs.zip` | Legacy aggregate policy CSV archive; current optimization consumes sparse predictions generated in `build/`. |
| 1,871,227,626 | `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_control.zip` | Legacy dense Control prediction archive. |
| 1,543,606,165 | `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_bracelet.zip` | Legacy dense Bracelet prediction archive. |
| 359,257,696 | `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_ink.zip` | Legacy dense Ink prediction archive. |
| 352,677,723 | `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_calendar.zip` | Legacy dense Calendar prediction archive. |

## Explicitly retained replacements and inputs

The deletion above must not include these active inputs:

| Bytes | Path | Purpose |
|---:|---|---|
| 8,013,008 | `optim/data/full-many-pots-experiment.rds` | Community--candidate-site geography for sparse prediction. |
| 12,213,428 | `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv` | Fixed experimental coverage target used by optimization. |
| 324,689 | `replication/inputs/policy/policy-bootstrap-parameters.csv` | Compact cluster-weighted policy parameters. |
| varies | `build/structural-fit/assigned/` | Current slim benchmark posterior chains and manifest. |
| varies | `temp-data/struct-postprocess/` excluding the single 12 GB file above | Focused structural summaries used by the current renderer. |

## Out of scope for this tranche

Do **not** delete the large fit-102/103 cluster-error objects, old fit-001
objects, files in `data/stan_analysis_data`, or named backup/archive directories
under this manifest. They require a separate appendix-robustness and provenance
audit even though several are likely removable later.
