# Workspace storage cleanup log

Date: 2026-08-28
Workspace: `/workspace/takeup`
Status: completed locally
Total removed: **136,694,945,742 bytes (127.31 GiB / 136.69 GB)**

This log records the cleanup performed after the benchmark structural workflow
was replaced by parameter-only Stan fits, compact generated quantities, focused
paper summaries, and sparse policy prediction. It is intended to support the
same cleanup on another server. Paths on that server must be resolved and
validated before deletion; do not assume its checkout or result transfers are
identical.

## Tranche 1: obsolete benchmark multiplier and dense policy data

Removed **9 files**, totaling **29,933,619,985 bytes (27.88 GiB)**:

- `temp-data/struct-postprocess/rvar_processed_dist_fit105_sm_draws_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS_1-4.rds`
- `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_2.zip`
- `optim/data/pred-demand-dist-fit86-static-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv`
- `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86.zip`
- `optim/data/agg-log-full-many-csvs.zip`
- `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_control.zip`
- `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_bracelet.zip`
- `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_ink.zip`
- `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/pred_demand_dist_fit86_calendar.zip`

Retained active policy inputs:

- `optim/data/full-many-pots-experiment.rds`
- `optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots/summ-agg-identity-experiment-target-constraint.csv`
- `replication/inputs/policy/policy-bootstrap-parameters.csv`

The more detailed first-tranche manifest is
`docs/storage-cleanup-delete-manifest.md`.

## Tranche 2: superseded untracked Stan outputs

Removed **136 files**, totaling **64,421,371,581 bytes (60.00 GiB)**, from
`data/stan_analysis_data/`.

The deleted set was defined as every **untracked** result file with extension
`.csv`, `.RData`, `.rds`, or `.zip`, except paths matching one of the retained
legacy specifications below. Git-tracked data-submodule files and the helper
script were not removed.

Retained legacy match patterns:

- `dist_fit102`
- `dist_fit103`
- `INDIV_DIST_COMMUNITY_FP_INDIV_VIS`
- `INDIV_DIST_INDIV_FP`
- `NO_OUTLIERS`
- `_SOB`

These correspond to four active specifications that have not yet been
streamlined:

1. Individual travel costs
2. Individual distance observed by peers
3. Excluding geographically dispersed communities
4. Perceived community observability

The exact retained files and migration condition are documented in
`data/stan_analysis_data/README.md`. On another server, first verify that its
current slim fits, compact GQs, exact multiplier draws, focused summaries, and
policy-model robustness bundle have been returned successfully. Preserve Git-
tracked inputs even if their extensions match the result-file extensions.

## Tranche 3: large legacy postprocessing and presentation caches

Removed **5 files and 3 cache directories**, totaling
**34,795,193,774 bytes (32.41 GiB)**.

Legacy derived files:

- `temp-data/rvar_processed_dist_fit102_cluster_error_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_1-2.rds`
- `temp-data/rvar_processed_dist_fit103_cluster_error_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS_1-2.rds`
- `temp-data/rvar_processed_dist_fit103_cluster_error_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS_1-2.rds`
- `temp-data/rvar_processed_dist_fit001_sm_draws_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_BELIEFS_SUBMODEL_1.rds`
- `temp-data/rvar_processed_dist_fit001_sm_draws_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_BELIEFS_CONSTANT_1.rds`

Rebuildable presentation caches:

- `presentations/takeup-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB95-figure-cache/`
- `presentations/takeup-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB95-table-cache/`
- `presentations/takeup-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIER_FOB95-figure-cache/`

## Verification and recovery notes

- The workspace fell from approximately 153 GB to approximately 32 GB.
- `.git/` was deliberately left untouched.
- Current slim benchmark chains, compact GQs, focused result bundles, paper
  artifacts, and sparse policy inputs were retained.
- The former untracked canonical workspace
  `data/stan_analysis_data/dist_fit104.RData` was included in tranche 2. The
  pipeline was subsequently refactored to rebuild the required `models` and
  `stan_data` objects directly from the current loaders with
  `make structural-workspace`; restoration is no longer required.
- Deleted files were untracked or rebuildable caches and cannot be recovered
  from local Git. Recovery requires regeneration or another server/archive.
- Before repeating this cleanup elsewhere, compare file sizes, confirm current
  result-bundle manifests, and generate a dry-run list. Do not use broad globs
  without inspecting the resolved targets.

## Midway2 active-robustness migration

On 2026-08-28, four active robustness specifications were removed from the
legacy-fit exception list after a separate server-side migration:

- Individual travel costs: lossless parameter-only subset, 800 draws.
- Individual distance observed by peers: lossless parameter-only subset, 800
  draws.
- Excluding geographically dispersed communities: latest main-core refit,
  4,000 draws.
- Perceived community observability: latest main-core refit, 4,000 draws.

All four have four compact GQ chains and fit diagnostics under
`data/stan_analysis_data/streamlined-active-robustness/`. The main-core refits
had zero divergences and zero maximum-treedepth transitions. Multiplier checks
passed against the legacy draws, except that the historical SOB compact GQ was
diagnosed as incorrectly hard-coding first-order belief coefficients. The SOB
gate instead used the legacy posterior with a corrected second-order GQ and
then passed.

After downstream policy fit paths were changed to the replacements, the
manifest-gated deletion removed exactly 16 active legacy CmdStan CSVs totaling
**11,094,305,455 bytes (10.33 GiB)**. The server audit directory contains the
deletion manifest, pre-deletion SHA-256 hashes, deletion timestamps, and
post-deletion existence checks. Older fit generations and archive files were
not included and remain untouched.

The downstream policy check used all 4,000 posterior draws for both new
main-core fits and solved all five policy scenarios. Slurm job families
`48878594`--`48878597` and `48878598`--`48878601` completed successfully; the
provenance-only summary refreshes were `48879246` and `48879247`. The promoted
policy bundles passed `scripts/policy/audit-streamlined-policy-replacements.R`.
The previous policy directories were retained in the dated superseded-policy
archive rather than deleted.

## Tranche 4: superseded backups, transfer archives, and debug intermediates

Removed **5 files and 3 directories**, totaling **1,619,928,183 bytes
(1.51 GiB)**:

- `temp-data/karim-fits.zip`
- `temp-data/backup_processed_dist_fit66.zip`
- `temp-data/tidy_processed_dist_fit.zip`
- `temp-data/debug-analysis-data.csv`
- `temp-data/debug-clean-census-data.csv`
- `temp-data/old-75-stuff/`
- `temp-data/archive/`
- `temp-data/backup/`

The named backup/archive directories included copies of fit-71 and fit-75
workspaces also present at the top level. No active Make or manuscript pipeline
reference was found for the transfer archive or debug intermediates.

## Tranche 5: migrated top-level processed outputs

Removed **194 files**, totaling **5,924,832,219 bytes (5.52 GiB)**, from the
top level of `temp-data/`.

The validated selection used basenames beginning with `processed_dist`,
`rf_prior_processed_dist`, `rvar_processed_dist`, `tidy_processed_dist`,
`dist_fit`, or `backup_processed`. The sole active exception was retained:

- `temp-data/tidy_processed_dist_prior95_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB_1-4.rds`

That fit-95 prior-predictive object remains an explicit input to the gray prior
curve in Figure 4. All deleted objects had no exact basename reference in the
active Make, R, scripts, appendix, or replication pipeline after excluding
archived notebooks, review notes, and old presentation sources. See
`docs/temp-data-legacy-processed-audit.md` for the audit and migration status.
