# Audit of top-level legacy processed outputs

Date: 2026-08-28
Status: **completed; the 194 superseded candidates were deleted and the active
prior-predictive exception was verified present**
Scope: top-level files in `temp-data/` whose basenames begin with
`processed_dist`, `rf_prior_processed_dist`, `rvar_processed_dist`,
`tidy_processed_dist`, `dist_fit`, or `backup_processed`.

## Result

- Files audited: **195**
- Total size: **5,983,248,528 bytes (5.57 GiB)**
- Required by the active Make/paper pipeline: **1 file, 58,416,309 bytes**
- Superseded candidates deleted: **194 files, 5,924,832,219 bytes (5.52 GiB)**

The audit searched exact basenames across `Makefile`, `_targets.R`, `scripts/`,
`R/`, `replication/`, `appendix/`, and the root READMEs. References confined to
`archive/`, old presentation sources, or review notes were treated as historical
rather than active pipeline dependencies.

## Active exception to retain

`temp-data/tidy_processed_dist_prior95_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB_1-4.rds`

`scripts/structural/render-paper.R` still reads this 55.7 MiB fit-95 prior-
predictive object to draw the gray prior curve in the social-multiplier figure.
It is also marked required for Figure 4 in `replication/data-manifest.csv`.
The current slim posterior and analytical multiplier renderer have replaced the
posterior side of this workflow, but the prior-predictive curve has not yet been
migrated to a compact main-core prior-GQ result.

## Superseded set

The other 194 files have no exact active-source basename reference. They consist
of historical full/lite processed workspaces, old reduced-form Stan outputs,
fit-95/96/98/100/101 posterior objects, obsolete fit-001 submodel outputs, and
intermediate rate-of-change, treatment-effect, and cluster-error objects.

Current paper-facing replacements are the Make-generated RF artifacts, slim
benchmark chains and compact GQs, focused `temp-data/struct-postprocess/`
summaries, exact 1.25 km multiplier-draw bundles, and complete policy-model
robustness summaries.

Before deleting this set on another machine, reproduce the selection by exact
top-level basename prefixes, exclude the active prior-predictive file above,
and inspect the resolved list and total size. Do not recursively apply this rule
inside active `temp-data/main-core-*`, `assigned-distance-comparison`,
`asymmetric-observability-comparison`, or `struct-postprocess` directories.
