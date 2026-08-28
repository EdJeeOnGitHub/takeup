# Legacy structural fit exceptions

The active benchmark and most structural robustness models use the current
parameter-only Stan fit plus a separate compact generated-quantities program.
Their paper-facing draws and summaries live in focused `temp-data/main-core-*`,
`temp-data/assigned-distance-comparison`, and policy-robustness bundles; large
historical Stan CSVs and processed workspaces are not retained.

The former benchmark input `data/stan_analysis_data/dist_fit104.RData` is no
longer required. The Make pipeline rebuilds the minimal `models` and `stan_data`
objects directly from the current loaders into
`build/structural-workspace/main-core-input.RData`. Fit-cache manifests include
the workspace hash, so posterior draws produced from an older data snapshot
cannot be silently reused.

Four active robustness specifications have not yet been refitted in the new
streamlined format. Their legacy sources remain under
`data/stan_analysis_data/`:

1. **Individual travel costs** (`INDIV_DIST_COMMUNITY_FP_INDIV_VIS`)
2. **Individual distance observed by peers** (`INDIV_DIST_INDIV_FP`)
3. **Excluding geographically dispersed communities** (`NO_OUTLIERS`)
4. **Perceived community observability** (`SOB`)

The retained legacy artifacts are:

- `dist_fit102.RData` and `dist_fit102.zip`
- `dist_fit103.RData` and `dist_fit103.zip`
- `dist_fit105_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS-[1-4].csv`
- `dist_fit106_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_SOB.rds`
- the four small `rvar_processed_*_SOB_1-4.rds` summaries

The fit-102/103 archives are retained conservatively because they contain the
legacy individual-distance variants; the server-returned policy provenance and
focused multiplier-draw bundles are the authoritative paper-facing records.

Remove these exceptions only after the four specifications are refitted with
the parameter-only model and compact GQ, and their multiplier and policy
summaries have been checked against the current paper artifacts.
