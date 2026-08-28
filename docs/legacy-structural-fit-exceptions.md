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

The four formerly active exceptions were migrated on Midway2 on 2026-08-28:

1. **Individual travel costs** (`INDIV_DIST_COMMUNITY_FP_INDIV_VIS`)
2. **Individual distance observed by peers** (`INDIV_DIST_INDIV_FP`)
3. **Excluding geographically dispersed communities** (`NO_OUTLIERS`)
4. **Perceived community observability** (`SOB`)

Their canonical server fits now live beneath
`data/stan_analysis_data/streamlined-active-robustness/<specification>/fits/`.
The two individual-distance fits are lossless parameter-column subsets of the
validated legacy posterior draws (800 draws each). `NO_OUTLIERS` and `SOB` were
refitted with the latest main-core model (4,000 draws each), with zero
divergences and zero maximum-treedepth transitions. Compact 1.25 km GQs and the
full audit are stored alongside the fits.

All multiplier equivalence gates passed. The historical compact SOB GQ was
found to hard-code first-order belief coefficients even for the second-order
model. The deletion audit therefore compares the new SOB fit against the old
posterior evaluated with the corrected second-order main-core GQ. Those
reader-facing medians differ by at most 0.009 across arms. The original and
streamlined posterior parameter summaries also agree closely.

The active policy launcher now reads these streamlined paths. On Midway2, the
manifest-gated cleanup removed the 16 superseded active fit CSVs listed in
`streamlined-active-robustness/audit/legacy-deletion-manifest.csv`: 10.33 GiB
(11,094,305,455 bytes). Each file was SHA-256 hashed before deletion; the
per-file receipt is `legacy-deletion-receipt.csv` in the same audit directory.

The policy downstream was then regenerated with all 4,000 draws for
`exclude-dispersed` and `second-order-observability`. The latter was necessary
because its superseded policy bundle had also extracted first-order rather than
second-order belief coefficients. Both replacements passed the policy migration
audit: five scenarios and 4,000 allocations per scenario, no prediction errors,
and explicit assigned-distance/fit-path provenance. Correcting SOB changes the
median endogenous Control-versus-Bracelet saving from 21 to 23 PoTs (95% interval
11--33), leaving the qualitative policy conclusion intact. The canonical server
policy directories now contain the streamlined versions; the old policy bundles
were moved, not deleted, to `policy-model-robustness-superseded-20260828/`.

The following older archives and fit generations remain exceptions and were
not part of that deletion manifest:

- `dist_fit102.RData` and `dist_fit102.zip`
- `dist_fit103.RData` and `dist_fit103.zip`
- older fit-95/98/102/103/104 posterior generations and transfer archives
- the four small `rvar_processed_*_SOB_1-4.rds` summaries

These are retained conservatively as historical/archive material rather than
active paper inputs. Do not remove them with a family-name glob: audit and
manifest any later archive cleanup separately.
