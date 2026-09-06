# Corrected policy rerun execution record

The main exercise is 75 model/cap/inference/scenario combinations (152,995 draw-level evaluations), plus cost/pooling sensitivities and short-cap median diagnostics. Existing assigned-distance fits are reused; no structural refits or manuscript edits are part of this run.

## Inputs and execution

- Census: staged `takeup_census.RData`, 39,301 unique adults across the canonical 144 communities.
- Benchmark: assigned-distance slim chains from `build/structural-fit/assigned`; 400 draws per chain, retaining the exact existing 50-per-chain comparison subset. The older 1,600-row cost input is a different fit and is not reused.
- Weighted modes: `/project/akaring/takeup-data/candidate-hpc-cd5f295-assigned/work/cluster-weight/modes`, with 999 complete assigned-distance exponential statuses.
- Common runtime: Gurobi 9.5.2 staged under `/project/akaring/takeup-data/scratch/policy-adult-population-20260906/runtime/gurobi952`, verified on both clusters.
- Jobs use eight draw workers with one solver thread each, seed zero, initially at most eight active jobs per cluster. The controller records exact submitted scripts, job IDs, and observations in its JSON ledger. Only pending jobs can be migrated, after confirmed cancellation.

## Superseded runs: geography mismatch

The shared default `/project/akaring/takeup-data/optim/data/full-many-pots-experiment.rds` contains 1,092 candidate sites, despite its canonical-looking filename. The first smoke run (`smoke-v1`) completed, and its assignments were internally consistent with that input, but its audit did not check the required candidate-site inventory. The subsequent cost smoke's explicit 1,451-site assertion caught the mismatch.

`production-v1` was stopped and its ledger-owned jobs were cancelled on both clusters. Its results must not be used for the requested exercise. The `smoke-v1-gate.json` gate was revoked. This is an input-validation failure, not an optimizer equivalence result.

The correct local `optim/data/full-many-pots-experiment.rds` has `candidate_site_mode="all"`, 1,451 unique candidate sites, and 208,944 community/site edges. It was staged separately at `/project/akaring/takeup-data/scratch/policy-adult-population-20260906/inputs/full-many-pots-experiment-1451.rds`. The adult loader now rejects a missing or incomplete geography; the independent audit also checks that the feasible-edge set includes every source edge within the cap.

The replacement `code-v3` snapshot has revision `6c94c44d1687f5e9809335955a185ecc3f61aa50e81472e1c9a4883e1155efb1`. All 117 replacement smoke tasks completed; the independent audit passed all 130 draw/scenario results, including source-edge completeness, population/target joins, achieved take-up, infeasibility classifications, and integer-optimality bounds. A Gurobi objective residue of approximately 0.000000615 required an objective comparison tolerance of 0.00001; the independently checked integer bound rules out one fewer site.

Corrected production is controlled by `temp-data/policy-adult-rerun/production-v3.json`, with up to 16 active jobs per cluster. Its output root is `/project/akaring/takeup-data/optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/policy-adult-population-20260906/production-v3`. The gate and independent audit evidence are saved beside the ledger. The controller supports `tune --ledger=... --max-active=N` without stopping its watch process, and `status` is read-only and does not acquire its writer lock.

After 107 tasks completed without failures and both partitions reported more than 1,700 idle CPUs, the limit was raised to 32 active eight-core jobs per cluster. Scheduler snapshots showed prompt starts and only brief resource waits. The reporting DAG is tracked separately in `temp-data/policy-adult-rerun/reports-v6.json`, using snapshot `code-report-v6` (manifest revision `7ed46979ae3f149b3d0e1d5fff2e266821a31e8e79e951880ec2851aeec96a48`). Its full-draw median diagnostics and figures completed and the maps and community-weighted distance density were visually reviewed. Cost jobs depend on canonical allocation collection; report assembly depends on both cost jobs, median diagnostics, and every main scenario summary.

Reporting smoke evidence: job 48979953 completed the two-draw benchmark cost reuse, equal-weight comparison, and pooling tests before the median adapter failed on missing draw IDs; those IDs were fixed. Job 48979962 then completed all five median-cap diagnostics, but rendering encountered an unavailable spatial library. The renderer was changed to use plain stored coordinates; job 48979978 rendered all figures successfully. Job 57900015 on Midway3 completed the two-draw bootstrap cost reuse and 0/50/100% pooling tests. The corrected renderer is staged separately as `render-paper-figures-v5.R`; the reporting worker and median fix are in `code-report-v4`. Final reporting should stage these tested changes together, with a fresh manifest, rather than run the older renderer embedded in `code-v3`.

## Local regression evidence

The cost-accounting extraction was compared with the committed pre-extraction solver on two three-community fixtures, with and without pooling and with nonzero signal/travel prices. Both returned identical complete allocation data frames and cost summaries. Full output validation and cost reuse tests subsequently passed on the corrected geography.

The accounting regression is reproducible with `Rscript scripts/checks/check-policy-cost-accounting.R`, with results in `temp-data/policy-adult-rerun/accounting-equivalence.csv`. Six controller tests cover missing observations, nonzero exit codes, confirmed cancellation before migration, the pending-to-running cancellation race, lost-submission reconciliation by owned job name, and enforcing the other cluster as the migration destination; run `python3 scripts/checks/test_policy_rerun_controller.py`.

## Production recovery notes

Finite-mixture Control draw 1327 initially missed its target by 0.0010118 expected adults after extracting the near-binary Gurobi solution. Snapshot `code-precision-v7` adds `IntFeasTol=1e-9` and `FeasibilityTol=1e-9`; its revision is `28f7bc181766f4f55ed09f28542eca103ad34b1beb773075854d707c2b0cdad0`. Smoke job 57900681 proved the same 103-site optimum with objective/bound 103 and achieved welfare 12778.12 versus target 12765.61. The attempted pending-job cancellation raced with successful completion; no duplicate smoke job was submitted. Retry provenance is recorded per task, and the target-attainment tolerance was not relaxed.

At 416 completed tasks, the project reached its file-count quota: `df -i` reported 887,920 used file slots and zero free, despite 916 GB of free byte capacity. Twenty-one tasks were then classified as failed, mainly late 10 km shards and collection stages. All allocation writers were confirmed terminal before archiving logs. `archive-policy-solver-logs.py` creates a tar.gz archive containing a SHA-256 manifest, verifies every archived log, checks that each source is unchanged, and only then removes individual log files. The initial archive attempt on the project failed because even one new file exceeded quota; recovery uses `/tmp/policy-solver-logs-v3-edjee-20260906.tar.gz` on Midway2 first. It was copied back and SHA-256 verified after file slots were reclaimed. Original assignment RDS files remain preserved.

The input manifest has hashed 1,141 unique source/derived inputs and recorded all 30,599 model/cap draw selections. The completed 1,600-draw benchmark cost run reports 25 median paired sites saved under corrected adult weights (central 95% interval 15.975–31), versus 27 under matched equal-community weights (16–35). Median corrected break-even site costs are 103.9633 without travel costs and 192.6256 at 0.10 per participant round-trip kilometre. Historical 260.5801/482.5247 thresholds use incorrect weights and older fit/draw provenance; this before/after change is not solely a weighting effect. `write-policy-weight-comparison.R` produces explicit matched-versus-historical comparison rows.

## Final completion evidence

All 445 production tasks, four reporting stages, 15 independent panel audits, and both assignment-level cost audits completed successfully. The assembled audit verifies 152,995 unique scenario/draw outcomes: 123,578 target-preserving optima, 29,385 target-infeasible outcomes, and 32 undefined equilibria. All 1,600 benchmark parameter rows are byte-identical across caps; adult target fields are identical. The unused community-welfare diagnostic differs by at most 1.07e-13 across hardware executions.

The verified log archives contain 123,430 pre-retry logs and 11,570 post-retry logs. The first archive SHA-256 is `22efaec4676a21480397f8e7d12c03643ae3c718c1762e390dfe3aeed527fa75`. The complete numerical archive, produced by job 48980248, has SHA-256 `e1e28caca5734a9e5a295593e5af39a075805d7ed9f093c41f2b02a2d9c94754`; its local transfer was checked against the server checksum. The final input manifest covers 1,144 unique files and all 30,599 model/cap draw selections.

The return directory is `temp-data/policy-adult-rerun/return/`; `review.pdf` compiles without layout warnings and its model/cap tables and figures were visually reviewed. The exact executed code snapshots, literal run ledgers/scripts, historical reference archive, full numerical archive, and independent audit reports are included. See `policy-adult-population-rerun-results-2026-09-06.md` for the substantive results and artifact index.
