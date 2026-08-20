# Minimal main model: cluster heterogeneity and cluster-weighted inference

## Implemented design

The minimal main model now supports two deliberately separate exercises.

1. `use_core_cluster_shock=1` adds the non-centred additive utility shock
   `core_cluster_shock_sd * core_cluster_shock_raw[c]` before the community
   fixed-point solve. The primary half-normal scale is 0.1 and the sensitivity
   scale is 0.25. Counterfactual generated quantities hold the realized shock
   fixed across treatment and distance.
2. `core_cluster_weight[c]` weights every cluster-linked contribution to the
   joint likelihood: monitored take-up, beliefs, randomized distance, and WTP.
   Priors are never weighted. Unit weights and a disabled shock recover the
   baseline target.

The data builder now retains `wtp_cluster_id`. Old workspaces can still run
unweighted fits because the WTP mapping is irrelevant with unit weights, but
weighted refits reject an old workspace rather than silently reweight WTP
incorrectly.

## Midway workflow

From `~/projects/takeup`:

```bash
# Submit data preparation, the unweighted initializer, cluster-weighted
# refits, and paper-output generation with dependencies.
bash hpc/structural/submit_main_core_cluster_bootstrap.sh
```

`hpc/structural/slurm_main_core_cluster_bootstrap.sh` is the only production worker and uses
the `STAGE` variable internally. `hpc/structural/submit_main_core_cluster_bootstrap.sh` is the
small submission wrapper. The generic `hpc/structural/slurm_main_core.sh` remains available
for a full posterior run of the minimal main model when needed.

Optimizer JSON files are sanitized before HMC warm starts. In particular, an
exact optimizer solution on a nonnegative boundary is moved to `1e-4` while
singleton Stan vectors retain their JSON array shape. This affects only the
initial state; it does not constrain or modify the sampled posterior.

## Validation gates

- Both sampling and compact-GQ Stan programs must pass full C++ compilation.
- The existing frozen-generic versus minimal analytic-gradient validation must
  pass with unit weights and the shock disabled.
- Every generated multinomial weight vector must preserve its county cluster
  count; weights may be zero because omitted clusters are part of ordinary
  resampling with replacement.
- No weighted production array should start until the regenerated workspace
  contains one valid `wtp_cluster_id` for every WTP row and the unweighted mode
  has completed.

## Completed production checks

- The regenerated main workspace contains 144 clusters and valid cluster IDs
  for all 998 WTP observations.
- Unit-weight/no-shock validation on the regenerated workspace gave maximum
  analytic-gradient differences of exactly zero at all three validation
  points. Direct Bernoulli/binomial kernel residuals were at most
  `9.09e-13`.
- The primary cluster-shock run used four chains with 1,000 warmup and 1,000
  retained draws. It had zero divergences and zero maximum-treedepth hits. The
  cluster-shock SD posterior had median 0.347, 5th--95th percentiles
  0.304--0.396, R-hat 1.00, and bulk ESS 1,810.
- The preferred paper-facing workflow uses 999 successful exponential
  cluster-weighted refits. All 999 attempted refits completed, all 999 weight
  vectors were distinct, and no refit was excluded. Conventional multinomial
  resampling remains an archived internal sensitivity check.
- All 16 short-MCMC audit cases completed. The combined audit passed the
  prespecified threshold: 289 of 304 reported cells (95.1%) differed from the
  short-MCMC posterior median by less than 10% of the corresponding refit
  interval width. The method-specific rates were 93.4% for the multinomial
  bootstrap and 96.7% for exponential weights; this slight bootstrap-specific
  shortfall should be disclosed when interpreting the mode-based percentile
  intervals.

## Paper-facing outputs

- `presentations/tables/fit105/main-core-cluster-robustness.tex` compares the
  social multipliers, with the exponential weighted-likelihood bootstrap
  immediately after baseline.
- `presentations/tables/fit105/main-core-exponential-cluster-weight-overall-te-table.tex`
  reports the weighted structural take-up results in the same Combined,
  Close, Far, and Far-minus-Close layout as the main structural table.
- `presentations/misc-figures/main-core-cluster-robustness.pdf` is the matching
  multiplier figure.
- Complete summaries and the audit cells are stored under
  `temp-data/main-core-cluster-robustness/` locally and under the corresponding
  `main-core-cluster-robustness` analysis directory on Midway.

## Output interpretation

- Baseline and cluster-shock intervals are posterior credible intervals.
- Multinomial-weight intervals are cluster-bootstrap percentile intervals.
- Exponential-weight intervals are a generalized-Bayes sensitivity exercise.
- The cluster shock changes the structural model; weighted refits adjust the
  baseline model's sampling uncertainty. They should not be described as the
  same robustness check.
- The exponential cluster weighted-likelihood bootstrap is the preferred
  generalized-Bayesian cluster adjustment. Conventional multinomial
  resampling is retained internally but is not promoted in the paper-facing
  appendix. The cluster random intercept remains a separate structural
  robustness specification.
