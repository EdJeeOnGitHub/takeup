# Compact multiplier prior grid

## Purpose

This workflow answers the request for broad prior sensitivity of the level of
the structural social multiplier. It is deliberately separate from the
exponential community-weight posterior-mode bootstrap: the prior grid changes
priors under the ordinary posterior, while the bootstrap changes likelihood
weights.

## Prespecified grid

The primary estimand is the finite multiplier over 0--2.5 km. The paper-facing
table reports Control, Ink, Calendar, and Bracelet levels with 95% equal-tailed
credible intervals, plus `Pr(M_No Signal > M_Signal)` and
`Pr(M_Calendar > M_Bracelet)`.

There are 13 specifications:

1. Baseline.
2. Tight and diffuse social-image-weight priors.
3. Tight and diffuse direct-distance-cost priors.
4. Tight and diffuse idiosyncratic-heterogeneity priors, holding the
   inverse-gamma prior mean fixed.
5. Tight and diffuse first-order visibility-schedule priors.
6. Tight and diffuse private-utility/WTP priors.
7. Joint tight and joint diffuse stress tests.

The exact numerical settings and seeds are written before estimation to
`prior-grid-manifest.csv` by
`scratch/generate-main-core-prior-grid-manifest.R`.

## Estimation and audit

The initial run uses four independent chains with 400 warmup and 400 retained
draws per specification, `adapt_delta=0.999`, maximum treedepth 12, and eight
threads per chain. Chains use the common-lambda posterior mode only as a
numerically stable initial point; independent HMC seeds determine their
trajectories. Compact generated quantities retain the three finite spans
0--2.5 km, 0--1.5 km, and 1.5--2.5 km without saving ancillary latent arrays.

`scratch/summarize-main-core-prior-grid.R` requires all four fit and GQ files,
saves parameter and estimand R-hat/ESS diagnostics, divergences, and treedepth
hits, and writes `prior-grid-needs-rerun.csv`. A specification is flagged for
an 800/800 rerun if R-hat exceeds 1.01, bulk or tail ESS is below 400, or any
divergence or treedepth hit occurs.

After the initial summary completes,
`submit_main_core_prior_grid_reruns.sh` reads that audit and submits only the
flagged specifications at 800/800. Rerun fits and GQs are isolated under
`fits-rerun/` and `gq-rerun/`; the final summary automatically prefers a
complete rerun pair and records the source used for every row.

## Midway3 production run

The isolated code staging directory is:

`/project/akaring/takeup-data/code/main-core-prior-grid`

The shared output directory is:

`/project/akaring/takeup-data/data/stan_analysis_data/main-core-prior-grid`

Submitted on 2026-08-10:

- prepare: 53216644
- compile: 53216645
- sample array: 53216646 (52 chain jobs, at most 16 concurrent)
- compact GQ array: 53216647 (13 jobs, at most 8 concurrent)
- summary: 53216648

The compile completed successfully. The first sample array failed before Stan
started because the R 4.2 linear-algebra library raised an illegal-instruction
fault in `eigen()` on Midway3's AMD partition. The unused dependent GQ and
summary jobs were cancelled. Sampling was resumed on the known-compatible
Caslake partition without recompiling:

- sample array: 53260438
- compact GQ array: 53260443
- summary: 53260444

All fit and compact-GQ tasks completed. The first summary attempt exposed a
vector-length bug in the directional-probability grouping code; no estimates
were affected. The condition was corrected to use the group's scalar estimand
label and the summary was rerun from the existing 52 GQ files.

The shared data helper was then hardened so Gaussian and estimated finite-
mixture fits no longer invoke the unused Student-t quadrature eigensolver.
Student-t fits still compute their required quadrature exactly. This removes
the AMD-specific pre-launch failure for future ordinary core-model runs.

The initial audit flagged ten specifications under the prespecified strict
whole-parameter rule. Their multiplier estimands were already well mixed
(maximum estimand R-hat 1.008; minimum estimand bulk ESS 643), but the longer
reruns were submitted as planned:

- 800/800 sample array: 53264536
- rerun compact GQ array: 53264537
- final combined summary: 53264545

After 800/800, only the joint-diffuse stress test remained flagged: 4 of 3,200
transitions (0.125%) reached treedepth 12. It had no divergences, maximum
parameter R-hat 1.004, minimum parameter tail ESS 905, maximum multiplier
R-hat 1.002, and minimum multiplier tail ESS 2,096. We nevertheless preserve
the strict rule and rerun only this row at 800/800 with maximum treedepth 15.
Those outputs live separately under `fits-rerun-td15/` and
`gq-rerun-td15/`; the summarizer prefers them only when all four fit and GQ
files are present.

Final targeted jobs:

- treedepth-15 sample array: 53295444
- treedepth-15 compact GQ: 53295445
- final summary: 53295446

All 13 specifications pass the final automated audit.

No files in the main manuscript or Overleaf tree are changed by this workflow.
