# Asymmetric observability and noisy social inference

## Objective

The baseline structural model treats the first-order definite-answer rate as
the probability that a peer's action is revealed perfectly. This robustness
exercise instead estimates how true administrative take-up is mapped into a
respondent's report (Yes, No, or Don't Know), and lets agents and observers
internalize that noisy information technology in equilibrium.

This is the first stage of the broader observability response. It addresses
asymmetric false positives and false negatives before Lee-bound and MAR/MNAR
selection corrections. The latter remain separate because missing outcomes do
not identify outcome-dependent selection without an additional assumption.

## Measurement model

Use the exact Table-A and no-SMS sample rules used by the main structural
belief likelihood. Link each named peer to administrative take-up. Peers whose
true status cannot be linked are excluded from the response-matrix likelihood
and retained in an explicit coverage audit.

For peer `j` evaluated by respondent `i`, estimate a sequential model

\[
 R_{ij}\sim\operatorname{Bernoulli}(r_{yzd}),\qquad
 S_{ij}\mid R_{ij}=1,y_j\sim\operatorname{Categorical}(q_{yzd}),
\]

where `R` denotes recognition, `S` is Yes/No/Don't Know, `y` is true
administrative take-up, `z` is treatment, and `d` is standardized community
distance. Recognition and conditional report logits have truth-specific
global intercepts and slopes, sum-to-zero treatment deviations, and
treatment-specific distance deviations. Priors are `N(0,1.5)` for global
intercepts, `N(0,0.5)` for global distance slopes and treatment deviations,
and `N(0,0.25)` for treatment-specific slope deviations.

The primary information matrix conditions on recognition and contains the
three report categories. The secondary matrix is unconditional and appends
unrecognized as a fourth null signal:

\[
 Q^U_{yzd}=\{r q_{Yes},r q_{No},r q_{DK},1-r\}.
\]

The exercise identifies the response matrix conditional on true status. It
does not claim to distinguish genuine knowledge from lucky guessing.

## Equilibrium mapping

For take-up probability `pi` and signal probabilities `q_sy`, replace the
baseline `p_observed * Delta(w*)` with

\[
 A(\pi,Q)\Delta(w^*),\qquad
 A(\pi,Q)=\pi(1-\pi)\sum_s
 \frac{(q_{s1}-q_{s0})^2}
 {\pi q_{s1}+(1-\pi)q_{s0}}.
\]

The fixed point remains scalar:

\[
 w^*+b_z-cd+\lambda A(\pi(w^*),Q_{zd})\Delta(w^*)=0.
\]

This nests the current model exactly. If the action is revealed correctly
with probability `p` and otherwise produces a common null signal, then
`A(pi,Q)=p` for every `pi`.

Generated quantities use the implicit derivative

\[
 F_w=1+\lambda(A_w\Delta+A\Delta_w),\qquad
 F_d=-\beta_d+\lambda A_d\Delta
\]

to calculate the total distance response and social multiplier. Derivatives
of the softmax response probabilities and `A` are analytic and checked against
central finite differences.

## Model modes and outputs

The minimal main model exposes one switch:

1. existing perfect-observation baseline;
2. asymmetric reports conditional on recognition (primary);
3. asymmetric reports including unrecognized as a null signal (sensitivity).

The first production comparison keeps the common social-image weight and
Gaussian intrinsic type distribution fixed so that only the information
technology changes. Compact generated quantities report response matrices,
sensitivity, specificity, recognition, the information factor, image returns,
fixed points, structural treatment effects, and multipliers.

Paper-facing output consists of response-matrix diagnostics, the existing
Combined/Close/Far/Far-minus-Close ATE layout, multipliers at 0.5, 1.5, and
2.5 km, and fast sparse policy results for dynamic Control and Bracelet,
their 0.5 km fixed-information counterparts, and the no-image case.

## Validation gates

- Reproduce the existing respondent-level recognition and definite-answer
  counts from the peer records; audit administrative-link coverage by arm and
  randomized-distance cell.
- Unit-test perfect, uninformative, symmetric-error, and null-signal response
  matrices. Every matrix must be finite and row-stochastic.
- With the new mode disabled, reproduce the frozen minimal likelihood and
  gradient exactly. With the restricted null-signal matrix, reproduce the old
  fixed point, multiplier, and policy demand.
- Compare analytic information-factor and multiplier derivatives with central
  finite differences at baseline and random parameter points.
- Run recovery designs with perfect, symmetric-error, asymmetric-error, and
  uninformative signals; an uninformative signal must not generate an image
  multiplier.
- Production begins with four chains of 400 warmup and 400 retained draws. If
  any scalar has split R-hat above 1.01, bulk or tail ESS below 400, any
  divergence, or any maximum-treedepth transition, rerun with 1,000 warmup and
  1,000 retained draws.

## Deferred selection stage

After asymmetric inference is validated, propagate the existing pairwise Lee
bounds through the structural model, add covariate-standardized MAR schedules,
and report an MNAR pattern-mixture tipping-point grid. The combined adverse
specification will apply the least favorable supported missingness schedule to
the unconditional asymmetric response matrix.

## Implementation status (2026-08-04)

Implemented in the minimal main-model sampler, compact generated quantities,
and fast sparse policy adapter. The linked-data audit recovers all 1,098 frozen
FOB/drop respondents and all 10,980 Table-A records; 10,079 records (91.8
percent) link to the named peer's administrative take-up. Link coverage is
written by treatment and randomized-distance cell for every asymmetric run.

The disabled switch has zero-dimensional measurement parameters and retains
the old likelihood code path. Direct CmdStan gradient checks against the
frozen generic sampler agree exactly at the zero initialization and two seeded
random initializations. Both asymmetric modes pass data initialization and
gradient evaluation; the conditional mode also passes an HMC and compact-GQ
smoke run. The unit test is
`scratch/test-main-core-asymmetric-observability.R`.

Production commands use `CORE_OBSERVATION_MODEL=1` (primary) or `2`
(unrecognized null signal) with `slurm_main_core.sh`. Start at 400 warmup and
400 retained draws in four array tasks, then run
`scratch/generate-main-core-compact-gq.R` with the matching observation-model
switch. The summary script
`scratch/summarize-main-core-asymmetric-observability.R` produces diagnostics,
response-channel, ATE, and multiplier CSVs plus the two panelled LaTeX tables.
The production fits themselves are intentionally not marked complete until all
four chains pass the stated diagnostic gates.
