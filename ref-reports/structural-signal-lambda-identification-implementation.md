# Signal-specific lambda identification workflow

This workflow implements the common, Any-Signal/No-Signal, and arm-specific
social-image-weight specifications requested by the structural robustness
plan. The flexible specifications use a common geometric mean and a prior SD
defined on pairwise log-lambda ratios.

## Production workflow

On Midway, run:

```bash
bash submit_main_core_lambda_identification.sh core
```

The submission helper reuses the four completed common-lambda production
chains, fits the grouped and arm-specific prior grid, computes compact
generated quantities, evaluates a regularized profile on
`log(lambda_signal/lambda_no_signal)`, and runs the hybrid recovery design.
Arrays are throttled to 20 profile jobs by default. Because Midway counts every
pending array element against the per-user submission limit, recovery is
submitted only after the core/profile queue drains, in these phases:

```bash
bash submit_main_core_lambda_identification.sh recovery-generate
bash submit_main_core_lambda_identification.sh recovery-fit-a
# after fit-a completes
bash submit_main_core_lambda_identification.sh recovery-fit-b
# after fit-b completes
bash submit_main_core_lambda_identification.sh recovery-hmc
# after the HMC audit completes
bash submit_main_core_lambda_identification.sh recovery-finish
```

Recovery modes and HMC chains are throttled to 40 and 20 concurrent jobs.

All raw and summarized results live below
`main-core-lambda-identification` in the shared Stan analysis directory. The
tracked assembly script copies completed paper artifacts into
`presentations/structural-lambda-identification` without modifying the main
paper file.

## Interpretation guardrails

- The profile retains baseline priors on nuisance parameters but removes the
  prior on the fixed signal/no-signal contrast. It is therefore labeled a
  regularized profile, not a formal likelihood-ratio confidence interval.
- Posterior existence is not treated as identification. Summaries report
  posterior/prior contraction and correlations with private utility and
  distance-cost parameters.
- Failed convergence, flat profiles, or poor simulation recovery remain in
  the manifests and are reported rather than silently discarded.
- Recovery uses the coherent posterior draw nearest the componentwise median
  of the structural parameters. This preserves simplex and other joint
  constraints that a literal vector of componentwise medians can violate.

## Validation

- `scratch/validate-main-core-target.R` compares the common-lambda minimal
  model with the original fit-105 target at zero and random parameter points.
  The maximum target-gradient difference is zero; the Bernoulli-to-binomial
  aggregation identity also has zero residual.
- `scratch/test-main-core-lambda-transform.R` checks the grouped geometric-mean
  parameterization and the arm-specific Helmert basis, including its stated
  pairwise log-ratio prior standard deviation.
- All three Stan programs (sampler, compact generated quantities, and hybrid
  simulator) compile under CmdStan 2.36. The recovery path is additionally
  gated on one complete optimize-plus-generated-quantities fit before either
  large mode array is submitted.
