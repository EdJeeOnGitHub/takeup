# Server handoff: full-sample two-stage accuracy gradients

## Objective

Estimate whether allowing reporting accuracy to change with distance improves
fit, and assess sensitivity to pooling accuracy gradients. Do not select a model
because it produces preferred multipliers. No manuscript replacement is authorized
by this handoff; return results for review.

Use ALL 144 communities, 9,805 take-up observations and 1,141 belief respondents.
Do not use NO_OUTLIERS or the historical main-core-asym-input workspace. The old
F3/U3 runs excluded five dispersed communities. Local verification of the full
workspace recovers 11,410 peer rows, 10,456 linked rows, 4,962 recognized linked
rows and 3,926 definite linked rows. See README.md for the measurement-only audit.

## Exact specification

Primary comparison conditions on recognition: observation_model=1,
recognition_structure=2. Keep the Gaussian type distribution, common image weight,
private-payoff/WTP priors, direct distance cost and all other behavioral assumptions
fixed across the new two-stage runs. Do not widen the valuation prior in this
exercise. Use the main model STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP, realized
community-centroid distance, and the same workspace for every run.

For true status y (1=nonparticipant, 2=participant), treatment z, and existing
standardized distance x:

    logit D_yz(x) = alpha_y + C_z b_y + (g_y + h * PublicSignal_z) x

This definite-answer equation is UNCHANGED from F3. PublicSignal=Ink/Bracelet;
h is shared across truth states, exactly as before. C is the existing orthonormal
sum-to-zero treatment contrast basis, canonical order Control, Ink, Calendar,
Bracelet. Do not replace it with reference-cell dummies.

Conditional accuracy is:

    mode 2: logit A_yz(x) = a_y + C_z c_y
    mode 3: logit A_yz(x) = a_y + C_z c_y + k_y x
    mode 4: logit A_yz(x) = a_y + C_z c_y + (k_y + tau_y C_z u_y) x

Priors added in code:

    k_y ~ Normal(0, 0.5)
    u_yj ~ Normal(0, 1)
    tau_y ~ HalfNormal(0, 0.25)

The two tau values are separate by truth; each scales three orthonormal arm
contrasts. These priors use the EXISTING standardized distance units, not km.
Record the actual distance scaling from the workspace in the results. Existing
accuracy intercept priors remain Normal(0,1.5) globally and Normal(0,0.5) on
arm contrasts. Existing definite-answer priors are unchanged.

The Yes/No/DK response matrix is:

    participant:    (D*A, D*(1-A), 1-D)
    nonparticipant: (D*(1-A), D*A, 1-D)

Use this matrix in the existing Bayesian noisy-information factor A(pi,Q),
fixed point and multiplier. This is NOT substituting correct-answer rates into
benchmark perfect observability. The GQ derivative includes BOTH the definite
and accuracy gradients. Mode 2 and the older derivative overload are retained.

## Files changed locally: sync all together

- stan_models/takeup_struct_main_core.stan
- stan_models/takeup_struct_main_core_compact_gq.stan
- stan_models/core_asymmetric_observability_functions.stan
- R/structural/main-core-data.R
- scripts/structural/sample-main-core.R
- tests/smoke/test-two-stage-accuracy-gradient.R
- this directory (audit and handoff)

These edits are included in the accuracy-gradient implementation commit. Sync
that commit and separately transfer the full workspace if absent on the server.
Record source hashes, workspace hash and exact commands. Do not reset unrelated changes. Force compilation of the changed sampler and compact GQ;
never reuse an executable built before these edits. Compile once before parallel
array submission to avoid concurrent writes to the shared executable.

## Preflight and compute environment

Prefer Midway3/caslake, using the same working R/CmdStan setup as the recent
valuation audit. The existing hpc/structural/slurm_main_core.sh contains old
broadwl/module defaults: adapt environment/module/CmdStan paths for Midway3 if
needed; do not alter model settings. Run commands from the server repository root.
The WORKSPACE below must be the full local workspace synced to that relative path,
or an identical full-sample server workspace verified by counts and priors.

    Rscript tests/smoke/test-main-core-asymmetric-observability.R
    Rscript tests/smoke/test-two-stage-accuracy-gradient.R

Before production, compile both complete Stan models; check sampler/GQ parameter
schemas match. Run 10 warmup + 10 retained HMC smoke draws for each new mode and
compact GQ on those draws. Require finite probabilities, fixed points and gradients,
row sums=1, and information factors in [0,1] up to rounding. Check finite-difference
TOTAL take-up distance derivatives against GQ multipliers at 0.5/1.25/2.5 km.
Local tests cover reporting derivatives; full equilibrium smoke remains server work.

## Production fits

Four chains per variant, 1,000 warmup and 1,000 retained per chain initially;
adapt_delta=0.99, max_treedepth=12, 8 threads/chain. First time a chain through
warmup to estimate remaining wall time; do not extrapolate from the tiny R audit.
Run constant accuracy on the full sample too: historical F3 is not a valid
same-sample comparator. Benchmark perfect observation is a fourth comparator.

Bash template after setting the compatible module environment in the wrapper:

```bash
export WORKSPACE="$PWD/build/structural-workspace/main-core-input.RData"
export MODEL=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
export OUTPUT_ROOT=/project/akaring/takeup-data/data/stan_analysis_data/accuracy-gradients-full-20260907
export DISTANCE_DEFINITION=realized
export CORE_LAMBDA_STRUCTURE=0 USE_CORE_CLUSTER_SHOCK=0
export CORE_REPORT_ARM_DIST_HIERARCHICAL=0
export THREADS_PER_CHAIN=8 ITER_WARMUP=1000 ITER_SAMPLING=1000
export ADAPT_DELTA=0.99 MAX_TREEDEPTH=12 SEED=20260907
mkdir -p temp/log "$OUTPUT_ROOT"
for mode in 2 3 4; do
  sbatch --partition=caslake --array=1-4 \
    --export="ALL,OUTPUT_PATH=$OUTPUT_ROOT/mode-$mode,CORE_OBSERVATION_MODEL=1,CORE_RECOGNITION_STRUCTURE=2,CORE_REPORT_STRUCTURE=$mode" \
    hpc/structural/slurm_main_core.sh
done
sbatch --partition=caslake --array=1-4 \
  --export="ALL,OUTPUT_PATH=$OUTPUT_ROOT/benchmark,CORE_OBSERVATION_MODEL=0,CORE_RECOGNITION_STRUCTURE=0,CORE_REPORT_STRUCTURE=0" \
  hpc/structural/slurm_main_core.sh
```

Set CMDSTAN_PATH to the actual supported server version before these commands.
Ensure no inherited INIT_FILE, CLUSTER_WEIGHT_FILE, STAN_FILE or other overrides
change the intended comparison. Explicitly record effective settings in manifests.
Do not combine chains generated with different model/workspace versions.

Generate compact GQ using scripts/structural/generate-compact-gq.R with the SAME
workspace, model, distance-definition and observation/recognition/report modes.
Pass the four actual chain CSVs via --fit-csvs as a comma-separated list. Example:

```bash
Rscript scripts/structural/generate-compact-gq.R \
  --workspace="$WORKSPACE" --model="$MODEL" \
  --fit-csvs="$FOUR_CHAIN_CSVS" --output-path="$OUTPUT_ROOT/mode-3/gq" \
  --core-observation-model=1 --core-recognition-structure=2 \
  --core-report-structure=3 --distance-definition=realized
```

Repeat matching settings for mode 2, mode 4 and benchmark. Never use mismatched
report-mode parameter schemas to generate quantities. New modes have additional
parameters; old CSVs cannot stand in for their new fits.

## Required return outputs

1. Source/workspace hashes, run times and exact settings, sample/link counts,
   priors including distance scaling and valuation conversion prior.
2. Sampling diagnostics: max Rhat, min bulk/tail ESS, divergences and treedepth.
   Target Rhat<=1.01, ESS>=400, zero divergences/treedepth hits. Extend/revise
   sampling as necessary; flag unconverged results instead of presenting as final.
3. Posterior accuracy slopes in standardized units AND per km, by truth/arm;
   common slopes, tau and within-arm deviations, intervals and correlations.
4. Observed and predicted Yes/No/DK frequencies, conditional accuracy, definite
   rates and information factors by truth, arm and distance. Include bin counts.
5. Take-up levels and ATEs by arm and Close/Far, including Bracelet-Calendar and
   pooled signal contrasts, using the SAME empirical distances and weights for
   all models. Include comparable observed moments. Legacy fitted-lognormal ATEs
   may be supplied separately but must not be compared to differently weighted RF.
6. Multiplier curves, with uncertainty, at 0.5/1.25/1.5/2.5 km and finite change
   0.5-2.5 km. Paired No-Signal minus Any-Signal contrasts within each draw.
   Report failed/nonfinite draws explicitly. Do not censor negative multipliers.
7. Separate measurement fit and take-up fit; do not select using preferred
   multipliers. Held-out-community validation is a subsequent comparison if
   the candidates pass these basic checks; in-sample fit is not held-out evidence.

Do not spend time on policy allocations yet. Return compact summaries and CSVs
into ref-reports/accuracy-gradient-audit/server-results/ and append conclusions
and any unresolved issues to this handoff. Keep raw chains on the server with
explicit paths. Do not silently change priors, drop communities, or tighten
slopes to obtain a preferred sign.

## Local validation completed

- Both full models passed stanc 2.39 syntax/type checks.
- Sampler and compact-GQ parameter declarations match exactly (ignoring comments/whitespace).
- Existing asymmetric-observability smoke tests passed.
- Compiled Stan reporting derivative agrees with central finite differences for
  both truth states; zero accuracy slope exactly nests the old derivative.
- prepare_main_core_data accepts modes 2/3/4 and confirms full-sample counts.
- git diff --check passed for edited production code.
- Full-model HMC and equilibrium-derivative validation have NOT been run locally;
  those are explicit server preflight gates above.
