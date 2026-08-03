# Minimal main structural model: target validation and benchmark

Date: 2026-08-02

## Target-equivalence tests

Test script: `scratch/validate-main-core-target.R`

Midway validation job: `48723593`

The frozen fit-105 generic model and `takeup_struct_main_core.stan` were
evaluated on the same 139-community no-outlier fit-105 data and at the same
three unconstrained parameter points (zero plus two small random
initializations).

- The maximum absolute difference between the complete models' analytic
  gradients was exactly zero at all three points (with 18-digit diagnostic
  output).
- The only intentional likelihood rewrite is aggregation of individual
  Bernoulli terms into the unnormalized cluster-binomial kernel using
  `binomial_lupmf`. This explicitly omits the normalized binomial likelihood's
  data-only constant, `5368.15976983222`, and gives the same Stan target as the
  repeated Bernoulli terms.
- Direct kernel checks at constant, deterministic heterogeneous, and random
  heterogeneous cluster probabilities had residual exactly zero in double
  precision.

CmdStan 2.33 rounds the diagnostic `Log probability=` summary line to six
significant figures even when `output sig_figs=18` is requested. The automated
test therefore uses the full-precision analytic gradients for the complete
model and a separate direct evaluation of the Bernoulli/binomial kernel.
Together these confirm that the unnormalized Stan sampling target is unchanged.

## Matched Midway benchmark

Benchmark script: `slurm_main_core_benchmark.sh`

Midway benchmark job: `48723540`

Both models ran sequentially on the same node with the same data,
initialization, seed, one CPU, diagonal metric, `adapt_delta=0.9999`, maximum
treedepth 12, 100 warmup iterations, and 100 sampling iterations. CmdStan's
reported sampling times exclude compilation.

| Model | Warmup (s) | Sampling (s) | Total (s) | CSV size | Divergences | Max-depth transitions |
|---|---:|---:|---:|---:|---:|---:|
| Frozen generic main model | 215.172 | 854.322 | 1069.490 | 426,896,694 bytes | 0 | 0 |
| Minimal main model | 129.088 | 546.139 | 675.227 | 36,421 bytes | 0 | 0 |

The minimal model was:

- `1.584x` faster overall (`36.9%` less total Stan time),
- `1.667x` faster in warmup,
- `1.564x` faster in retained sampling, and
- `11,721x` smaller on disk (`99.9915%` reduction).

The complete sequential Slurm job took 38 minutes 15 seconds and peaked at
4.22 GB RSS, but that peak combines both compilation processes and both model
runs; it is not a clean per-model memory comparison.

## Historical main-fit baseline

The four existing fit-105 main-model chains each used 400 warmup and 400
sampling iterations on one thread. Their total reported times were 985.417,
955.151, 875.807, and 1223.820 seconds (mean 1009.999 seconds), and each CSV was
approximately 1.73 GB. Those historical runs used a different initialization
and compute environment, so the matched benchmark above is the appropriate
speed comparison. Applying its speed ratio mechanically would imply roughly
10.6 minutes per historical-length chain, but a production pilot is needed
before treating that projection as a scheduling estimate.

## Cluster-extension regression validation

Midway job `48723660` reran the complete target-equivalence test after adding
the optional cluster shock and joint-likelihood cluster weights. The test used
unit weights, disabled the shock, and supplied a harmless placeholder WTP
cluster mapping for the legacy workspace. The maximum absolute analytic-gradient
difference from the frozen generic model remained exactly zero at all three
test points. The three direct Bernoulli/binomial kernel residuals also remained
exactly zero. Thus the new interface does not alter the baseline target.

The same validation was repeated with the regenerated 144-cluster production
workspace in Midway job `48724944`. Maximum analytic-gradient differences were
again exactly zero at all three points. The direct likelihood residuals were
`9.09e-13`, zero, and zero, confirming target equivalence with the production
WTP-to-cluster mapping as well as with the legacy compatibility mapping.
