# Student-t latent-type robustness

This check replaces the baseline standard-normal intrinsic-motivation type
with a variance-normalized Student-t distribution with five degrees of
freedom. The social-image inference term remains the difference in posterior
mean intrinsic motivation between takers and non-takers. Gaussian decision
noise is retained, so both the take-up probability and the inference term use
the convolution of the Student-t type and the estimated normal noise.

The Student-t is evaluated as a deterministic generalized Gauss-Laguerre
normal-variance mixture with 12 components. Over the tested range, the
approximation error is below `3e-4` for the t5 CDF and below `2e-4` for its
density. The baseline Gaussian branch is selected by a data switch and has
been re-audited against the frozen original target with zero gradient
difference.

The robustness workflow is:

```bash
bash hpc/structural/submit_main_core_student_t.sh
```

It first computes a mode, then runs four 400-warmup/400-sampling chains,
compact generated quantities, diagnostics, and a Gaussian-versus-t5 social
multiplier table. Outputs are stored in
`main-core-student-t-robustness` under the shared Stan analysis directory.

This is a parametric heavy-tail robustness check, not a nonparametric estimate
of the latent-type distribution. Degrees of freedom are fixed at five so that
the exercise changes tail thickness without introducing another weakly
identified shape parameter.

The paper-facing discussion and harmonized Gaussian-versus-Student-t table are
in `appendix/structural-robustness/sections/student-t-types.tex` and
`appendix/structural-robustness/tables/main-core-student-t-multipliers.tex`.
They form part of the portable `docmute` appendix bundle and are not yet
included in the main manuscript.
