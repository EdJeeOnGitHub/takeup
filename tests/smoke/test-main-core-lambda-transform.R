#!/usr/bin/env Rscript

source("R/structural/main-core-data.R")

helmert_basis <- matrix(
  c(
    1 / sqrt(2), 1 / sqrt(6), 1 / sqrt(12),
    -1 / sqrt(2), 1 / sqrt(6), 1 / sqrt(12),
    0, -2 / sqrt(6), 1 / sqrt(12),
    0, 0, -3 / sqrt(12)
  ),
  nrow = 4,
  byrow = TRUE
)
stopifnot(
  max(abs(colSums(helmert_basis))) < 1e-12,
  max(abs(crossprod(helmert_basis) - diag(3))) < 1e-12
)

base_lambda <- 0.8
theta <- 0.4
group_lambda <- base_lambda * exp(c(-theta / 2, theta / 2, -theta / 2, theta / 2))
stopifnot(
  abs(exp(mean(log(group_lambda))) - base_lambda) < 1e-12,
  abs(mean(log(group_lambda[c(2, 4)])) -
      mean(log(group_lambda[c(1, 3)])) - theta) < 1e-12
)

set.seed(20260810)
prior_sd <- 0.25
raw <- matrix(rnorm(3 * 200000), ncol = 3)
deviation <- prior_sd / sqrt(2) * raw %*% t(helmert_basis)
stopifnot(
  max(abs(rowSums(deviation))) < 1e-12,
  abs(sd(deviation[, 1] - deviation[, 2]) - prior_sd) < 0.002,
  abs(sd(deviation[, 2] - deviation[, 4]) - prior_sd) < 0.002
)

# The deterministic normal-mixture rule should closely reproduce a
# variance-one t_5 CDF/density, while retaining the Gaussian decision-noise
# convolution needed by the structural type-inference term.
mixture <- main_core_student_t_mixture(df = 5, components = 12)
evaluation_grid <- seq(-4, 4, by = 0.05)
mixture_cdf <- vapply(evaluation_grid, function(cutoff) {
  conditional_sd <- sqrt(mixture$scale_sq / mixture$precision)
  sum(mixture$weight * pnorm(cutoff / conditional_sd))
}, numeric(1))
mixture_pdf <- vapply(evaluation_grid, function(cutoff) {
  conditional_sd <- sqrt(mixture$scale_sq / mixture$precision)
  sum(mixture$weight * dnorm(cutoff / conditional_sd) / conditional_sd)
}, numeric(1))
exact_cdf <- pt(evaluation_grid / sqrt(mixture$scale_sq), df = mixture$df)
exact_pdf <- dt(evaluation_grid / sqrt(mixture$scale_sq), df = mixture$df) /
  sqrt(mixture$scale_sq)
stopifnot(
  max(abs(mixture_cdf - exact_cdf)) < 3e-4,
  max(abs(mixture_pdf - exact_pdf)) < 2e-4,
  abs(sum(mixture$weight) - 1) < 1e-12,
  all(mixture$weight > 0),
  all(mixture$precision > 0)
)
cat("Main-core lambda transforms: OK\n")
