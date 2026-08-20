#!/usr/bin/env Rscript

# Diagnose the deliberately weakly identified signal-specific-lambda
# sensitivity. In particular, report whether treatment contrasts contract
# relative to their prior and whether they trade off with private utilities.

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

fit_csvs <- strsplit(
  option_value("--fit-csvs", ""), ",", fixed = TRUE
)[[1L]]
output_path <- option_value(
  "--output-path", "temp-data/main-core-signal-lambda"
)
lambda_structure <- as.integer(option_value("--core-lambda-structure", "2"))
log_ratio_sd_prior <- as.numeric(
  option_value("--core-lambda-log-ratio-sd-prior", "0.25")
)
if (length(fit_csvs) < 2L || any(!file.exists(fit_csvs))) {
  stop("--fit-csvs must list at least two existing chain CSVs.", call. = FALSE)
}
if (!lambda_structure %in% 1:2) {
  stop("--core-lambda-structure must be 1 or 2.", call. = FALSE)
}
if (!is.finite(log_ratio_sd_prior) || log_ratio_sd_prior <= 0) {
  stop("The log-lambda ratio prior SD must be positive.", call. = FALSE)
}

fit <- read_cmdstan_csv(fit_csvs)
raw_names <- if (lambda_structure == 1L) {
  "core_lambda_group_log_ratio_raw[1]"
} else {
  paste0("core_lambda_arm_log_ratio_raw[", 1:3, "]")
}
required <- c("base_mu_rep", raw_names)
available <- dimnames(fit$post_warmup_draws)[[3]]
if (!all(required %in% available)) {
  stop("Fit does not contain the signal-specific-lambda parameters.", call. = FALSE)
}

draws <- as_draws_df(fit$post_warmup_draws)
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
treatment <- c("Control", "Ink", "Calendar", "Bracelet")
raw <- as.matrix(draws[, raw_names, drop = FALSE])
if (lambda_structure == 1L) {
  group_log_ratio <- log_ratio_sd_prior * raw[, 1]
  log_deviation <- outer(
    group_log_ratio,
    c(-0.5, 0.5, -0.5, 0.5)
  )
  marginal_log_prior_sd <- log_ratio_sd_prior / 2
} else {
  log_deviation <-
    log_ratio_sd_prior / sqrt(2) * raw %*% t(helmert_basis)
  marginal_log_prior_sd <- log_ratio_sd_prior * sqrt(3 / 8)
}
lambda <- exp(log(draws$base_mu_rep) + log_deviation)
signal_vs_no_signal <- rowMeans(log(lambda[, c(2, 4), drop = FALSE])) -
  rowMeans(log(lambda[, c(1, 3), drop = FALSE]))

summarize_columns <- function(values, labels, quantity, prior_sd = NA_real_) {
  data.frame(
    treatment = labels,
    quantity = quantity,
    mean = apply(values, 2, mean),
    sd = apply(values, 2, sd),
    q025 = apply(values, 2, quantile, probs = 0.025),
    q50 = apply(values, 2, quantile, probs = 0.5),
    q975 = apply(values, 2, quantile, probs = 0.975),
    prior_sd = prior_sd,
    posterior_to_prior_sd = if (is.na(prior_sd)) {
      NA_real_
    } else {
      apply(values, 2, sd) / prior_sd
    }
  )
}

lambda_summary <- rbind(
  summarize_columns(lambda, treatment, "lambda"),
  summarize_columns(
    log_deviation,
    treatment,
    "log_lambda_deviation",
    marginal_log_prior_sd
  )
)
contrast_summary <- data.frame(
  contrast = "Any Signal - No Signal",
  mean = mean(signal_vs_no_signal),
  sd = sd(signal_vs_no_signal),
  q025 = quantile(signal_vs_no_signal, 0.025),
  q50 = quantile(signal_vs_no_signal, 0.5),
  q975 = quantile(signal_vs_no_signal, 0.975),
  prior_sd = if (lambda_structure == 1L) {
    log_ratio_sd_prior
  } else {
    log_ratio_sd_prior / sqrt(2)
  }
)
contrast_summary$posterior_to_prior_sd <-
  contrast_summary$sd / contrast_summary$prior_sd

diagnostics <- as.data.frame(summarise_draws(
  subset_draws(fit$post_warmup_draws, variable = required)
))
private_names <- intersect(
  c(
    "beta_intercept", "beta_ink_effect", "beta_bracelet_effect",
    "wtp_value_utility", "hyper_wtp_mu", "dist_beta_v[1]", "base_mu_rep"
  ),
  names(draws)
)
correlations <- do.call(rbind, lapply(seq_along(treatment), function(index) {
  data.frame(
    treatment = treatment[index],
    parameter = private_names,
    correlation = vapply(
      private_names,
      function(parameter) cor(log_deviation[, index], draws[[parameter]]),
      numeric(1)
    )
  )
}))

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(
  lambda_summary,
  file.path(output_path, "signal-lambda-summary.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics,
  file.path(output_path, "signal-lambda-sampler-diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  contrast_summary,
  file.path(output_path, "signal-lambda-contrast-summary.csv"),
  row.names = FALSE
)
write.csv(
  correlations,
  file.path(output_path, "signal-lambda-confounding-correlations.csv"),
  row.names = FALSE
)
print(lambda_summary, row.names = FALSE)
print(diagnostics, row.names = FALSE)
print(contrast_summary, row.names = FALSE)
print(correlations, row.names = FALSE)
