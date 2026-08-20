#!/usr/bin/env Rscript

source("R/structural/main-core-data.R")
source("scripts/appendix/main-core-multiplier-contrasts.R")
source("R/policy/bootstrap.R")

phi_approx <- function(x) plogis(0.07056 * x^3 + 1.5976 * x)
information_factor <- function(cutoff, total_sd, q1, q0) {
  probability <- phi_approx(-cutoff / total_sd)
  signal_probability <- probability * q1 + (1 - probability) * q0
  probability * (1 - probability) * sum((q1 - q0)^2 / signal_probability)
}

# The response-channel information factor nests the old perfect-observation
# schedule and vanishes for an uninformative channel.
stopifnot(
  abs(information_factor(0.4, 1.3, c(1, 0), c(0, 1)) - 1) < 1e-12,
  abs(information_factor(0.4, 1.3, c(0.3, 0.7), c(0.3, 0.7))) < 1e-12
)

# The design-pooled basis separates the Any-Signal/No-Signal contrast from
# the two within-group contrasts and still spans every centered four-arm slope.
group_basis <- c(-0.5, 0.5, -0.5, 0.5)
within_basis <- rbind(
  c(-1 / sqrt(2), 0), c(0, -1 / sqrt(2)),
  c(1 / sqrt(2), 0), c(0, 1 / sqrt(2))
)
design_basis <- cbind(group_basis, within_basis)
stopifnot(
  max(abs(colSums(design_basis))) < 1e-12,
  max(abs(crossprod(group_basis, within_basis))) < 1e-12,
  qr(design_basis)$rank == 3L
)
unrestricted_slope <- c(-0.4, 0.2, 0.1, 0.1)
coefficient <- qr.solve(design_basis, unrestricted_slope)
stopifnot(max(abs(design_basis %*% coefficient - unrestricted_slope)) < 1e-12)
pooled_slope <- 0.6 * group_basis
stopifnot(
  pooled_slope[1] == pooled_slope[3],
  pooled_slope[2] == pooled_slope[4],
  mean(pooled_slope[c(2, 4)]) - mean(pooled_slope[c(1, 3)]) == 0.6
)

# Paired contrasts must be computed within draws rather than by subtracting
# marginal interval endpoints.
synthetic_multiplier <- cbind(
  Control = c(1.5, 1.7), Ink = c(0.9, 1.1),
  Calendar = c(1.2, 1.4), Bracelet = c(0.8, 1.0)
)
synthetic_contrast <- main_core_multiplier_contrasts(synthetic_multiplier)
stopifnot(all.equal(
  unname(synthetic_contrast[, "No Signal - Any Signal"]), c(0.5, 0.5)
))
reveal_probability <- 0.63
null_signal_taker <- c(reveal_probability, 0, 1 - reveal_probability)
null_signal_nontaker <- c(0, reveal_probability, 1 - reveal_probability)
symmetric_taker <- c(0.8, 0.2)
symmetric_nontaker <- c(0.2, 0.8)
asymmetric_taker <- c(0.75, 0.15, 0.10)
asymmetric_nontaker <- c(0.10, 0.55, 0.35)
stopifnot(
  abs(information_factor(
    0.4, 1.3, null_signal_taker, null_signal_nontaker
  ) - reveal_probability) < 1e-12,
  information_factor(0.4, 1.3, symmetric_taker, symmetric_nontaker) > 0,
  information_factor(0.4, 1.3, asymmetric_taker, asymmetric_nontaker) > 0,
  abs(sum(symmetric_taker) - 1) < 1e-12,
  abs(sum(asymmetric_nontaker) - 1) < 1e-12
)

# Check the analytic cutoff and distance derivatives used by compact GQ.
softmax <- function(value) {
  shifted <- value - max(value)
  exp(shifted) / sum(exp(shifted))
}

# Structure 0 is the exact pre-ladder response likelihood: the same full
# recognition and multinomial logits, with no altered scaling or priors.
set.seed(20260805)
contrast <- matrix(rnorm(12), 4, 3)
recognition_intercept <- rnorm(2)
recognition_dist <- rnorm(2)
recognition_arm <- matrix(rnorm(6), 2, 3)
recognition_arm_dist <- matrix(rnorm(6), 2, 3)
report_intercept <- matrix(rnorm(4), 2, 2)
report_dist <- matrix(rnorm(4), 2, 2)
report_arm <- matrix(rnorm(12), 2, 6)
report_arm_dist <- matrix(rnorm(12), 2, 6)
distance_test <- -0.37
treatment_test <- 3L
truth_test <- 2L
old_recognition_eta <- recognition_intercept[truth_test] +
  sum(recognition_arm[truth_test, ] * contrast[treatment_test, ]) +
  (recognition_dist[truth_test] +
   sum(recognition_arm_dist[truth_test, ] * contrast[treatment_test, ])) *
  distance_test
old_report_eta <- vapply(1:2, function(category) {
  columns <- (category - 1L) * 3L + 1:3
  report_intercept[truth_test, category] +
    sum(report_arm[truth_test, columns] * contrast[treatment_test, ]) +
    (report_dist[truth_test, category] +
     sum(report_arm_dist[truth_test, columns] * contrast[treatment_test, ])) *
    distance_test
}, numeric(1))
new_recognition_eta <- recognition_intercept[truth_test]
if (0L == 0L) {
  new_recognition_eta <- new_recognition_eta +
    sum(recognition_arm[truth_test, ] * contrast[treatment_test, ]) +
    (recognition_dist[truth_test] +
     sum(recognition_arm_dist[truth_test, ] * contrast[treatment_test, ])) *
    distance_test
}
new_report_eta <- vapply(1:2, function(category) {
  columns <- (category - 1L) * 3L + 1:3
  arm_slope <- if (0L == 0L) {
    sum(report_arm_dist[truth_test, columns] * contrast[treatment_test, ])
  } else 0
  report_intercept[truth_test, category] +
    sum(report_arm[truth_test, columns] * contrast[treatment_test, ]) +
    (report_dist[truth_test, category] + arm_slope) * distance_test
}, numeric(1))
recognized <- 17L
total <- 23L
report_count <- c(8L, 6L, 3L)
old_target <- dbinom(recognized, total, plogis(old_recognition_eta), log = TRUE) +
  dmultinom(report_count, prob = softmax(c(old_report_eta, 0)), log = TRUE)
new_target <- dbinom(recognized, total, plogis(new_recognition_eta), log = TRUE) +
  dmultinom(report_count, prob = softmax(c(new_report_eta, 0)), log = TRUE)
stopifnot(
  identical(old_recognition_eta, new_recognition_eta),
  identical(old_report_eta, new_report_eta),
  identical(old_target, new_target)
)
channel <- function(distance) {
  make_truth <- function(intercept, slope, recognition_intercept,
                         recognition_slope) {
    conditional <- softmax(c(intercept + slope * distance, 0))
    conditional_derivative <- conditional *
      (c(slope, 0) - sum(conditional * c(slope, 0)))
    recognition <- plogis(recognition_intercept + recognition_slope * distance)
    recognition_derivative <- recognition * (1 - recognition) * recognition_slope
    list(
      q = c(recognition * conditional, 1 - recognition),
      qd = c(
        recognition_derivative * conditional +
          recognition * conditional_derivative,
        -recognition_derivative
      )
    )
  }
  list(
    taker = make_truth(c(1.0, -0.3), c(0.2, -0.1), 0.4, -0.2),
    nontaker = make_truth(c(-0.4, 0.7), c(-0.1, 0.15), -0.2, 0.1)
  )
}
analytic_information <- function(cutoff, total_sd, response) {
  z <- -cutoff / total_sd
  probability <- phi_approx(z)
  probability_w <- -probability * (1 - probability) *
    (1.5976 + 3 * 0.07056 * z^2) / total_sd
  difference <- response$taker$q - response$nontaker$q
  difference_d <- response$taker$qd - response$nontaker$qd
  signal_probability <- probability * response$taker$q +
    (1 - probability) * response$nontaker$q
  signal_probability_d <- probability * response$taker$qd +
    (1 - probability) * response$nontaker$qd
  h <- sum(difference^2 / signal_probability)
  h_probability <- -sum(difference^3 / signal_probability^2)
  c(
    value = probability * (1 - probability) * h,
    cutoff = probability_w * ((1 - 2 * probability) * h +
      probability * (1 - probability) * h_probability),
    distance = probability * (1 - probability) * sum(
      2 * difference * difference_d / signal_probability -
        difference^2 * signal_probability_d / signal_probability^2
    )
  )
}
cutoff <- 0.25
distance <- -0.35
total_sd <- 1.4
epsilon <- 1e-6
analytic <- analytic_information(cutoff, total_sd, channel(distance))
cutoff_fd <- (
  information_factor(cutoff + epsilon, total_sd,
                     channel(distance)$taker$q, channel(distance)$nontaker$q) -
  information_factor(cutoff - epsilon, total_sd,
                     channel(distance)$taker$q, channel(distance)$nontaker$q)
) / (2 * epsilon)
distance_fd <- (
  information_factor(cutoff, total_sd,
                     channel(distance + epsilon)$taker$q,
                     channel(distance + epsilon)$nontaker$q) -
  information_factor(cutoff, total_sd,
                     channel(distance - epsilon)$taker$q,
                     channel(distance - epsilon)$nontaker$q)
) / (2 * epsilon)
stopifnot(
  abs(analytic["cutoff"] - cutoff_fd) < 1e-7,
  abs(analytic["distance"] - distance_fd) < 1e-7
)

# A correct-or-null channel must exactly reproduce the old p-observed fixed
# point and policy demand when old mu_rep=lambda*p.
nesting_channel <- list(
  taker = matrix(rep(null_signal_taker, 3L), nrow = 3L, byrow = TRUE),
  nontaker = matrix(rep(null_signal_nontaker, 3L), nrow = 3L, byrow = TRUE)
)
nesting_benefit <- c(-0.2, 0.1, 0.4)
nesting_lambda <- 0.8
old_cutoff <- solve_policy_fixedpoint(
  nesting_benefit, rep(nesting_lambda * reveal_probability, 3L), 0.5
)
new_cutoff <- solve_policy_noisy_fixedpoint(
  nesting_benefit, nesting_lambda, 0.5, nesting_channel
)
stopifnot(max(abs(old_cutoff - new_cutoff)) < 1e-8)

# The fast policy adapter consumes the same channel and returns valid roots.
policy_parameter <- as.list(setNames(rep(0, 48), paste0(
  "observation_parameter_", seq_len(48)
)))
policy_parameter$model_family <- "asymmetric_unconditional"
policy_parameter$base_mu_rep <- 0.7
policy_parameter$u_sd <- 0.5
policy_parameter$total_error_sd <- sqrt(1.25)
policy_parameter$observation_parameter_17 <- -1
policy_parameter$observation_parameter_18 <- 1
policy_channel <- policy_noisy_channel(
  c(-0.5, 0, 0.5), policy_parameter, "control"
)
policy_channel_scalar <- policy_noisy_channel(0, policy_parameter, "bracelet")
policy_cutoff <- solve_policy_noisy_fixedpoint(
  benefit = c(0.1, 0.2, 0.3), lambda = policy_parameter$base_mu_rep,
  u_sd = policy_parameter$u_sd, channel = policy_channel
)
stopifnot(
  max(abs(rowSums(policy_channel$taker) - 1)) < 1e-12,
  max(abs(rowSums(policy_channel$nontaker) - 1)) < 1e-12,
  nrow(policy_channel_scalar$taker) == 1L,
  abs(sum(policy_channel_scalar$taker) - 1) < 1e-12,
  all(is.finite(policy_cutoff)),
  max(abs(policy_noisy_residual(
    policy_cutoff, c(0.1, 0.2, 0.3), policy_parameter$base_mu_rep,
    policy_parameter$u_sd, policy_channel
  ))) < 1e-8
)

# Reconstruct and audit the actual no-outlier peer sample when it is available.
workspace <- "data/stan_analysis_data/dist_fit105.RData"
if (file.exists(workspace)) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(rlang)
  })
  baseline <- prepare_main_core_data(
    workspace, "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS",
    observation_model = 0L
  )
  asymmetric <- prepare_main_core_data(
    workspace, "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS",
    observation_model = 1L
  )
  stopifnot(
    baseline$core_num_peer_response_rows == 0L,
    asymmetric$core_num_peer_response_rows > 0L,
    sum(asymmetric$core_peer_report_count) ==
      sum(asymmetric$core_peer_recognized),
    all(rowSums(asymmetric$core_peer_report_count) ==
      asymmetric$core_peer_recognized)
  )
}

cat("Main-core asymmetric observability tests: OK\n")
