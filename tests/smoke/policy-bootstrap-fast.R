#!/usr/bin/env Rscript

source("R/policy/bootstrap.R")

parameter <- data.frame(
  draw = 1L,
  replicate = 17L,
  beta_control = 0.35,
  beta_ink = 0.1,
  beta_calendar = 0.2,
  beta_bracelet = 0.15,
  dist_beta = 0.22,
  mu_control = -0.25,
  mu_ink = 0.15,
  mu_calendar = 0.05,
  mu_bracelet = 0.35,
  mu_dist_control = -0.10,
  mu_dist_ink = -0.05,
  mu_dist_calendar = -0.03,
  mu_dist_bracelet = -0.08,
  base_mu_rep = 0.6,
  u_sd = 1.1,
  total_error_sd = sqrt(1 + 1.1^2),
  sd_of_dist = 664
)
distances <- seq(25, 3500, length.out = 200)
distance_sd <- distances / parameter$sd_of_dist
benefit <- parameter$beta_control - parameter$dist_beta * distance_sd
mu <- policy_mu_rep(distance_sd, parameter, "control")

fast <- solve_policy_fixedpoint(benefit, mu, parameter$u_sd)
reference <- vapply(seq_along(benefit), function(index) {
  uniroot(
    function(cutoff) policy_fixedpoint_residual(
      cutoff, benefit[index], mu[index], parameter$u_sd
    ),
    interval = c(-20, 20), tol = 1e-12
  )$root
}, numeric(1))
if (max(abs(fast - reference)) > 1e-9) {
  stop("Vectorized fixed-point solver does not match scalar reference.", call. = FALSE)
}

predictions <- predict_policy_draw(parameter, distances)
if (nrow(predictions) != length(distances) * nrow(policy_scenarios) ||
    any(!is.finite(predictions$demand)) ||
    any(predictions$demand < 0 | predictions$demand > 1)) {
  stop("Invalid policy predictions.", call. = FALSE)
}

no_rep <- predictions[predictions$scenario == "suppress-reputation", ]
closed_form <- pnorm(benefit / parameter$total_error_sd)
if (max(abs(no_rep$demand - closed_form)) > 1e-12) {
  stop("No-reputation scenario is not using the exact closed form.", call. = FALSE)
}

for (scenario in c("static-control", "static-bracelet")) {
  value <- predictions[predictions$scenario == scenario, ]
  implied_cutoff <- qnorm(1 - value$demand) * parameter$total_error_sd
  if (max(abs(diff(implied_cutoff + benefit))) > 1e-10) {
    stop("Static scenario does not have a distance-invariant signal: ", scenario,
         call. = FALSE)
  }
}

message("Fast policy fixed-point and closed-form tests passed.")
