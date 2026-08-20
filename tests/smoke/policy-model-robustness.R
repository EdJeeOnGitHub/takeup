#!/usr/bin/env Rscript

source("R/policy/bootstrap.R")

mixture <- policy_student_t_mixture(df = 5, components = 12)
grid <- seq(-4, 4, by = 0.05)
mixture_cdf <- policy_student_t_moments(grid, u_sd = 0, mixture)$probability_below
exact_cdf <- pt(grid / sqrt(3 / 5), df = 5)
stopifnot(max(abs(mixture_cdf - exact_cdf)) < 3e-4)

draw <- data.frame(
  beta_intercept = 0, beta_ink_effect = 0, beta_bracelet_effect = 0,
  wtp_value_utility = 0, hyper_wtp_mu = 0, base_mu_rep = 0.8,
  raw_u_sd.1 = 0.4, dist_beta_v.1 = 0.2
)
for (index in 1:4) {
  draw[[paste0("hyper_beta_1ord.", index)]] <- 0
  draw[[paste0("hyper_dist_beta_1ord.", index)]] <- 0
}
draw$core_lambda_group_log_ratio_raw.1 <- 1
grouped <- canonical_policy_draw(
  draw, 1, 1, "grouped", "Grouped", lambda_structure = "grouped",
  lambda_log_ratio_sd_prior = 0.25
)
stopifnot(
  abs(log(grouped$lambda_bracelet / grouped$lambda_control) - 0.25) < 1e-12,
  abs(sqrt(grouped$lambda_bracelet * grouped$lambda_control) - 0.8) < 1e-12
)

benefit <- seq(-1, 1, length.out = 31)
mu <- rep(0.4, length(benefit))
cutoff <- solve_policy_student_t_fixedpoint(benefit, mu, 0.4, mixture)
moments <- policy_student_t_moments(cutoff, 0.4, mixture)
stopifnot(
  all(is.finite(cutoff)),
  max(abs(cutoff + benefit + mu * moments$delta)) < 2e-5,
  all(moments$probability_below > 0 & moments$probability_below < 1)
)

# A saddle-node counterfactual has no exact fixed point. It must be marked
# undefined rather than replaced by the numerical minimum-residual point.
undefined <- solve_policy_fixedpoint(
  benefit = -0.7245538, mu_rep = 1.285234, u_sd = 0.379876
)
stopifnot(is.na(undefined), attr(undefined, "undefined_count") == 1L)

# Stan's fit-106 cluster index and the policy village row order contain the
# same 144 external clusters but are not positionally identical.
if (file.exists("data/stan_analysis_data/dist_fit106.RData")) {
  fit_environment <- new.env(parent = emptyenv())
  suppressWarnings(load(
    "data/stan_analysis_data/dist_fit106.RData", envir = fit_environment
  ))
  cluster_map <- unique(
    fit_environment$stan_data$analysis_data[, c("cluster_id", "cluster.id")]
  )
  cluster_map <- cluster_map[order(cluster_map$cluster_id), ]
  distance_object <- readRDS("optim/data/full-many-pots-experiment.rds")
  village_id <- distance_object$village_df$cluster.id
  stopifnot(
    identical(cluster_map$cluster_id, seq_len(144L)),
    setequal(cluster_map$cluster.id, village_id),
    !identical(cluster_map$cluster.id, village_id)
  )
}

cat("Policy model-robustness transforms and Student-t predictor passed.\n")
