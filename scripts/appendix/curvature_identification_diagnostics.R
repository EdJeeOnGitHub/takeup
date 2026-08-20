#!/usr/bin/env Rscript

args <- docopt::docopt(
  "
Curvature-based identification diagnostics for the structural model.

Usage:
  scripts/appendix/curvature_identification_diagnostics.R <fit-version> --model=<model> [options]

Options:
  --fit-path=<path>       Directory containing fitted CmdStan CSV files [default: data/stan_analysis_data]
  --analysis-file=<path>  Saved scripts/structural/run-model.R workspace [default: data/stan_analysis_data/dist_fit{fit-version}.RData]
  --output-path=<path>    Output directory [default: temp-data/bayesian-sensitivity]
  --chains=<chains>       Comma-separated chain numbers [default: 1,2,3,4]
  --num-clusters=<n>      Number of clusters [default: 144]
  --digits=<n>            Digits in formatted tables [default: 2]
  ",
  args = commandArgs(trailingOnly = TRUE)
)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(posterior)
  library(readr)
  library(tidyr)
})

fit_version <- args$fit_version
model_name <- args$model
chains <- as.integer(strsplit(args$chains, ",", fixed = TRUE)[[1]])
num_clusters <- as.integer(args$num_clusters)
digits <- as.integer(args$digits)
analysis_file <- stringr::str_replace_all(
  args$analysis_file,
  stringr::fixed("{fit-version}"),
  fit_version
)

fit_files <- file.path(
  args$fit_path,
  sprintf("dist_fit%s_%s-%d.csv", fit_version, model_name, chains)
)
missing_files <- c(
  fit_files[!file.exists(fit_files)],
  analysis_file[!file.exists(analysis_file)]
)
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

fit_environment <- new.env(parent = globalenv())
load(analysis_file, envir = fit_environment)
stan_data <- get("stan_data", envir = fit_environment)

if (stan_data$num_clusters != num_clusters) {
  stop(
    "--num-clusters=", num_clusters,
    " does not match saved Stan data: ", stan_data$num_clusters
  )
}

cluster_n <- as.numeric(tabulate(stan_data$obs_cluster_id, nbins = num_clusters))
cluster_y <- as.numeric(rowsum(
  as.integer(stan_data$takeup),
  stan_data$obs_cluster_id,
  reorder = FALSE
))
distance <- as.numeric(stan_data$cluster_standard_dist)
treatment <- factor(
  stan_data$cluster_assigned_treatment,
  levels = c("control", "ink", "calendar", "bracelet")
)
arm_matrix <- model.matrix(~ 0 + treatment)

weighted_residual <- function(y, x, weights) {
  x <- as.matrix(x)
  sqrt_w <- sqrt(weights)
  fit <- qr(sqrt_w * x, LAPACK = TRUE)
  y - as.numeric(x %*% qr.coef(fit, sqrt_w * y))
}

weighted_rms <- function(x, weights) {
  sqrt(sum(weights * x^2) / sum(weights))
}

# Level is the original take-up likelihood. Slope is within-arm centered
# distance. Curvature is squared distance residualized against arm-specific
# levels and slopes. Slope and curvature are scaled to weighted RMS one, so an
# epsilon perturbation has comparable observation-weight magnitude.
level_weight <- rep(1, num_clusters)
slope_weight <- weighted_residual(distance, arm_matrix, cluster_n)
slope_weight <- slope_weight / weighted_rms(slope_weight, cluster_n)
arm_slope_matrix <- arm_matrix * distance
curvature_weight <- weighted_residual(
  distance^2,
  cbind(arm_matrix, arm_slope_matrix),
  cluster_n
)
curvature_weight <-
  curvature_weight / weighted_rms(curvature_weight, cluster_n)

contrast_weights <- cbind(
  level = level_weight,
  slope = slope_weight,
  curvature = curvature_weight
)

parameter_map <- c(
  lambda = "base_mu_rep",
  delta = "dist_beta_v[1]",
  sigma_u = "u_sd[1]",
  beta_control = "beta[1]",
  beta_ink = "beta[2]",
  beta_calendar = "beta[3]",
  beta_bracelet = "beta[4]",
  visibility_level_control = "centered_cluster_beta_beliefs[1,1]",
  visibility_level_ink = "centered_cluster_beta_beliefs[1,2]",
  visibility_level_calendar = "centered_cluster_beta_beliefs[1,3]",
  visibility_level_bracelet = "centered_cluster_beta_beliefs[1,4]",
  visibility_distance_control = "centered_cluster_dist_beta_beliefs[1,1]",
  visibility_distance_ink = "centered_cluster_dist_beta_beliefs[1,2]",
  visibility_distance_calendar = "centered_cluster_dist_beta_beliefs[1,3]",
  visibility_distance_bracelet = "centered_cluster_dist_beta_beliefs[1,4]",
  psi = "hyper_wtp_mu"
)

derived_spec <- tidyr::crossing(
  measure = c("social_multiplier", "social_image_return"),
  distance_m = c(500L, 1500L, 2500L),
  treatment_id = seq_len(4)
) %>%
  mutate(
    distance_index = distance_m %/% 100L + 1L,
    treatment = c("Control", "Ink", "Calendar", "Bracelet")[treatment_id],
    stan_base = if_else(
      measure == "social_multiplier",
      "cluster_social_multiplier",
      "cluster_rep_return"
    ),
    quantity = paste(measure, tolower(treatment), distance_m, sep = "_"),
    quantity_label = if_else(
      measure == "social_multiplier",
      sprintf("M_{%s}(%0.1f\\,\\mathrm{km})", treatment, distance_m / 1000),
      sprintf("SI_{%s}(%0.1f\\,\\mathrm{km})", treatment, distance_m / 1000)
    )
  )

derived_variables <- unlist(lapply(seq_len(nrow(derived_spec)), function(i) {
  row <- derived_spec[i, ]
  sprintf(
    "%s[%d,%d,%d]",
    row$stan_base,
    row$distance_index,
    seq_len(num_clusters),
    row$treatment_id
  )
}))
probability_variables <- sprintf(
  "structural_cluster_takeup_prob[%d]",
  seq_len(num_clusters)
)

parameter_draws_by_chain <- vector("list", length(fit_files))
derived_draws_by_chain <- vector("list", length(fit_files))
probability_draws_by_chain <- vector("list", length(fit_files))
contrast_draws_by_chain <- vector("list", length(fit_files))

for (chain_index in seq_along(fit_files)) {
  fit_draws <- read_cmdstan_csv(
    fit_files[chain_index],
    variables = c(
      unname(parameter_map),
      probability_variables,
      derived_variables
    )
  )$post_warmup_draws %>%
    as_draws_matrix()

  parameter_draws_by_chain[[chain_index]] <-
    fit_draws[, unname(parameter_map), drop = FALSE]
  probability_draws <- fit_draws[, probability_variables, drop = FALSE]
  probability_draws_by_chain[[chain_index]] <- probability_draws

  probability_safe <- pmin(pmax(probability_draws, 1e-12), 1 - 1e-12)
  cluster_log_lik <-
    sweep(log(probability_safe), 2, cluster_y, "*") +
    sweep(log1p(-probability_safe), 2, cluster_n - cluster_y, "*")
  contrast_draws_by_chain[[chain_index]] <-
    cluster_log_lik %*% contrast_weights

  chain_derived <- matrix(
    NA_real_,
    nrow = nrow(fit_draws),
    ncol = nrow(derived_spec),
    dimnames = list(NULL, derived_spec$quantity)
  )
  for (row_index in seq_len(nrow(derived_spec))) {
    row <- derived_spec[row_index, ]
    cluster_variables <- sprintf(
      "%s[%d,%d,%d]",
      row$stan_base,
      row$distance_index,
      seq_len(num_clusters),
      row$treatment_id
    )
    cluster_average <-
      rowMeans(fit_draws[, cluster_variables, drop = FALSE])
    chain_derived[, row_index] <- if (row$measure == "social_multiplier") {
      -cluster_average / fit_draws[, "dist_beta_v[1]"]
    } else {
      cluster_average
    }
  }
  derived_draws_by_chain[[chain_index]] <- chain_derived
  rm(fit_draws, probability_draws, probability_safe, cluster_log_lik)
  gc(verbose = FALSE)
}

parameter_draws <- do.call(rbind, parameter_draws_by_chain)
colnames(parameter_draws) <- names(parameter_map)
derived_draws <- do.call(rbind, derived_draws_by_chain)
probability_draws <- do.call(rbind, probability_draws_by_chain)
contrast_draws <- do.call(rbind, contrast_draws_by_chain)
colnames(contrast_draws) <- colnames(contrast_weights)

selected_draws <- cbind(
  parameter_draws[, c("sigma_u", "lambda", "delta"), drop = FALSE],
  derived_draws
)
selected_labels <- c(
  sigma_u = "\\sigma_u",
  lambda = "\\lambda",
  delta = "\\delta",
  setNames(derived_spec$quantity_label, derived_spec$quantity)
)

calculate_sensitivity <- function(draws, contrasts) {
  posterior_mean <- colMeans(draws)
  posterior_sd <- apply(draws, 2, sd)
  posterior_variance <- posterior_sd^2

  mean_sensitivity <- cov(draws, contrasts)
  mean_sensitivity <- sweep(mean_sensitivity, 1, posterior_sd, "/")
  colnames(mean_sensitivity) <- paste0("mean_", colnames(mean_sensitivity))

  centered_squared <- sweep(draws, 2, posterior_mean, "-")^2
  log_sd_sensitivity <- cov(centered_squared, contrasts)
  log_sd_sensitivity <- sweep(
    log_sd_sensitivity,
    1,
    2 * posterior_variance,
    "/"
  )
  colnames(log_sd_sensitivity) <-
    paste0("log_sd_", colnames(log_sd_sensitivity))

  bind_cols(
    tibble(
      quantity = colnames(draws),
      quantity_label = unname(selected_labels[colnames(draws)]),
      posterior_mean = unname(posterior_mean),
      posterior_sd = unname(posterior_sd)
    ),
    as.data.frame(mean_sensitivity),
    as.data.frame(log_sd_sensitivity)
  )
}

sensitivity <- calculate_sensitivity(selected_draws, contrast_draws)

# Posterior-local Jacobian: regress each cluster's predicted take-up on the 16
# standardized structural parameters. This estimates the model-moment
# derivative in the neighborhood explored by the posterior.
standardized_parameters <- scale(parameter_draws)
svd_inverse <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  decomposition <- svd(x)
  keep <- decomposition$d > max(decomposition$d) * tolerance
  decomposition$v[, keep, drop = FALSE] %*%
    diag(1 / decomposition$d[keep], nrow = sum(keep)) %*%
    t(decomposition$u[, keep, drop = FALSE])
}
jacobian <- t(
  svd_inverse(standardized_parameters) %*%
    scale(probability_draws, center = TRUE, scale = FALSE)
)
colnames(jacobian) <- colnames(parameter_draws)

mean_probability <- colMeans(probability_draws)
information_weight <- cluster_n /
  pmax(mean_probability * (1 - mean_probability), 1e-8)

sigma_derivative <- jacobian[, "sigma_u"]
other_derivatives <- jacobian[, setdiff(colnames(jacobian), "sigma_u"), drop = FALSE]
sqrt_information_weight <- sqrt(information_weight)
projection_coefficient <- svd_inverse(
  sqrt_information_weight * other_derivatives
) %*% (sqrt_information_weight * sigma_derivative)
conditional_sigma_derivative <- as.numeric(
  sigma_derivative - other_derivatives %*% projection_coefficient
)

weighted_project <- function(y, x, weights) {
  x <- as.matrix(x)
  coefficient <- svd_inverse(sqrt(weights) * x) %*% (sqrt(weights) * y)
  projected <- x %*% coefficient
  if (is.null(dim(y))) as.numeric(projected) else projected
}
weighted_norm_sq <- function(x, weights) sum(weights * x^2)

level_basis <- arm_matrix
slope_basis <- arm_slope_matrix -
  weighted_project(arm_slope_matrix, level_basis, information_weight)
level_slope_basis <- cbind(level_basis, slope_basis)
curvature_basis_raw <- arm_matrix * distance^2
curvature_basis <- curvature_basis_raw -
  weighted_project(
    curvature_basis_raw,
    level_slope_basis,
    information_weight
  )

level_component <- weighted_project(
  conditional_sigma_derivative,
  level_basis,
  information_weight
)
residual_after_level <- conditional_sigma_derivative - level_component
slope_component <- weighted_project(
  residual_after_level,
  slope_basis,
  information_weight
)
residual_after_slope <- residual_after_level - slope_component
curvature_component <- weighted_project(
  residual_after_slope,
  curvature_basis,
  information_weight
)
remainder_component <- residual_after_slope - curvature_component

total_information <- weighted_norm_sq(
  conditional_sigma_derivative,
  information_weight
)
unconditional_information <- weighted_norm_sq(
  sigma_derivative,
  information_weight
)
decomposition <- tibble(
  component = c("Level", "Slope", "Curvature", "Higher-order/remainder"),
  conditional_information = c(
    weighted_norm_sq(level_component, information_weight),
    weighted_norm_sq(slope_component, information_weight),
    weighted_norm_sq(curvature_component, information_weight),
    weighted_norm_sq(remainder_component, information_weight)
  )
) %>%
  mutate(
    share = conditional_information / total_information,
    total_conditional_information = total_information,
    unconditional_sigma_information = unconditional_information,
    conditional_to_unconditional_ratio =
      total_information / unconditional_information
  )

weight_diagnostics <- tibble(
  contrast = colnames(contrast_weights),
  weighted_mean = colSums(cluster_n * contrast_weights) / sum(cluster_n),
  weighted_rms = sqrt(
    colSums(cluster_n * contrast_weights^2) / sum(cluster_n)
  )
)

dir.create(args$output_path, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(
  args$output_path,
  sprintf("curvature-identification-fit%s-%s", fit_version, model_name)
)
write_csv(sensitivity, paste0(stem, "-sensitivity.csv"))
write_csv(decomposition, paste0(stem, "-jacobian-decomposition.csv"))
write_csv(weight_diagnostics, paste0(stem, "-weight-diagnostics.csv"))

format_number <- function(x) formatC(x, format = "f", digits = digits)
write_markdown <- function(data, columns, labels, path) {
  displayed <- data[, columns, drop = FALSE]
  displayed[] <- lapply(displayed, format_number)
  displayed <- cbind(Quantity = paste0("$", data$quantity_label, "$"), displayed)
  lines <- c(
    paste0("| ", paste(c("Quantity", labels), collapse = " | "), " |"),
    paste0("|:--|", paste(rep("--:", length(labels)), collapse = "|"), "|"),
    apply(displayed, 1, function(row) {
      paste0("| ", paste(row, collapse = " | "), " |")
    })
  )
  writeLines(lines, path)
}
write_markdown(
  sensitivity,
  c("mean_level", "mean_slope", "mean_curvature"),
  c("Level", "Slope", "Curvature"),
  paste0(stem, "-mean-sensitivity.md")
)
write_markdown(
  sensitivity,
  c("log_sd_level", "log_sd_slope", "log_sd_curvature"),
  c("Level", "Slope", "Curvature"),
  paste0(stem, "-uncertainty-sensitivity.md")
)

decomposition_markdown <- c(
  "| Component | Conditional information | Share of conditional information |",
  "|:--|--:|--:|",
  apply(decomposition, 1, function(row) {
    paste0(
      "| ", row[["component"]],
      " | ", format(
        as.numeric(row[["conditional_information"]]),
        scientific = TRUE,
        digits = 3
      ),
      " | ", formatC(
        as.numeric(row[["share"]]),
        format = "f",
        digits = 3
      ),
      " |"
    )
  }),
  "",
  sprintf(
    "Conditional-to-unconditional sigma information ratio: **%.3e**.",
    decomposition$conditional_to_unconditional_ratio[1]
  )
)
writeLines(
  decomposition_markdown,
  paste0(stem, "-jacobian-decomposition.md")
)

message("Wrote curvature sensitivity and Jacobian diagnostics to ", args$output_path)
