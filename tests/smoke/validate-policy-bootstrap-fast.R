#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/policy/bootstrap.R")
suppressPackageStartupMessages(library(nleqslv))

parameter_csv <- policy_option_value(args, "--parameter-csv")
distance_data <- policy_option_value(args, "--distance-data", "optim/data/full-many-pots-experiment.rds")
num_replicates <- as.integer(policy_option_value(args, "--num-replicates", "5"))
tolerance <- as.numeric(policy_option_value(args, "--tolerance", "1e-8"))
output_csv <- policy_option_value(args, "--output-csv")
if (is.null(parameter_csv)) stop("--parameter-csv is required.", call. = FALSE)

parameters <- read.csv(parameter_csv, stringsAsFactors = FALSE)
parameters <- parameters[seq_len(min(num_replicates, nrow(parameters))), ]
distance_object <- readRDS(distance_data)
distances <- sort(unique(
  distance_object$long_distance_mat$dist[
    distance_object$long_distance_mat$dist <= 3500
  ]
))
parameters$sd_of_dist <- distance_object$sd_of_dist

checks <- list()
check_index <- 1L
for (row in seq_len(nrow(parameters))) {
  parameter <- parameters[row, ]
  distance_sd <- distances / parameter$sd_of_dist
  benefit <- parameter$beta_control - parameter$dist_beta * distance_sd
  for (visibility in c("control", "bracelet")) {
    mu <- policy_mu_rep(distance_sd, parameter, visibility)
    started_fast <- proc.time()[[3L]]
    fast <- solve_policy_fixedpoint(benefit, mu, parameter$u_sd)
    fast_seconds <- proc.time()[[3L]] - started_fast
    started_legacy <- proc.time()[[3L]]
    legacy_fit <- lapply(seq_along(benefit), function(index) {
      nleqslv(
        x = -benefit[index],
        fn = function(cutoff) policy_fixedpoint_residual(
          cutoff, benefit[index], mu[index], parameter$u_sd
        )
      )
    })
    legacy_seconds <- proc.time()[[3L]] - started_legacy
    term_codes <- vapply(legacy_fit, `[[`, numeric(1), "termcd")
    if (any(term_codes > 2)) {
      stop("Legacy nleqslv failed for retained validation points.", call. = FALSE)
    }
    legacy <- vapply(legacy_fit, function(fit) fit$x, numeric(1))
    fast_demand <- 1 - pnorm(fast / parameter$total_error_sd)
    legacy_demand <- 1 - pnorm(legacy / parameter$total_error_sd)
    checks[[check_index]] <- data.frame(
      draw = parameter$draw,
      replicate = parameter$replicate,
      visibility = visibility,
      distances = length(distances),
      max_cutoff_difference = max(abs(fast - legacy)),
      max_demand_difference = max(abs(fast_demand - legacy_demand)),
      fast_seconds = fast_seconds,
      legacy_seconds = legacy_seconds,
      speedup = legacy_seconds / max(fast_seconds, .Machine$double.eps),
      stringsAsFactors = FALSE
    )
    check_index <- check_index + 1L
  }
}
checks <- do.call(rbind, checks)
if (max(checks$max_demand_difference) > tolerance) {
  stop("Fast demand differs from legacy nleqslv by more than ", tolerance, call. = FALSE)
}
if (!is.null(output_csv)) {
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(checks, output_csv, row.names = FALSE)
}
print(checks)
message(
  "Validated ", nrow(parameters), " draws at ", length(distances),
  " distances; median solver speedup = ",
  round(median(checks$speedup), 1), "x."
)
