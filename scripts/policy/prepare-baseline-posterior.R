#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/policy/bootstrap.R")
value <- function(name, default = NULL) policy_option_value(args, name, default)
fit_path <- value("--fit-path", "build/structural-fit/assigned")
output_path <- value("--output-path", "build/policy/posterior")
model <- value("--model", "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP")
draws_per_chain <- as.integer(value("--draws-per-chain", "50"))
distance_data <- value("--distance-data", "optim/data/full-many-pots-experiment.rds")
files <- file.path(fit_path, sprintf("dist_fit105_%s-%d.csv", model, 1:4))
if (!all(file.exists(files))) stop("Missing slim baseline chain(s).")
if (!file.exists(distance_data)) stop("Missing policy distance data.")
distance_sd <- readRDS(distance_data)$sd_of_dist

selected <- list()
for (chain in seq_along(files)) {
  x <- read.csv(files[[chain]], comment.char = "#", check.names = FALSE)
  if (nrow(x) < draws_per_chain) stop("Chain has fewer retained draws than requested.")
  comparison <- unique(round(seq(1, nrow(x), length.out = min(50L, nrow(x)))))
  index <- if (draws_per_chain > 50L) {
    remaining <- setdiff(seq_len(nrow(x)), comparison)
    sort(c(comparison, remaining[round(seq(1, length(remaining), length.out = draws_per_chain - length(comparison)))]))
  } else unique(round(seq(1, nrow(x), length.out = draws_per_chain)))
  if (length(index) != draws_per_chain) stop("Could not select balanced draws.")
  rows <- do.call(rbind, lapply(seq_along(index), function(j) {
    canonical_policy_draw(
      x[index[[j]], , drop = FALSE], draw_id = (chain - 1L) * draws_per_chain + j,
      chain = chain, model_id = "benchmark", model_label = "Benchmark",
      source_csv = normalizePath(files[[chain]])
    )
  }))
  rows$iteration <- index
  rows$comparison_200 <- index %in% comparison
  selected[[chain]] <- rows
}
parameters <- do.call(rbind, selected)
parameters$replicate <- parameters$draw
parameters$sd_of_dist <- distance_sd
numeric_parameter <- names(parameters)[vapply(parameters, is.numeric, logical(1))]
numeric_parameter <- setdiff(numeric_parameter, c("draw", "replicate", "chain", "iteration"))
median_row <- parameters[1, , drop = FALSE]
for (column in numeric_parameter) median_row[[column]] <- median(parameters[[column]])
median_row$draw <- 1L
median_row$replicate <- 1L
median_row$chain <- NA_integer_
median_row$iteration <- NA_integer_
median_row$mode_csv <- "componentwise posterior median"
median_row$source_csv <- "componentwise posterior median"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
write.csv(parameters, file.path(output_path, "policy-posterior-parameters.csv"), row.names = FALSE)
write.csv(median_row, file.path(output_path, "policy-posterior-median-parameters.csv"), row.names = FALSE)
write.csv(data.frame(
  chains = 4L, draws_per_chain = draws_per_chain, draws = nrow(parameters),
  selection = "evenly spaced retained draws within each chain"
), file.path(output_path, "policy-posterior-manifest.csv"), row.names = FALSE)
message("Wrote ", nrow(parameters), " balanced posterior policy draws.")

write.csv(parameters, file.path(output_path, "policy-model-parameters.csv"), row.names = FALSE)
saveRDS(list(parameters = parameters, cluster_shock = NULL, cluster_external_id = NULL),
        file.path(output_path, "policy-model-parameters.rds"), compress = FALSE)
write.csv(data.frame(model_id = "benchmark", model_label = "Benchmark", model_family = "gaussian",
                    draws = nrow(parameters), chains = 4L),
          file.path(output_path, "policy-model-parameter-status.csv"), row.names = FALSE)
write.csv(data.frame(path = normalizePath(files), md5 = unname(tools::md5sum(files))),
          file.path(output_path, "policy-source-fits.csv"), row.names = FALSE)
