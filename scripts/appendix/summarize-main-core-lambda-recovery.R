#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
manifest_path <- main_core_option_value(args, "--manifest")
output_path <- main_core_option_value(
  args, "--output-path", "temp-data/main-core-lambda-identification/recovery"
)
if (is.null(manifest_path) || !file.exists(manifest_path)) {
  stop("--manifest must name an existing CSV.", call. = FALSE)
}
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(ggplot2)
})

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
treatments <- c("Control", "Ink", "Calendar", "Bracelet")
rows <- list()
hmc_rows <- list()
for (index in seq_len(nrow(manifest))) {
  task <- manifest[index, ]
  run_path <- file.path(output_path, "fits", task$label)
  status_path <- file.path(run_path, "status.csv")
  if (!file.exists(status_path)) next
  status <- read.csv(status_path, stringsAsFactors = FALSE)
  if (status$status[1] != "complete" || !file.exists(status$gq_csv[1])) next
  gq <- tryCatch(
    read_cmdstan_csv(status$gq_csv[1]),
    error = function(error) stop(
      "Could not read recovery GQ for ", task$label, ": ",
      conditionMessage(error), call. = FALSE
    )
  )
  draw <- tryCatch(
    as_draws_df(gq$generated_quantities),
    error = function(error) stop(
      "Could not convert recovery GQ for ", task$label, ": ",
      conditionMessage(error), call. = FALSE
    )
  )
  estimate_lambda <- as.numeric(unlist(draw[1, paste0(
    "core_compact_signal_lambda[", 1:4, "]"
  )], use.names = FALSE))
  estimate_sm <- -as.numeric(unlist(draw[1, paste0(
    "core_compact_sm_rescaled[2,", 1:4, "]"
  )], use.names = FALSE))
  true_lambda <- as.numeric(unlist(task[paste0(
    "true_lambda_", tolower(treatments)
  )], use.names = FALSE))
  true_sm <- as.numeric(unlist(task[paste0(
    "true_sm_", tolower(treatments), "_1500m"
  )], use.names = FALSE))
  for (treatment_index in 1:4) {
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = task$scenario,
      replicate = task$replicate,
      fit_structure = task$fit_structure,
      treatment = treatments[treatment_index],
      true_lambda = true_lambda[treatment_index],
      estimate_lambda = estimate_lambda[treatment_index],
      lambda_error = estimate_lambda[treatment_index] - true_lambda[treatment_index],
      true_multiplier = true_sm[treatment_index],
      estimate_multiplier = estimate_sm[treatment_index],
      multiplier_error = estimate_sm[treatment_index] - true_sm[treatment_index],
      true_signal_log_ratio = task$true_signal_log_ratio,
      estimate_signal_log_ratio = as.numeric(
        draw$core_compact_signal_vs_no_signal_log_lambda[1]
      )
    )
  }

  if (task$replicate <= 2L) {
    hmc_path <- file.path(output_path, "hmc", task$label)
    hmc_csvs <- sort(list.files(hmc_path, pattern = "\\.csv$", full.names = TRUE))
    if (length(hmc_csvs) == 4L) {
      # A provisional summary may run while audit chains are still writing.
      # Skip that diagnostic cell until all four CSVs are readable together.
      hmc <- tryCatch(read_cmdstan_csv(hmc_csvs), error = function(error) NULL)
      if (!is.null(hmc)) {
        diagnostic_row <- tryCatch({
          summary <- summarise_draws(hmc$post_warmup_draws)
          sampler <- as_draws_df(hmc$post_warmup_sampler_diagnostics)
          data.frame(
            scenario = task$scenario,
            replicate = task$replicate,
            fit_structure = task$fit_structure,
            max_rhat = max(summary$rhat, na.rm = TRUE),
            min_ess_bulk = min(summary$ess_bulk, na.rm = TRUE),
            min_ess_tail = min(summary$ess_tail, na.rm = TRUE),
            divergences = sum(sampler$divergent__, na.rm = TRUE),
            max_treedepth = sum(sampler$treedepth__ >= 12, na.rm = TRUE)
          )
        }, error = function(error) NULL)
        if (!is.null(diagnostic_row)) {
          hmc_rows[[length(hmc_rows) + 1L]] <- diagnostic_row
        }
      }
    }
  }
}
message("Recovery files read: ", length(rows), " estimand rows")
results <- bind_rows(rows)
message("Recovery estimand rows combined: ", nrow(results))
if (nrow(results) == 0L) stop("No completed recovery fits.", call. = FALSE)
metrics <- results |>
  group_by(scenario, fit_structure) |>
  summarise(
    successful_replications = n_distinct(replicate),
    lambda_bias = mean(lambda_error),
    lambda_rmse = sqrt(mean(lambda_error^2)),
    multiplier_bias = mean(multiplier_error),
    multiplier_rmse = sqrt(mean(multiplier_error^2)),
    signal_contrast_bias = mean(
      estimate_signal_log_ratio - true_signal_log_ratio,
      na.rm = TRUE
    ),
    signal_contrast_rmse = sqrt(mean(
      (estimate_signal_log_ratio - true_signal_log_ratio)^2,
      na.rm = TRUE
    )),
    .groups = "drop"
  )
message("Recovery metrics computed: ", nrow(metrics), " cells")
hmc_diagnostics <- bind_rows(hmc_rows)
write.csv(results, file.path(output_path, "lambda-recovery-results.csv"), row.names = FALSE)
write.csv(metrics, file.path(output_path, "lambda-recovery-metrics.csv"), row.names = FALSE)
write.csv(
  hmc_diagnostics,
  file.path(output_path, "lambda-recovery-hmc-diagnostics.csv"),
  row.names = FALSE
)
message("Recovery CSV artifacts written")

plot_data <- results |>
  distinct(
    scenario, replicate, fit_structure,
    true_signal_log_ratio, estimate_signal_log_ratio
  ) |>
  filter(is.finite(true_signal_log_ratio))
message("Recovery plot data prepared")
p1 <- ggplot(
  plot_data,
  aes(true_signal_log_ratio, estimate_signal_log_ratio, color = fit_structure)
) +
  geom_abline(slope = 1, intercept = 0, color = "grey55") +
  geom_jitter(width = 0.012, height = 0, alpha = 0.45, size = 0.9) +
  facet_wrap(~scenario) +
  labs(
    x = "True signal/no-signal log-lambda contrast",
    y = "Recovered contrast", color = "Fitted structure"
  ) + theme_bw(base_size = 9) + theme(legend.position = "bottom")
p2 <- ggplot(
  results,
  aes(true_multiplier, estimate_multiplier, color = fit_structure)
) +
  geom_abline(slope = 1, intercept = 0, color = "grey55") +
  geom_point(alpha = 0.35, size = 0.75) +
  facet_grid(scenario ~ treatment, scales = "free") +
  labs(x = "True multiplier at 1500m", y = "Recovered multiplier", color = NULL) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom")
ggsave(file.path(output_path, "main-core-lambda-recovery-contrast.pdf"), p1,
       width = 7.2, height = 4.8)
ggsave(file.path(output_path, "main-core-lambda-recovery-multipliers.pdf"), p2,
       width = 7.2, height = 8.2)

fmt <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "--")
lines <- c(
  "\\begin{tabular}{llrrrr}", "\\toprule",
  "DGP & Fitted model & $N$ & Bias($\\lambda$) & RMSE($\\lambda$) & RMSE(multiplier) \\\\",
  "\\midrule"
)
for (index in seq_len(nrow(metrics))) {
  row <- metrics[index, ]
  lines <- c(lines, paste0(
    row$scenario, " & ", row$fit_structure, " & ",
    row$successful_replications, " & ", fmt(row$lambda_bias), " & ",
    fmt(row$lambda_rmse), " & ", fmt(row$multiplier_rmse), " \\\\"
  ))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(output_path, "main-core-lambda-recovery-table.tex"))
print(metrics)
