#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
option_value <- function(name, default = NULL) {
  main_core_option_value(args, name, default)
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(ggplot2)
})

input_path <- option_value(
  "--input-path", "temp-data/main-core-benchmark-recovery"
)
table_path <- option_value(
  "--table-path",
  "appendix/structural-robustness/tables/main-core-benchmark-recovery.tex"
)
figure_path <- option_value(
  "--figure-path",
  "appendix/structural-robustness/figures/main-core-benchmark-likelihood.pdf"
)
expected_replicates <- as.integer(option_value("--expected-replicates", "50"))

manifest_file <- file.path(input_path, "benchmark-recovery-manifest.csv")
truth_file <- file.path(input_path, "truth-values.csv")
candidate_file <- file.path(input_path, "likelihood-candidate-manifest.csv")
if (!all(file.exists(c(manifest_file, truth_file, candidate_file)))) {
  stop("Benchmark-recovery preparation files are incomplete.")
}
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
truth <- read.csv(truth_file, stringsAsFactors = FALSE)
candidates <- read.csv(candidate_file, stringsAsFactors = FALSE)
if (nrow(manifest) != expected_replicates) {
  stop("Expected ", expected_replicates, " datasets; found ", nrow(manifest), ".")
}

labels <- c(
  delta = "Distance cost ($\\delta$)",
  lambda = "Social-image weight ($\\lambda$)",
  sigma_u = "Decision noise ($\\sigma_u$)",
  psi = "Private-value mapping ($\\psi$)",
  multiplier_control = "Control", multiplier_ink = "Ink",
  multiplier_calendar = "Calendar", multiplier_bracelet = "Bracelet",
  contrast_no_signal_minus_any_signal = "No Signal $-$ Any Signal"
)

recovery_rows <- list()
likelihood_rows <- list()
audit_rows <- list()
for (index in seq_len(nrow(manifest))) {
  task <- manifest[index, ]
  run_root <- file.path(input_path, "fits", sprintf("rep-%03d", task$replicate))
  status_file <- file.path(run_root, "status.csv")
  if (!file.exists(status_file)) stop("Missing status for replicate ", task$replicate)
  status <- read.csv(status_file, stringsAsFactors = FALSE)
  if (status$status[1] != "complete") stop("Incomplete replicate ", task$replicate)
  fit_csvs <- strsplit(status$fit_csvs[1], ";", fixed = TRUE)[[1L]]
  gq_csvs <- strsplit(status$gq_csvs[1], ";", fixed = TRUE)[[1L]]
  required <- c(fit_csvs, gq_csvs, status$likelihood_gq[1])
  if (any(!file.exists(required))) stop("Missing output for replicate ", task$replicate)

  fit_draws <- as_draws_df(read_cmdstan_csv(fit_csvs)$post_warmup_draws)
  gq_draws <- as_draws_df(read_cmdstan_csv(gq_csvs)$generated_quantities)
  multiplier_matrix <- -as.matrix(gq_draws[, paste0(
    "core_compact_sm_rescaled[1,", 1:4, "]"
  )])
  for (truth_index in seq_len(nrow(truth))) {
    item <- truth[truth_index, ]
    values <- if (item$object == "contrast_no_signal_minus_any_signal") {
      (multiplier_matrix[, 1] + multiplier_matrix[, 3] -
       multiplier_matrix[, 2] - multiplier_matrix[, 4]) / 2
    } else if (startsWith(item$object, "multiplier_")) {
      -as.numeric(gq_draws[[item$variable]])
    } else {
      as.numeric(fit_draws[[item$variable]])
    }
    recovery_rows[[length(recovery_rows) + 1L]] <- data.frame(
      replicate = task$replicate, object = item$object, truth = item$truth,
      posterior_mean = mean(values), posterior_median = median(values),
      lower_50 = quantile(values, 0.25), upper_50 = quantile(values, 0.75),
      lower_80 = quantile(values, 0.10), upper_80 = quantile(values, 0.90),
      lower_90 = quantile(values, 0.05), upper_90 = quantile(values, 0.95),
      lower_95 = quantile(values, 0.025), upper_95 = quantile(values, 0.975),
      stringsAsFactors = FALSE
    )
  }

  final_audit <- file.path(
    run_root, paste0("attempt-", status$final_attempt[1]), "audit.csv"
  )
  audit <- read.csv(final_audit, stringsAsFactors = FALSE)
  audit$replicate <- task$replicate
  audit_rows[[length(audit_rows) + 1L]] <- audit

  ll_draws <- as_draws_df(
    read_cmdstan_csv(status$likelihood_gq[1])$generated_quantities
  )
  if (nrow(ll_draws) != nrow(candidates)) {
    stop("Likelihood grid row mismatch for replicate ", task$replicate)
  }
  likelihood_rows[[length(likelihood_rows) + 1L]] <- data.frame(
    replicate = task$replicate, candidates,
    takeup = ll_draws$core_compact_log_lik_takeup,
    beliefs = ll_draws$core_compact_log_lik_beliefs,
    wtp = ll_draws$core_compact_log_lik_wtp,
    multiplier_control = -ll_draws$`core_compact_sm_rescaled[1,1]`,
    multiplier_ink = -ll_draws$`core_compact_sm_rescaled[1,2]`,
    multiplier_calendar = -ll_draws$`core_compact_sm_rescaled[1,3]`,
    multiplier_bracelet = -ll_draws$`core_compact_sm_rescaled[1,4]`,
    stringsAsFactors = FALSE
  )
}

recovery <- bind_rows(recovery_rows) |>
  mutate(
    covered_50 = lower_50 <= truth & upper_50 >= truth,
    covered_80 = lower_80 <= truth & upper_80 >= truth,
    covered_90 = lower_90 <= truth & upper_90 >= truth,
    covered_95 = lower_95 <= truth & upper_95 >= truth,
    median_error = posterior_median - truth,
    mean_error = posterior_mean - truth,
    width_95 = upper_95 - lower_95,
    direction_recovered = case_when(
      object %in% c("multiplier_control", "multiplier_calendar") ~ posterior_median > 1,
      object %in% c("multiplier_ink", "multiplier_bracelet") ~ posterior_median < 1,
      object == "contrast_no_signal_minus_any_signal" ~ posterior_median > 0,
      TRUE ~ NA
    )
  )
metrics <- recovery |>
  group_by(object, truth) |>
  summarise(
    replications = n(), mean_posterior_median = mean(posterior_median),
    median_bias = median(median_error), mean_bias = mean(mean_error),
    rmse = sqrt(mean(median_error^2)),
    coverage_50 = mean(covered_50), coverage_80 = mean(covered_80),
    coverage_90 = mean(covered_90), coverage_95 = mean(covered_95),
    mean_width_95 = mean(width_95),
    direction_recovery = mean(direction_recovered, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(label = unname(labels[object]))

likelihood <- bind_rows(likelihood_rows) |>
  mutate(total = takeup + beliefs + wtp) |>
  group_by(replicate, object) |>
  mutate(
    truth_row = which.min(abs(factor - 1))[1],
    total_difference_from_truth = total - total[truth_row],
    multiplier_difference = pmax(
      abs(multiplier_control - multiplier_control[truth_row]),
      abs(multiplier_ink - multiplier_ink[truth_row]),
      abs(multiplier_calendar - multiplier_calendar[truth_row]),
      abs(multiplier_bracelet - multiplier_bracelet[truth_row])
    )
  ) |>
  ungroup()
likelihood_summary <- likelihood |>
  group_by(object, candidate_id, factor, value, truth) |>
  summarise(
    mean_log_lik_difference = mean(total_difference_from_truth),
    mcse = ifelse(n() > 1, sd(total_difference_from_truth) / sqrt(n()), 0),
    max_multiplier_difference = max(multiplier_difference), .groups = "drop"
  ) |>
  group_by(object) |>
  mutate(
    log_lik_difference = mean_log_lik_difference - max(mean_log_lik_difference),
    maximizing_factor = factor[which.max(mean_log_lik_difference)][1]
  ) |>
  ungroup()

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
write.csv(recovery, file.path(input_path, "benchmark-recovery-estimates.csv"),
          row.names = FALSE)
write.csv(metrics, file.path(input_path, "benchmark-recovery-metrics.csv"),
          row.names = FALSE)
write.csv(bind_rows(audit_rows), file.path(input_path, "benchmark-recovery-hmc-audit.csv"),
          row.names = FALSE)
write.csv(likelihood_summary,
          file.path(input_path, "benchmark-recovery-likelihood-slices.csv"),
          row.names = FALSE)
scale_invariance <- likelihood_summary |>
  filter(object == "raw_scale_path") |>
  transmute(
    factor, log_lik_difference = mean_log_lik_difference,
    max_multiplier_difference
  )
write.csv(scale_invariance,
          file.path(input_path, "benchmark-recovery-scale-invariance.csv"),
          row.names = FALSE)

fmt <- function(x) sprintf("%.3f", x)
table_lines <- c(
  "\\begin{tabular}{lrrrrr}", "\\toprule",
  "Multiplier object & Truth & Mean estimate & Bias & RMSE & 95\\% coverage \\\\",
  "\\midrule"
)
for (object in c(
  paste0("multiplier_", c("control", "ink", "calendar", "bracelet")),
  "contrast_no_signal_minus_any_signal"
)) {
  row <- metrics[metrics$object == object, ]
  table_lines <- c(table_lines, paste0(
    row$label, " & ", fmt(row$truth), " & ", fmt(row$mean_posterior_median),
    " & ", fmt(row$median_bias), " & ", fmt(row$rmse), " & ",
    sprintf("%.2f", row$coverage_95), " \\\\"
  ))
}
table_lines <- c(table_lines, "\\bottomrule", "\\end{tabular}")
writeLines(table_lines, table_path)

plot_labels <- c(
  normalized_delta = "Normalized distance cost",
  normalized_lambda = "Normalized social-image weight",
  raw_scale_path = "Raw utility-scale ridge"
)
p <- ggplot(likelihood_summary, aes(factor, log_lik_difference)) +
  geom_ribbon(aes(
    ymin = log_lik_difference - 1.96 * mcse,
    ymax = log_lik_difference + 1.96 * mcse
  ), fill = "grey80") +
  geom_line(linewidth = 0.55, color = "#1f4e79") +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45) +
  facet_wrap(~object, scales = "free_y", labeller = as_labeller(plot_labels)) +
  labs(
    x = "Candidate value relative to generating truth",
    y = "Average log likelihood relative to maximum"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())
ggsave(figure_path, p, width = 7.2, height = 4.8)

maxima <- likelihood_summary |>
  filter(object != "raw_scale_path") |>
  distinct(object, maximizing_factor)
if (any(abs(maxima$maximizing_factor - 1) > 1 / (length(unique(candidates$factor)) - 1) + 1e-8)) {
  warning("At least one average likelihood slice peaks more than one grid step from truth.")
}
if (any(!is.finite(scale_invariance$log_lik_difference)) ||
    any(!is.finite(scale_invariance$max_multiplier_difference)) ||
    any(abs(scale_invariance$log_lik_difference) > 1e-5) ||
    any(scale_invariance$max_multiplier_difference > 1e-5)) {
  stop("The raw-scale invariance audit failed tolerance.")
}
message("Wrote benchmark recovery table and likelihood figure.")
