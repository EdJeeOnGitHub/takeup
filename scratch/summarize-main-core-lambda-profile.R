#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
source("R/structural/main-core-data.R")
manifest_path <- main_core_option_value(args, "--manifest")
output_path <- main_core_option_value(
  args, "--output-path", "temp-data/main-core-lambda-identification/profile"
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
treatment <- c("Control", "Ink", "Calendar", "Bracelet")
rows <- list()
multiplier_rows <- list()
for (index in seq_len(nrow(manifest))) {
  run_path <- file.path(output_path, "fits", manifest$label[index])
  status_path <- file.path(run_path, "status.csv")
  if (!file.exists(status_path)) next
  status <- read.csv(status_path, stringsAsFactors = FALSE)
  if (status$status[1] != "complete" || !file.exists(status$gq_csv[1])) next
  gq <- read_cmdstan_csv(status$gq_csv[1])
  draw <- as_draws_df(gq$generated_quantities)
  takeup_ll <- draw$core_compact_log_lik_takeup[1]
  beliefs_ll <- draw$core_compact_log_lik_beliefs[1]
  wtp_ll <- draw$core_compact_log_lik_wtp[1]
  distance_ll <- draw$core_compact_log_lik_distance[1]
  rows[[length(rows) + 1L]] <- data.frame(
    theta = status$theta[1], objective = status$lp[1],
    takeup_log_lik = takeup_ll, beliefs_log_lik = beliefs_ll,
    wtp_log_lik = wtp_ll, distance_log_lik = distance_ll,
    total_log_lik = takeup_ll + beliefs_ll + wtp_ll + distance_ll
  )
  for (treatment_index in 1:4) {
    multiplier_rows[[length(multiplier_rows) + 1L]] <- data.frame(
      theta = status$theta[1], treatment = treatment[treatment_index],
      multiplier_1500m = -draw[[sprintf(
        "core_compact_sm_rescaled[2,%d]", treatment_index
      )]][1]
    )
  }
}
profile <- bind_rows(rows) |> arrange(theta)
multipliers <- bind_rows(multiplier_rows)
if (nrow(profile) == 0L) stop("No completed profile points.", call. = FALSE)
profile <- profile |>
  mutate(
    objective_drop = max(objective) - objective,
    data_log_lik_drop = max(total_log_lik) - total_log_lik,
    support_2 = objective_drop <= 2,
    support_4 = objective_drop <= 4
  )
write.csv(profile, file.path(output_path, "lambda-profile-summary.csv"), row.names = FALSE)
write.csv(
  multipliers,
  file.path(output_path, "lambda-profile-multipliers.csv"),
  row.names = FALSE
)
support <- profile |>
  summarise(
    theta_hat = theta[which.max(objective)],
    support_2_low = min(theta[support_2]),
    support_2_high = max(theta[support_2]),
    support_4_low = min(theta[support_4]),
    support_4_high = max(theta[support_4]),
    left_endpoint_drop = first(objective_drop),
    right_endpoint_drop = last(objective_drop)
  )
write.csv(support, file.path(output_path, "lambda-profile-support.csv"), row.names = FALSE)

p1 <- ggplot(profile, aes(theta, objective_drop)) +
  geom_hline(yintercept = c(2, 4), linetype = c("solid", "dashed"), color = "grey55") +
  geom_line(size = 0.5) + geom_point(size = 0.8) +
  coord_cartesian(ylim = c(0, min(10, max(profile$objective_drop)))) +
  labs(
    x = expression(log(lambda[Signal] / lambda[No~Signal])),
    y = "Drop in regularized profile objective"
  ) + theme_bw(base_size = 9)
p2 <- ggplot(multipliers, aes(theta, multiplier_1500m, color = treatment)) +
  geom_hline(yintercept = 1, color = "grey65", size = 0.3) +
  geom_line(size = 0.55) +
  labs(
    x = expression(log(lambda[Signal] / lambda[No~Signal])),
    y = "Social multiplier at 1500m", color = NULL
  ) + theme_bw(base_size = 9) + theme(legend.position = "bottom")
pdf(file.path(output_path, "main-core-lambda-profile.pdf"), width = 7.2, height = 6.4)
print(p1)
print(p2)
dev.off()
print(support)
