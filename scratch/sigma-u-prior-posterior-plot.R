#!/usr/bin/env Rscript

library(tidyverse)
library(cmdstanr)
library(posterior)
library(ggthemes)

script_options <- docopt::docopt(
  "Usage:
  sigma-u-prior-posterior-plot.R [options]

Options:
  --input-path=<path>   Stan CSV input path [default: data/stan_analysis_data]
  --output-path=<path>  Figure output path [default: presentations/figures]
  --fit-version=<num>   Fit version [default: 95]
  --model=<model>       Structural model name [default: STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB]
",
  args = if (interactive()) "" else commandArgs(trailingOnly = TRUE)
)

read_sigma_u <- function(files) {
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop("Missing Stan CSV(s): ", paste(missing_files, collapse = ", "))
  }

  fit <- cmdstanr::read_cmdstan_csv(files, variables = "u_sd")

  posterior::as_draws_df(fit$post_warmup_draws) %>%
    transmute(sigma_u = `u_sd[1]`) %>%
    filter(!is.na(sigma_u))
}

fit_files <- file.path(
  script_options$input_path,
  str_glue("dist_fit{script_options$fit_version}_{script_options$model}-{1:4}.csv")
)

prior_files <- file.path(
  script_options$input_path,
  str_glue("dist_prior{script_options$fit_version}_{script_options$model}-{1:4}.csv")
)

sigma_u_draws <- bind_rows(
  read_sigma_u(prior_files) %>% mutate(distribution = "Prior"),
  read_sigma_u(fit_files) %>% mutate(distribution = "Posterior")
) %>%
  mutate(distribution = factor(distribution, levels = c("Prior", "Posterior")))

summary_df <- sigma_u_draws %>%
  group_by(distribution) %>%
  summarise(
    mean = mean(sigma_u),
    median = median(sigma_u),
    q025 = quantile(sigma_u, 0.025),
    q975 = quantile(sigma_u, 0.975),
    .groups = "drop"
  )

dir.create(script_options$output_path, recursive = TRUE, showWarnings = FALSE)

plot_file <- file.path(script_options$output_path, "sigma-u-prior-posterior.pdf")
summary_file <- file.path(script_options$output_path, "sigma-u-prior-posterior-summary.csv")

canva_palette_vibrant <- "Primary colors with a vibrant twist"

sigma_u_plot <- ggplot(sigma_u_draws, aes(x = sigma_u, color = distribution, fill = distribution)) +
  geom_density(alpha = 0.28, linewidth = 0.9, adjust = 1.1) +
  geom_vline(
    data = summary_df,
    aes(xintercept = median, color = distribution),
    linetype = "dashed",
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  ggthemes::scale_color_canva(palette = canva_palette_vibrant) +
  ggthemes::scale_fill_canva(palette = canva_palette_vibrant) +
  coord_cartesian(xlim = quantile(sigma_u_draws$sigma_u, c(0, 0.995), na.rm = TRUE)) +
  labs(
    x = expression(sigma[u]),
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(plot_file, sigma_u_plot, width = 5.8, height = 3.6)
readr::write_csv(summary_df, summary_file)

message("Wrote ", plot_file)
message("Wrote ", summary_file)
