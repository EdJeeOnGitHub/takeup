#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (!length(hit)) default else substring(hit[[1L]], nchar(name) + 2L)
}
specification <- value("--spec", "assigned")
source("R/distance/spec.R")
specification <- takeup_distance_spec(specification)
output_root <- value("--output-root", file.path("build", "candidate-components", specification))
figure_root <- file.path(output_root, "figures")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

sms_source <- file.path("build", "work", specification, "temp-data", "p-sms-tes.pdf")
if (!file.exists(sms_source)) stop("Run the reduced-form target before auxiliary rendering.")
file.copy(sms_source, file.path(figure_root, "sms-TE-by-dist-incentive.pdf"),
          overwrite = TRUE)

ri_path <- file.path("build", specification, "balance", "balance-ri-fe.rds")
if (!file.exists(ri_path)) stop("Run the balance fit-ri section before auxiliary rendering.")
ri <- readRDS(ri_path)
required <- c("plot_perm_fit_df", "ri_p_val_df", "realised_fit_df")
if (!all(required %in% names(ri))) stop("Unexpected balance RI object.")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
clean <- function(data) data |>
  mutate(clean_name = case_when(
    lhs == "lhs: age.census" ~ "Age",
    lhs == "lhs: ethnicity_luhya" ~ "Main ethnicity/Luhya",
    lhs == "lhs: have_phone_lgl" ~ "Phone owner",
    lhs == "lhs: n_per_cluster" ~ "Number of individuals per community",
    lhs == "lhs: religion_christianity" ~ "Christian",
    TRUE ~ clean_name
  )) |>
  filter(!is.na(clean_name),
         !clean_name %in% c("Distance to PoT", "Dewormed in the past year"))
draws <- clean(ri$plot_perm_fit_df)
labels <- clean(ri$ri_p_val_df)
plot <- ggplot(draws, aes(x = statistic)) +
  geom_histogram(color = "grey", linewidth = 0.1) +
  facet_wrap(~ clean_name, scales = "free", labeller = label_wrap_gen(width = 20)) +
  geom_vline(data = labels, aes(xintercept = realised_statistic)) +
  geom_label(data = labels, aes(x = x, y = Inf, label = p_val), vjust = 1.2) +
  theme_bw() + labs(x = "t statistic", y = "count")
ggsave(file.path(figure_root, "dist-ri-plot.pdf"), plot, width = 10, height = 10)
utils::write.csv(data.frame(
  distance_definition = specification,
  sms_source = sms_source, ri_source = ri_path,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output_root, "auxiliary-rf-manifest.csv"), row.names = FALSE)
