#!/usr/bin/Rscript
# structural_tables.R
# Standalone script: generates structural model tables (overall ATEs,
# signal/private decomposition, belief ATEs) from quick-postprocess RDS output.
#
# Usage:
#   Rscript --no-save --no-restore structural_tables.R \
#     --fit-version=104 \
#     --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
#
# Requires quick_roc_postprocess.R --sm to have been run first (for sm-decomp figure).
#
# Outputs (to presentations/tables/fit<VERSION>/):
#   struct-overall-te-table.tex
#   private-signal-te-table.tex
#   fob-beliefs-table.tex
#   indiv-dist-community-fp-indiv-vis-robust-struct-overall-te-table.tex  [--write-robustness]
#   indiv-dist-indiv-fp-robust-struct-overall-te-table.tex                [--write-robustness]
#   struct-robustness-nooutliers-overall-te-table.tex                     [--write-robustness]
#
# Outputs (to presentations/figures/fit<VERSION>/):
#   sm-decomp-annotated.pdf

script_options <- docopt::docopt(
  stringr::str_glue(
"Usage:
  structural_tables.R [options]

Options:
  --fit-version=<v>       Fit version number [default: 104]
  --model=<m>             Structural model name [default: STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP]
  --input-path=<path>     Path to Stan analysis data [default: temp-data/struct-postprocess]
  --output-path=<path>    Path to postprocessed RDS files [default: temp-data/struct-postprocess]
  --table-output=<path>   Path to write .tex tables [default: presentations/tables/fit<VERSION>]
  --figure-output=<path>  Path to write figures [default: presentations/figures/fit<VERSION>]
  --width=<w>             Credible interval width [default: 0.95]
  --write-robustness      Also write appendix robustness tables
  "),
  args = if (interactive()) "--fit-version=105 --write-robustness" else commandArgs(trailingOnly = TRUE)
)

library(tidyverse)
library(posterior)
library(tidybayes)
library(knitr)
library(kableExtra)
library(magrittr)
library(stringr)
library(ggthemes)

options(knitr.kable.NA = '')  

fit_version_int <- as.integer(script_options$fit_version)
params <- list(
  fit_version      = fit_version_int,
  struct_models    = script_options$model,
  input_path       = script_options$input_path,
  output_path      = script_options$output_path,
  table_output_path = if (script_options$table_output == "presentations/tables/fit<VERSION>")
    str_glue("presentations/tables/fit{fit_version_int}")
  else
    script_options$table_output,
  figure_output_path = if (script_options$figure_output == "presentations/figures/fit<VERSION>")
    str_glue("presentations/figures/fit{fit_version_int}")
  else
    script_options$figure_output,
  width            = as.numeric(script_options$width)
)

cat(str_glue("fit_version:   {params$fit_version}\n"))
cat(str_glue("struct_models: {params$struct_models}\n"))
cat(str_glue("input_path:    {params$input_path}\n"))
cat(str_glue("output_path:   {params$output_path}\n"))
cat(str_glue("table_output:  {params$table_output_path}\n\n"))

dir.create(params$table_output_path, showWarnings = FALSE, recursive = TRUE)
dir.create(params$figure_output_path, showWarnings = FALSE, recursive = TRUE)
dir.create(params$output_path, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load analysis data
# ---------------------------------------------------------------------------
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path("dist_structural_util.R"))
source(file.path("multilvlr", "multilvlr_util.R"))

standardize   <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

load(file.path("data", "takeup_village_pot_dist.RData"))
load(file.path("data", "analysis.RData"))

monitored_nosms_data <- analysis.data %>%
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>%
  left_join(
    village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
    by = "cluster.id"
  ) %>%
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>%
  group_by(cluster.id) %>%
  mutate(cluster_id = cur_group_id()) %>%
  ungroup()

analysis_data <- monitored_nosms_data

# ---------------------------------------------------------------------------
# Helper: save latex table (strips \begin{table}...\end{table} wrapper)
# ---------------------------------------------------------------------------
custom_save_latex_table <- function(table, table_name) {
  table_conn <- file(file.path(params$table_output_path, paste0(table_name, ".tex")))
  attr(table, "kable_meta")$contents <-
    str_replace_all(attr(table, "kable_meta")$contents, "removeme12345", " ")
  table[1] <- str_replace_all(table[1], "removeme12345", " ")
  clean_table <- table %>%
    str_remove(., fixed("\\begin{table}")) %>%
    str_remove(., "\\\\caption\\{.*\\}") %>%
    str_remove(., "\\\\end\\{table\\}")
  clean_table %>% writeLines(table_conn)
  close(table_conn)
  invisible(table)
}

create_robustness_tbl <- function(data) {
  data %>%
    kbl(
      col.names = c("Dependent variable: Take-up", paste0("(", 1:4, ")")),
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      align = "lcccc",
      caption = "Overall Results"
    ) %>%
    kable_styling(latex_options = c("scale_down")) %>%
    add_header_above(
      c(" ", "Combined", "Close", "Far", "Far - Close"),
      line = FALSE
    ) %>%
    add_header_above(c(" " = 1, "Structural" = 4)) %>%
    row_spec(c(3), hline_after = TRUE)
}

# ---------------------------------------------------------------------------
# Load structural ATE RDS (no RF data needed for structural-only tables)
# ---------------------------------------------------------------------------
ate_rvar_struct <- read_rds(
  file.path(
    params$input_path,
    str_glue("rvar_processed_dist_fit{params$fit_version}_ates_{params$struct_models}_1-4.rds")
  )
)

# ate_rvar_df contains only structural rows; spread_rf() will produce only
# struct_* columns, and the downstream select(-contains("rf_")) is a no-op.
ate_rvar_df <- ate_rvar_struct

# ---------------------------------------------------------------------------
# Helper functions (adapted from tables.Rmd)
# ---------------------------------------------------------------------------
create_cis <- function(.data, .width = 0.95) {
  med_fun <- function(x) {
    mean_x   <- mean(x) %>% round(3)
    conf.low  <- quantile(x, (1 - .width) / 2) %>% round(3)
    conf.high <- quantile(x, 1 - (1 - .width) / 2) %>% round(3)
    linebreak(
      paste0(mean_x, "\n", str_glue("({conf.low}, {conf.high})")),
      align = "c"
    )
  }
  .data %>% mutate(across(where(is_rvar), med_fun))
}

create_ate_table <- function(.data, .estimand, group_var = treatment) {
  .data %>%
    filter(estimand == .estimand) %>%
    select(model, {{ group_var }}, dist_group, value) %>%
    pivot_wider(names_from = dist_group, values_from = value) %>%
    select(model, {{ group_var }}, any_of(c("combined", "close", "far"))) %>%
    arrange({{ group_var }}) %>%
    bind_rows(
      # bracelet minus calendar row
      .data %>%
        filter({{ group_var }} %in% c("Bracelet", "Calendar")) %>%
        filter(estimand == .estimand) %>%
        pivot_wider(names_from = {{ group_var }}, values_from = value) %>%
        mutate(bracelet_minus_calendar = Bracelet - Calendar) %>%
        select(model, dist_group, bracelet_minus_calendar) %>%
        pivot_wider(names_from = dist_group, values_from = bracelet_minus_calendar) %>%
        mutate("{{ group_var }}" := "bracelet_minus_calendar"),
      # signal minus no-signal row
      .data %>%
        filter(estimand == .estimand) %>%
        pivot_wider(names_from = dist_group, values_from = value) %>%
        select(model, treatment, any_of(c("combined", "close", "far"))) %>%
        mutate(signal = case_when(
          treatment %in% c("Bracelet", "Ink") ~ "Signal",
          treatment %in% c("Calendar")         ~ "No Signal",
          TRUE ~ NA
        )) %>%
        filter(!is.na(signal)) %>%
        group_by(model, signal) %>%
        summarise(across(where(is_rvar), rvar_mean)) %>%
        group_by(model) %>%
        summarise(
          across(where(is_rvar), ~ .x[signal == "Signal"] - .x[signal == "No Signal"])
        ) %>%
        mutate("{{ group_var }}" := "signal_minus_no_signal")
    )
}

ate_tbl_levels <- c(
  "Bracelet", "Calendar", "Ink", "Control",
  "signal_minus_no_signal", "bracelet_minus_calendar"
)

recode_control_mean <- function(data) {
  data %>% mutate(
    treatment = fct_recode(
      treatment,
      "Control mean"       = "Control",
      "$\\Delta$ (Bracelet - Calendar)" = "bracelet_minus_calendar",
      "$\\Delta$ (Signal - No Signal)"  = "signal_minus_no_signal"
    )
  )
}

spread_rf <- function(data) {
  data %>%
    mutate(model_type = case_when(
      model == params$struct_models ~ "struct",
      TRUE ~ "rf"
    )) %>%
    select(-model) %>%
    pivot_wider(
      names_from  = model_type,
      id_cols     = treatment,
      names_glue  = "{model_type}_{.value}",
      values_from = any_of(c("combined", "close", "far", "far_minus_close"))
    ) %>%
    select(treatment, starts_with("rf"), starts_with("struct"))
}

prep_robustness <- function(data) {
  data %>%
    filter(variable != "far_minus_close") %>%
    filter(treatment != "bracelet_minus_calendar") %>%
    rename(dist_group = variable) %>%
    create_ate_table(.estimand = "overall", group_var = treatment) %>%
    mutate(treatment = factor(treatment, levels = ate_tbl_levels)) %>%
    mutate(far_minus_close = far - close) %>%
    arrange(treatment) %>%
    create_cis(.width = params$width) %>%
    recode_control_mean() %>%
    mutate(estimand = "overall") %>%
    select(-model, -estimand)
}

write_robustness_tables <- function() {
  cat("Writing robustness appendix tables ...\n")

  robustness_models <- tribble(
    ~model, ~table_name, ~is_indiv,
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_COMMUNITY_FP_INDIV_VIS",
    "indiv-dist-community-fp-indiv-vis-robust-struct-overall-te-table",
    TRUE,
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP",
    "indiv-dist-indiv-fp-robust-struct-overall-te-table",
    TRUE,
    "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_NO_OUTLIERS",
    "struct-robustness-nooutliers-overall-te-table",
    FALSE
  )

  for (i in seq_len(nrow(robustness_models))) {
    robustness_row <- robustness_models[i, ]
    rds_name <- if (robustness_row$is_indiv) {
      str_glue(
        "rvar_processed_dist_fit{params$fit_version}_INDIV_incentive_tes_{robustness_row$model}_1-4.rds"
      )
    } else {
      str_glue(
        "rvar_processed_dist_fit{params$fit_version}_ates_{robustness_row$model}_1-4.rds"
      )
    }

    robust_rvar <- read_rds(file.path(params$input_path, rds_name))

    robust_tbl <- if (robustness_row$is_indiv) {
      prep_robustness(robust_rvar)
    } else {
      robust_rvar %>%
        create_ate_table(.estimand = "overall", group_var = treatment) %>%
        mutate(treatment = factor(treatment, levels = ate_tbl_levels)) %>%
        mutate(far_minus_close = far - close) %>%
        arrange(treatment) %>%
        create_cis(.width = params$width) %>%
        recode_control_mean() %>%
        mutate(estimand = "overall") %>%
        select(-model, -estimand)
    }

    robust_tbl %>%
      create_robustness_tbl() %>%
      custom_save_latex_table(robustness_row$table_name)
  }
}

# ---------------------------------------------------------------------------
# Build ATE data frames
# ---------------------------------------------------------------------------
incentive_ate_df <- ate_rvar_df %>%
  create_ate_table(.estimand = "overall", group_var = treatment)

signal_ate_df <- ate_rvar_df %>%
  filter(model == params$struct_models) %>%
  mutate(mu_treatment = str_to_title(mu_treatment)) %>%
  create_ate_table(.estimand = "signal", group_var = mu_treatment)

private_ate_df <- ate_rvar_df %>%
  filter(model == params$struct_models) %>%
  create_ate_table(.estimand = "private", group_var = treatment)

incentive_ate_tbl <- incentive_ate_df %>%
  mutate(treatment = factor(treatment, levels = ate_tbl_levels)) %>%
  mutate(far_minus_close = far - close) %>%
  arrange(treatment) %>%
  create_cis(.width = params$width) %>%
  recode_control_mean() %>%
  spread_rf() %>%
  mutate(estimand = "overall")

private_ate_tbl <- private_ate_df %>%
  mutate(treatment = factor(treatment, levels = ate_tbl_levels)) %>%
  arrange(treatment) %>%
  create_cis(.width = params$width) %>%
  recode_control_mean() %>%
  spread_rf() %>%
  mutate(estimand = "private") %>%
  filter(treatment %in% c("Bracelet", "Calendar", "Ink"))

signal_ate_tbl <- signal_ate_df %>%
  rename(treatment = mu_treatment) %>%
  mutate(treatment = factor(treatment, levels = ate_tbl_levels)) %>%
  mutate(far_minus_close = far - close) %>%
  arrange(treatment) %>%
  create_cis(.width = params$width) %>%
  recode_control_mean() %>%
  spread_rf() %>%
  mutate(estimand = "signal") %>%
  filter(treatment %in% c("Bracelet", "Calendar", "Ink"))

all_ate_tbl <- bind_rows(incentive_ate_tbl, signal_ate_tbl, private_ate_tbl) %>%
  mutate(estimand = factor(estimand, levels = c("private", "signal", "overall")) %>% fct_rev) %>%
  select(estimand, treatment, everything())

# ---------------------------------------------------------------------------
# Table 1: Structural overall ATEs
# ---------------------------------------------------------------------------
cat("Writing struct-overall-te-table.tex ...\n")

struct_overall_tbl <- all_ate_tbl %>%
  arrange(estimand, treatment) %>%
  filter(estimand == "overall") %>%
  select(-estimand) %>%
  select(-contains("rf_")) %>%
  kbl(
    col.names = c("Dependent variable: Take-up", paste0("(", 1:4, ")")),
    format = "latex",
    linesep = "\\addlinespace ",
    booktabs = TRUE, escape = FALSE, align = "lcccc",
    caption = "Overall Results"
  ) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  add_header_above(
    c(" ", "Combined", "Close", "Far", "Far - Close"),
    line = FALSE
  ) %>%
  add_header_above(c(" " = 1, "Structural" = 4)) %>%
  row_spec(c(3), hline_after = TRUE)

struct_overall_tbl %>% custom_save_latex_table("struct-overall-te-table")

# ---------------------------------------------------------------------------
# Table 2: Signal + Private decomposition
# ---------------------------------------------------------------------------
cat("Writing private-signal-te-table.tex ...\n")

decomposed_te_kbl_df <- all_ate_tbl %>%
  arrange(estimand, treatment) %>%
  filter(estimand %in% c("signal", "private")) %>%
  select(estimand, treatment, contains("struct_")) %>%
  select(-estimand) %>%
  kbl(
    col.names = c("Dependent variable: Take-up", paste0("(", 1:4, ")")),
    format = "latex",
    linesep = "\\addlinespace ",
    booktabs = TRUE, escape = FALSE, align = "lcccc",
    caption = "Signal and Private Value"
  ) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  add_header_above(
    c(" ", "Combined", "Close", "Far", "Far - Close"),
    line = FALSE
  ) %>%
  add_header_above(c(" " = 1, "Structural" = 4)) %>%
  pack_rows(
    index = c("Panel A: Social Image" = 3, "Panel B: Private" = 3),
    italic = TRUE, escape = FALSE,
    hline_after = TRUE, hline_before = TRUE, bold = TRUE
  )

decomposed_te_kbl_df %>% custom_save_latex_table("private-signal-te-table")

# ---------------------------------------------------------------------------
# Table 3: Belief ATEs (first-order beliefs)
# ---------------------------------------------------------------------------
cat("Writing fob-beliefs-table.tex ...\n")

obs_rvar_struct <- read_rds(
  file.path(
    params$input_path,
    str_glue("rvar_processed_dist_fit{params$fit_version}_belief_ates_{params$struct_models}_1-4.rds")
  )
)

lvl_obs_rvar_struct <- read_rds(
  file.path(
    params$input_path,
    str_glue("rvar_processed_dist_fit{params$fit_version}_belief_probs_{params$struct_models}_1-4.rds")
  )
) %>%
  as.data.frame() %>%
  filter(treatment == "Control") %>%
  filter(variable == "prob_1ord") %>%
  as_tibble() %>%
  select(-dist_treat_idx)

obs_rvar_struct <- bind_rows(
  obs_rvar_struct %>% as.data.frame(),
  lvl_obs_rvar_struct %>% as.data.frame()
) %>% as_tibble()

create_ate_kbl <- function(data, y_var_name = "Dependent variable: Observability") {
  data %>%
    kbl(
      col.names = c(y_var_name, paste0("(", 1:4, ")")),
      format = "latex", booktabs = TRUE, escape = FALSE, align = "lcccc",
      caption = "Belief Results"
    ) %>%
    kable_styling(latex_options = c("scale_down")) %>%
    add_header_above(
      c(" ", "Combined", "Close", "Far", "Far - Close"),
      line = FALSE
    ) %>%
    add_header_above(c(" " = 1, "Structural" = 4)) %>%
    row_spec(c(3), hline_after = TRUE)
}

obs_rvar_df <- obs_rvar_struct %>%
  as.data.frame() %>%
  filter(variable == "ate_1ord" | variable == "prob_1ord") %>%
  as_tibble() %>%
  mutate(estimand = "overall") %>%
  create_ate_table(.estimand = "overall", group_var = treatment)

obs_rvar_tbl <- obs_rvar_df %>%
  mutate(treatment = factor(treatment, levels = ate_tbl_levels)) %>%
  mutate(far_minus_close = far - close) %>%
  arrange(treatment) %>%
  create_cis(.width = params$width) %>%
  recode_control_mean() %>%
  select(-model) %>%
  create_ate_kbl()

obs_rvar_tbl %>% custom_save_latex_table("fob-beliefs-table")

if (isTRUE(script_options$write_robustness)) {
  write_robustness_tables()
}

# ---------------------------------------------------------------------------
# Figure: Social Multiplier Decomposition (sm-decomp-annotated)
# Requires: quick_roc_postprocess.R --sm has been run first.
# ---------------------------------------------------------------------------
cat("Writing sm-decomp-annotated.pdf ...\n")

sm_summ_df <- read_rds(
  file.path(
    params$input_path,
    str_glue("rvar_processed_dist_fit{params$fit_version}_sm_summ_{params$struct_models}_1-4.rds")
  )
)

# Reuse fit-95 prior predictive (model identical; no prior run for fit 105).
# Use tidy_processed (not rvar_processed): sm_rescaled is NA in the rvar file.
prior_sm_path <- "temp-data/tidy_processed_dist_prior95_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB_1-4.rds"
if (!file.exists(prior_sm_path)) {
  warning("Prior predictive tidy sm_draws not found at: ", prior_sm_path,
          "\nSkipping prior predictive line in sm-decomp plot.")
  prior_sm_summ_df <- NULL
} else {
  prior_tidy <- read_rds(prior_sm_path)
  prior_sm_summ_df <- prior_tidy[prior_tidy$param == "sm_draws", ]$tidy_draws[[1]] %>%
    filter(variable == "sm_rescaled", roc_distance <= 2500, .width == 0.95) %>%
    mutate(across(c(value, conf.low, conf.high), ~ . * -1))
}

canva_palette_vibrant = "Primary colors with a vibrant twist"

theme_set(theme_minimal() +
            theme(legend.position = "bottom"))

col_df <- tibble(
  variable = c("Social Multiplier", "Social Multiplier - Type Inference Only", "Prior Predictive"),
  colour   = c(canva_pal(canva_palette_vibrant)(2)[1:2], "grey")
)

sm_plot_df <- sm_summ_df %>%
  filter(variable %in% c("sm_rescaled", "sm_delta_part_rescaled")) %>%
  mutate(
    variable_label = case_when(
      variable == "sm_rescaled"            ~ "Social Multiplier",
      variable == "sm_delta_part_rescaled" ~ "Social Multiplier - Type Inference Only"
    ),
    variable_label = factor(variable_label, levels = col_df$variable),
    variable       = factor(variable, levels = c("sm_rescaled", "sm_delta_part_rescaled"))
  ) %>%
  filter(roc_distance <= 2500) %>%
  mutate(across(c(value, conf.low, conf.high), ~ . * -1))

sm_plot_df %>%
  filter(fit_type == "fit") %>%
  select(
    treatment,
    distance = roc_distance,
    variable_label,
    value
  ) %>%
  pivot_wider(
    names_from = variable_label,
    values_from = value
  ) %>%
  write_csv("temp-data/social-multiplier-decomposition-values.csv")

if (!is.null(prior_sm_summ_df)) {
  sm_plot_df <- bind_rows(
    sm_plot_df,
    prior_sm_summ_df %>%
      mutate(
        treatment      = as.character(treatment),
        variable_label = factor("Prior Predictive", levels = col_df$variable),
        variable       = factor(variable, levels = c("sm_rescaled", "sm_delta_part_rescaled"))
      )
  )
}

sm_plot_df <- sm_plot_df %>%
  mutate(treatment = factor(treatment, levels = c("Control", "Ink", "Calendar", "Bracelet")))

sm_decomp_plot <- sm_plot_df %>%
  ggplot(aes(
    x        = roc_distance / 1000,
    y        = value,
    linetype = variable,
    colour   = variable_label
  )) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~treatment) +
  geom_hline(yintercept = 1, linetype = "longdash") +
  scale_color_manual("", values = col_df$colour, labels = col_df$variable) +
  guides(linetype = "none") +
  labs(x = "Distance [km]", y = "Social Multiplier") +
  annotate("text", x = 0.3,  y = 1.05, label = "Amplification", size = 3, alpha = 0.7) +
  annotate("text", x = 2.25, y = 0.95, label = "Mitigation",    size = 3, alpha = 0.7)


ggsave(
  file.path(params$figure_output_path, "sm-decomp-annotated.pdf"),
  plot   = sm_decomp_plot,
  width  = 8,
  height = 6
)
sm_decomp_plot
cat("\nDone. Tables written to:", params$table_output_path, "\n")
cat("Figures written to:", params$figure_output_path, "\n")
