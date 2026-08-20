library(tidyverse)
library(fixest)

set.seed(20260710)

params <- lst(
  table_output_path = "presentations/rf-tables/main-specs",
  output_path = "temp-data",
  n_sims = 1000
)

source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source("R/common/analysis.R")
source("R/structural/legacy-utils.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source(file.path("scratch", "reduced-form-setup.R"))

dir.create(params$table_output_path, showWarnings = FALSE, recursive = TRUE)
dir.create(params$output_path, showWarnings = FALSE, recursive = TRUE)

summ_know_A_df <- summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.A")) %>%
  select(KEY.individ, know.table.type, obs_know_person)

base_df <- endline_data %>%
  left_join(summ_know_A_df, by = "KEY.individ") %>%
  mutate(
    observed_table_a = !is.na(know.table.type) & obs_know_person > 0,
    missing_know_table = !in_know_table,
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group
  )

observed_a_df <- base_df %>%
  filter(observed_table_a) %>%
  mutate(attrit_know_table = FALSE)

missing_df <- base_df %>%
  filter(missing_know_table) %>%
  mutate(attrit_know_table = TRUE)

stopifnot(nrow(observed_a_df) == 1141)
stopifnot(nrow(missing_df) == 252)

coef_names <- c(
  "assigned_treatment::ink",
  "assigned_treatment::calendar",
  "assigned_treatment::bracelet",
  "assigned_treatment::ink:assigned_dist_group::close",
  "assigned_treatment::ink:assigned_dist_group::far",
  "assigned_treatment::calendar:assigned_dist_group::close",
  "assigned_treatment::calendar:assigned_dist_group::far",
  "assigned_treatment::bracelet:assigned_dist_group::close",
  "assigned_treatment::bracelet:assigned_dist_group::far",
  "femaleTRUE",
  "age.census",
  "mu_d"
)

extract_terms <- function(model, sim, model_name) {
  coefs <- coef(model)
  vc <- vcov(model)

  map_dfr(intersect(names(coefs), coef_names), function(term) {
    tibble(
      sim = sim,
      model = model_name,
      term = term,
      estimate = unname(coefs[term]),
      se = sqrt(vc[term, term])
    )
  })
}

contrast_row <- function(model, sim, model_name, label, lhs, rhs) {
  coefs <- coef(model)
  vc <- vcov(model)
  contrast <- rep(0, length(coefs))
  names(contrast) <- names(coefs)
  contrast[lhs] <- 1
  contrast[rhs] <- -1

  tibble(
    sim = sim,
    model = model_name,
    term = label,
    estimate = sum(contrast * coefs),
    se = sqrt(drop(t(contrast) %*% vc %*% contrast))
  )
}

fit_one_sim <- function(sim) {
  imputed_missing_a <- missing_df %>%
    slice_sample(n = nrow(missing_df) / 2)

  sim_df <- bind_rows(observed_a_df, imputed_missing_a)

  pooled_no_covs <- feols(
    attrit_know_table ~ i(assigned_treatment, ref = "control") | county,
    data = sim_df,
    cluster = ~cluster.id
  )

  dist_no_covs <- feols(
    attrit_know_table ~ i(assigned_treatment, assigned_dist_group, ref = "control") | county,
    data = sim_df,
    cluster = ~cluster.id
  )

  pooled_covs <- feols(
    attrit_know_table ~ i(assigned_treatment, ref = "control") + .[l_cov_vars] + mu_d | county,
    data = sim_df,
    cluster = ~cluster.id
  )

  dist_covs <- feols(
    attrit_know_table ~ i(assigned_treatment, assigned_dist_group, ref = "control") +
      .[l_cov_vars] + mu_d | county,
    data = sim_df,
    cluster = ~cluster.id
  )

  term_rows <- bind_rows(
    extract_terms(pooled_no_covs, sim, "pooled_no_covs"),
    extract_terms(dist_no_covs, sim, "dist_no_covs"),
    extract_terms(pooled_covs, sim, "pooled_covs"),
    extract_terms(dist_covs, sim, "dist_covs"),
    contrast_row(
      pooled_no_covs, sim, "pooled_no_covs",
      "bracelet_minus_calendar",
      "assigned_treatment::bracelet", "assigned_treatment::calendar"
    ),
    contrast_row(
      pooled_covs, sim, "pooled_covs",
      "bracelet_minus_calendar",
      "assigned_treatment::bracelet", "assigned_treatment::calendar"
    ),
    contrast_row(
      dist_no_covs, sim, "dist_no_covs",
      "bracelet_close_minus_calendar_close",
      "assigned_treatment::bracelet:assigned_dist_group::close",
      "assigned_treatment::calendar:assigned_dist_group::close"
    ),
    contrast_row(
      dist_covs, sim, "dist_covs",
      "bracelet_close_minus_calendar_close",
      "assigned_treatment::bracelet:assigned_dist_group::close",
      "assigned_treatment::calendar:assigned_dist_group::close"
    ),
    contrast_row(
      dist_no_covs, sim, "dist_no_covs",
      "bracelet_far_minus_calendar_far",
      "assigned_treatment::bracelet:assigned_dist_group::far",
      "assigned_treatment::calendar:assigned_dist_group::far"
    ),
    contrast_row(
      dist_covs, sim, "dist_covs",
      "bracelet_far_minus_calendar_far",
      "assigned_treatment::bracelet:assigned_dist_group::far",
      "assigned_treatment::calendar:assigned_dist_group::far"
    )
  )

  control_rows <- sim_df %>%
    filter(assigned_treatment == "control") %>%
    summarise(
      control_mean_pooled = mean(attrit_know_table),
      control_mean_close = mean(attrit_know_table[assigned_dist_group == "close"]),
      control_mean_far = mean(attrit_know_table[assigned_dist_group == "far"]),
      n = n(),
      .groups = "drop"
    ) %>%
    pivot_longer(everything(), names_to = "term", values_to = "estimate") %>%
    mutate(
      sim = sim,
      model = case_when(
        term %in% c("control_mean_pooled", "n") ~ "pooled_no_covs",
        TRUE ~ "dist_no_covs"
      ),
      se = 0
    ) %>%
    select(sim, model, term, estimate, se)

  bind_rows(term_rows, control_rows)
}

sim_results_path <- file.path(params$output_path, "simulated-table-A-attrition-results.csv")

if (file.exists(sim_results_path)) {
  sim_results <- read_csv(sim_results_path, show_col_types = FALSE)
} else {
  sim_results <- map_dfr(seq_len(params$n_sims), fit_one_sim, .progress = TRUE)

  sim_results %>%
    write_csv(sim_results_path)
}

combined_results <- sim_results %>%
  group_by(model, term) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    within_var = mean(se^2, na.rm = TRUE),
    between_var = var(estimate, na.rm = TRUE),
    total_se = sqrt(within_var + (1 + 1 / n()) * between_var),
    n_sims = n(),
    pval = 2 * pnorm(-abs(mean_estimate / total_se)),
    .groups = "drop"
  ) %>%
  rename(estimate = mean_estimate)

combined_results %>%
  write_csv(file.path(params$output_path, "simulated-table-A-attrition-combined.csv"))

fmt <- function(x) sprintf("%.3f", x)
stars <- function(p) case_when(
  is.na(p) ~ "",
  p < 0.01 ~ "$^{***}$",
  p < 0.05 ~ "$^{**}$",
  p < 0.1 ~ "$^{*}$",
  TRUE ~ ""
)
pval_fmt <- function(p) case_when(
  is.na(p) ~ "",
  p < 0.001 ~ "$<$0.001",
  TRUE ~ sprintf("%.3f", p)
)

cell <- function(model, term) {
  row <- combined_results %>% filter(model == !!model, term == !!term)
  if (nrow(row) == 0) return("")
  paste0(fmt(row$estimate), stars(row$pval), "\\\\{(", fmt(row$total_se), ")}")
}

stat_cell <- function(model, term) {
  row <- combined_results %>% filter(model == !!model, term == !!term)
  if (nrow(row) == 0) return("")
  fmt(row$estimate)
}

pval_cell <- function(model, term) {
  row <- combined_results %>% filter(model == !!model, term == !!term)
  if (nrow(row) == 0) return("")
  pval_fmt(row$pval)
}

rows <- tribble(
  ~label, ~term, ~type,
  "Ink", "assigned_treatment::ink", "coef_pooled",
  "Calendar", "assigned_treatment::calendar", "coef_pooled",
  "Bracelet", "assigned_treatment::bracelet", "coef_pooled",
  "Ink $\\times$ Close", "assigned_treatment::ink:assigned_dist_group::close", "coef_dist",
  "Ink $\\times$ Far", "assigned_treatment::ink:assigned_dist_group::far", "coef_dist",
  "Calendar $\\times$ Close", "assigned_treatment::calendar:assigned_dist_group::close", "coef_dist",
  "Calendar $\\times$ Far", "assigned_treatment::calendar:assigned_dist_group::far", "coef_dist",
  "Bracelet $\\times$ Close", "assigned_treatment::bracelet:assigned_dist_group::close", "coef_dist",
  "Bracelet $\\times$ Far", "assigned_treatment::bracelet:assigned_dist_group::far", "coef_dist",
  "Female", "femaleTRUE", "coef_cov",
  "Age", "age.census", "coef_cov",
  "Expected distance to PoT", "mu_d", "coef_cov"
)

tex_coef_rows <- rows %>%
  rowwise() %>%
  mutate(
    c1 = case_when(type == "coef_pooled" ~ cell("pooled_no_covs", term), TRUE ~ ""),
    c2 = case_when(type == "coef_dist" ~ cell("dist_no_covs", term), TRUE ~ ""),
    c3 = case_when(type %in% c("coef_pooled", "coef_cov") ~ cell("pooled_covs", term), TRUE ~ ""),
    c4 = case_when(type %in% c("coef_dist", "coef_cov") ~ cell("dist_covs", term), TRUE ~ ""),
    row = paste0(label, " & \\makecell[c]{", c1, "} & \\makecell[c]{", c2,
                 "} & \\makecell[c]{", c3, "} & \\makecell[c]{", c4, "} \\\\")
  ) %>%
  ungroup() %>%
  pull(row)

tex_table <- c(
  "\\centering",
  "\\caption{Simulated Table A Knowledge-Table Attrition}",
  "\\label{tab:simulated-table-a-attrition}",
  "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{",
  "\\begin{tabular}[t]{lcccc}",
  "\\toprule",
  " & Pooled & Treatment $\\times$ Distance & Pooled & Treatment $\\times$ Distance \\\\",
  "Dependent variable: Missing from Table A & (1) & (2) & (3) & (4) \\\\",
  "\\midrule",
  tex_coef_rows,
  "\\midrule",
  paste0("Control mean & ", stat_cell("pooled_no_covs", "control_mean_pooled"), " &  & ",
         stat_cell("pooled_no_covs", "control_mean_pooled"), " &  \\\\"),
  paste0("Control mean, Close &  & ", stat_cell("dist_no_covs", "control_mean_close"), " &  & ",
         stat_cell("dist_no_covs", "control_mean_close"), " \\\\"),
  paste0("Control mean, Far &  & ", stat_cell("dist_no_covs", "control_mean_far"), " &  & ",
         stat_cell("dist_no_covs", "control_mean_far"), " \\\\"),
  "County FEs & Yes & Yes & Yes & Yes \\\\",
  "RF controls & No & No & Yes & Yes \\\\",
  paste0("$H_0$: Bracelet = Calendar, $p$-value & ",
         pval_cell("pooled_no_covs", "bracelet_minus_calendar"), " &  & ",
         pval_cell("pooled_covs", "bracelet_minus_calendar"), " &  \\\\"),
  paste0("$H_0$: Bracelet $\\times$ Close = Calendar $\\times$ Close, $p$-value &  & ",
         pval_cell("dist_no_covs", "bracelet_close_minus_calendar_close"), " &  & ",
         pval_cell("dist_covs", "bracelet_close_minus_calendar_close"), " \\\\"),
  paste0("$H_0$: Bracelet $\\times$ Far = Calendar $\\times$ Far, $p$-value &  & ",
         pval_cell("dist_no_covs", "bracelet_far_minus_calendar_far"), " &  & ",
         pval_cell("dist_covs", "bracelet_far_minus_calendar_far"), " \\\\"),
  "Observations & 1,267 & 1,267 & 1,267 & 1,267 \\\\",
  "Simulations & 1,000 & 1,000 & 1,000 & 1,000 \\\\",
  "\\bottomrule",
  "\\end{tabular}}",
  "\\\\[0.5em]",
  "\\footnotesize{Each simulation randomly assigns exactly 126 of the 252 missing knowledge-table respondents to Table A and treats them as missing Table A observations. Point estimates average across simulations. Standard errors combine clustered regression variance and simulation variance.}"
)

writeLines(tex_table, file.path(params$table_output_path, "simulated_table_A_attrition_regression_tbl.tex"))
