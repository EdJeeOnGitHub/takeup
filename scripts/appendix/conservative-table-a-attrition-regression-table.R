library(tidyverse)
library(fixest)

params <- lst(
  table_output_path = "presentations/rf-tables/main-specs",
  output_path = "temp-data",
  show_probs = FALSE,
  width = 0.95,
  cache = FALSE,
  fit = FALSE,
  fit_version = "",
  stat = "std.error"
)

source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source("R/common/analysis.R")
source("R/structural/legacy-utils.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source(file.path("R", "reduced-form", "context.R"))
analysis_context <- takeup_get_analysis_context()
takeup_context_into_environment(analysis_context, environment())

dir.create(params$table_output_path, showWarnings = FALSE, recursive = TRUE)

summ_know_A_df <- summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.A")) %>%
  select(KEY.individ, know.table.type, obs_know_person)

full_sample_attrition_df <- all_endline_data %>%
  mutate(
    attrit_know_table = !in_know_table,
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group
  )

conservative_table_a_attrition_df <- endline_data %>%
  left_join(summ_know_A_df, by = "KEY.individ") %>%
  mutate(
    observed_table_a_main = !is.na(know.table.type) & obs_know_person > 0,
    missing_know_table = !in_know_table,
    conservative_table_a_frame = observed_table_a_main | missing_know_table,
    attrit_know_table = missing_know_table,
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group
  ) %>%
  filter(conservative_table_a_frame)

pval_fmt <- function(x) {
  case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "$<$0.001",
    TRUE ~ sprintf("%.3f", x)
  )
}

contrast_p <- function(model, lhs, rhs) {
  coefs <- coef(model)
  vc <- vcov(model)
  contrast <- rep(0, length(coefs))
  names(contrast) <- names(coefs)
  contrast[lhs] <- 1
  contrast[rhs] <- -1

  estimate <- sum(contrast * coefs)
  se <- sqrt(drop(t(contrast) %*% vc %*% contrast))
  pval <- 2 * pnorm(-abs(estimate / se))
  pval_fmt(pval)
}

wald_p <- function(model, terms) {
  terms <- intersect(terms, names(coef(model)))
  if (length(terms) == 0) return("")

  out <- wald(model, keep = paste(terms, collapse = "|"))
  pval_fmt(out$p)
}

mean_fmt <- function(x) sprintf("%.3f", x)

fit_full_pooled_no_covs <- feols(
  attrit_know_table ~ i(assigned_treatment, ref = "control") | county,
  data = full_sample_attrition_df,
  cluster = ~cluster.id
)

fit_full_dist_no_covs <- feols(
  attrit_know_table ~ assigned_dist_group + i(assigned_treatment, assigned_dist_group, ref = "control") | county,
  data = full_sample_attrition_df,
  cluster = ~cluster.id
)

fit_pooled_no_covs <- feols(
  attrit_know_table ~ i(assigned_treatment, ref = "control") | county,
  data = conservative_table_a_attrition_df,
  cluster = ~cluster.id
)

fit_dist_no_covs <- feols(
  attrit_know_table ~ assigned_dist_group + i(assigned_treatment, assigned_dist_group, ref = "control") | county,
  data = conservative_table_a_attrition_df,
  cluster = ~cluster.id
)

fit_pooled_covs <- feols(
  attrit_know_table ~ i(assigned_treatment, ref = "control") + .[l_cov_vars] + mu_d | county,
  data = conservative_table_a_attrition_df,
  cluster = ~cluster.id
)

fit_dist_covs <- feols(
  attrit_know_table ~ assigned_dist_group + i(assigned_treatment, assigned_dist_group, ref = "control") +
    .[l_cov_vars] + mu_d | county,
  data = conservative_table_a_attrition_df,
  cluster = ~cluster.id
)

control_mean_full_pooled <- full_sample_attrition_df %>%
  filter(assigned_treatment == "control") %>%
  summarise(mean = mean(attrit_know_table)) %>%
  pull(mean) %>%
  mean_fmt()

control_mean_full_close <- full_sample_attrition_df %>%
  filter(assigned_treatment == "control", assigned_dist_group == "close") %>%
  summarise(mean = mean(attrit_know_table)) %>%
  pull(mean) %>%
  mean_fmt()

control_mean_full_far <- full_sample_attrition_df %>%
  filter(assigned_treatment == "control", assigned_dist_group == "far") %>%
  summarise(mean = mean(attrit_know_table)) %>%
  pull(mean) %>%
  mean_fmt()

control_mean_pooled <- conservative_table_a_attrition_df %>%
  filter(assigned_treatment == "control") %>%
  summarise(mean = mean(attrit_know_table)) %>%
  pull(mean) %>%
  mean_fmt()

control_mean_close <- conservative_table_a_attrition_df %>%
  filter(assigned_treatment == "control", assigned_dist_group == "close") %>%
  summarise(mean = mean(attrit_know_table)) %>%
  pull(mean) %>%
  mean_fmt()

control_mean_far <- conservative_table_a_attrition_df %>%
  filter(assigned_treatment == "control", assigned_dist_group == "far") %>%
  summarise(mean = mean(attrit_know_table)) %>%
  pull(mean) %>%
  mean_fmt()

setFixest_dict(c(
  attrit_know_table = "Missing knowledge table",
  "assigned_treatment::ink" = "Ink",
  "assigned_treatment::calendar" = "Calendar",
  "assigned_treatment::bracelet" = "Bracelet",
  "assigned_treatment::ink:assigned_dist_group::close" = "Ink $\\times$ Close",
  "assigned_treatment::ink:assigned_dist_group::far" = "Ink $\\times$ Far",
  "assigned_treatment::calendar:assigned_dist_group::close" = "Calendar $\\times$ Close",
  "assigned_treatment::calendar:assigned_dist_group::far" = "Calendar $\\times$ Far",
  "assigned_treatment::bracelet:assigned_dist_group::close" = "Bracelet $\\times$ Close",
  "assigned_treatment::bracelet:assigned_dist_group::far" = "Bracelet $\\times$ Far",
  assigned_dist_groupfar = "Far",
  femaleTRUE = "Female",
  age.census = "Age",
  mu_d = "Expected distance to PoT"
))

pooled_terms <- c(
  "assigned_treatment::ink",
  "assigned_treatment::calendar",
  "assigned_treatment::bracelet"
)

close_terms <- c(
  "assigned_treatment::ink:assigned_dist_group::close",
  "assigned_treatment::calendar:assigned_dist_group::close",
  "assigned_treatment::bracelet:assigned_dist_group::close"
)

far_terms <- c(
  "assigned_treatment::ink:assigned_dist_group::far",
  "assigned_treatment::calendar:assigned_dist_group::far",
  "assigned_treatment::bracelet:assigned_dist_group::far"
)

tex_postprocessing <- function(tex) {
  tex %>%
    str_remove("\\\\begin\\{table\\}\\[htbp\\]") %>%
    str_remove("\\\\end\\{table\\}")
}

etable(
  fit_full_pooled_no_covs,
  fit_full_dist_no_covs,
  fit_pooled_no_covs,
  fit_dist_no_covs,
  fit_pooled_covs,
  fit_dist_covs,
  headers = c(
    "Full: Pooled",
    "Full: Treat. $\\times$ Dist.",
    "Cons.: Pooled",
    "Cons.: Treat. $\\times$ Dist.",
    "Cons.: Pooled",
    "Cons.: Treat. $\\times$ Dist."
  ),
  depvar = FALSE,
  fitstat = c("n", "r2"),
  se.below = TRUE,
  tex = TRUE,
  title = "Knowledge-Table Attrition",
  label = "tab:conservative-table-a-attrition-regression",
  notes = "Columns 1--2 use the full cleaned endline sample and define attrition as missing from any knowledge table. Columns 3--6 restrict to the conservative Table A frame, treating all 252 SMS-control respondents missing from any knowledge table as missing from Table A. Standard errors clustered by community.",
  file = file.path(params$table_output_path, "conservative_table_A_attrition_regression_tbl.tex"),
  replace = TRUE,
  postprocess.tex = tex_postprocessing,
  digits = 3,
  digits.stats = 3,
  drop.section = "fixef",
  style.df = style.df(depvar.title = "", fixef.title = "", var.title = "", stats.title = ""),
  extralines = list(
    "Control mean" = c(control_mean_full_pooled, "", control_mean_pooled, "", control_mean_pooled, ""),
    "Control mean, Close" = c("", control_mean_full_close, "", control_mean_close, "", control_mean_close),
    "Control mean, Far" = c("", control_mean_full_far, "", control_mean_far, "", control_mean_far),
    "County FEs" = rep("Yes", 6),
    "RF controls" = c("No", "No", "No", "No", "Yes", "Yes"),
    "$H_0$: Ink = Calendar = Bracelet = Control, $p$-value" = c(
      wald_p(fit_full_pooled_no_covs, pooled_terms),
      "",
      wald_p(fit_pooled_no_covs, pooled_terms),
      "",
      wald_p(fit_pooled_covs, pooled_terms),
      ""
    ),
    "$H_0$: All Close treatment cells = Control, $p$-value" = c(
      "",
      wald_p(fit_full_dist_no_covs, close_terms),
      "",
      wald_p(fit_dist_no_covs, close_terms),
      "",
      wald_p(fit_dist_covs, close_terms)
    ),
    "$H_0$: All Far treatment cells = Control, $p$-value" = c(
      "",
      wald_p(fit_full_dist_no_covs, far_terms),
      "",
      wald_p(fit_dist_no_covs, far_terms),
      "",
      wald_p(fit_dist_covs, far_terms)
    )
  )
)
