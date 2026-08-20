library(tidyverse)
library(fixest)

params <- lst(
  table_output_path = "presentations/tables",
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

selection_df <- endline_data %>%
  left_join(summ_know_A_df, by = "KEY.individ") %>%
  mutate(
    observed_table_a = !is.na(know.table.type) & obs_know_person > 0,
    missing_know_table = !in_know_table,
    conservative_table_a_frame = observed_table_a | missing_know_table,
    realized_first_responder = endline.type == "endline",
    census_first_responder = !(KEY.individ %in% census_data$KEY.individ[census_data$endline.backup == TRUE]),
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
  pval_fmt(2 * pnorm(-abs(estimate / se)))
}

mean_fmt <- function(x) sprintf("%.3f", x)

control_mean <- function(var, dist_group) {
  selection_df %>%
    filter(assigned_treatment == "control", assigned_dist_group == dist_group) %>%
    summarise(mean = mean(.data[[var]], na.rm = TRUE)) %>%
    pull(mean) %>%
    mean_fmt()
}

selection_spec <- function(outcome) {
  feols(
    as.formula(paste0(outcome, " ~ i(assigned_treatment, assigned_dist_group, ref = 'control') | county")),
    data = selection_df,
    cluster = ~cluster.id
  )
}

fit_table_a_observed <- selection_spec("observed_table_a")
fit_realized_first_responder <- selection_spec("realized_first_responder")
fit_census_first_responder <- selection_spec("census_first_responder")

setFixest_dict(c(
  observed_table_a = "Table A observed",
  realized_first_responder = "Realized first responder",
  census_first_responder = "Census first responder",
  "assigned_treatment::ink:assigned_dist_group::close" = "Ink $\\times$ Close",
  "assigned_treatment::ink:assigned_dist_group::far" = "Ink $\\times$ Far",
  "assigned_treatment::calendar:assigned_dist_group::close" = "Calendar $\\times$ Close",
  "assigned_treatment::calendar:assigned_dist_group::far" = "Calendar $\\times$ Far",
  "assigned_treatment::bracelet:assigned_dist_group::close" = "Bracelet $\\times$ Close",
  "assigned_treatment::bracelet:assigned_dist_group::far" = "Bracelet $\\times$ Far"
))

tex_postprocessing <- function(tex) {
  tex %>%
    str_remove("\\\\begin\\{table\\}\\[htbp\\]") %>%
    str_remove("\\\\end\\{table\\}")
}

etable(
  fit_table_a_observed,
  fit_realized_first_responder,
  fit_census_first_responder,
  headers = c("Table A observed", "Realized first responder", "Census first responder"),
  depvar = FALSE,
  fitstat = c("n", "r2"),
  se.below = TRUE,
  tex = TRUE,
  title = "Selection Margins for Table A Observability",
  label = "tab:first-responder-selection-diagnostics",
  notes = paste(
    "Sample treats all 252 SMS-control respondents missing from any knowledge table as missing from Table A.",
    "Realized first responder is defined from endline.type == endline.",
    "Census first responder is the non-backup flag used in the first-responder-only FOB table.",
    "Standard errors clustered by community."
  ),
  file = file.path(params$table_output_path, "first-responder-selection-diagnostics.tex"),
  replace = TRUE,
  postprocess.tex = tex_postprocessing,
  digits = 3,
  digits.stats = 3,
  drop.section = "fixef",
  style.df = style.df(depvar.title = "", fixef.title = "", var.title = "", stats.title = ""),
  extralines = list(
    "Control mean, Close" = c(
      control_mean("observed_table_a", "close"),
      control_mean("realized_first_responder", "close"),
      control_mean("census_first_responder", "close")
    ),
    "Control mean, Far" = c(
      control_mean("observed_table_a", "far"),
      control_mean("realized_first_responder", "far"),
      control_mean("census_first_responder", "far")
    ),
    "$H_0$: Bracelet $\\times$ Close = Calendar $\\times$ Close, $p$-value" = c(
      contrast_p(fit_table_a_observed, "assigned_treatment::bracelet:assigned_dist_group::close", "assigned_treatment::calendar:assigned_dist_group::close"),
      contrast_p(fit_realized_first_responder, "assigned_treatment::bracelet:assigned_dist_group::close", "assigned_treatment::calendar:assigned_dist_group::close"),
      contrast_p(fit_census_first_responder, "assigned_treatment::bracelet:assigned_dist_group::close", "assigned_treatment::calendar:assigned_dist_group::close")
    ),
    "$H_0$: Bracelet $\\times$ Far = Calendar $\\times$ Far, $p$-value" = c(
      contrast_p(fit_table_a_observed, "assigned_treatment::bracelet:assigned_dist_group::far", "assigned_treatment::calendar:assigned_dist_group::far"),
      contrast_p(fit_realized_first_responder, "assigned_treatment::bracelet:assigned_dist_group::far", "assigned_treatment::calendar:assigned_dist_group::far"),
      contrast_p(fit_census_first_responder, "assigned_treatment::bracelet:assigned_dist_group::far", "assigned_treatment::calendar:assigned_dist_group::far")
    )
  )
)

selection_df %>%
  group_by(assigned_treatment, assigned_dist_group) %>%
  summarise(
    n = n(),
    table_a_observed = mean(observed_table_a),
    realized_first_responder = mean(realized_first_responder),
    census_first_responder = mean(census_first_responder),
    .groups = "drop"
  ) %>%
  write_csv(file.path(params$output_path, "first-responder-selection-diagnostics-cell-means.csv"))
