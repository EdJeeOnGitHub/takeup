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
source(file.path("scratch", "reduced-form-setup.R"))

dir.create(params$table_output_path, showWarnings = FALSE, recursive = TRUE)

summ_know_A_df <- summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.A")) %>%
  select(KEY.individ, know.table.type, obs_know_person)

first_responder_balance_df <- endline_data %>%
  left_join(summ_know_A_df, by = "KEY.individ") %>%
  mutate(
    observed_table_a = !is.na(know.table.type) & obs_know_person > 0,
    census_first_responder = !(KEY.individ %in% census_data$KEY.individ[census_data$endline.backup == TRUE]),
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group,
    female_num = as.numeric(female),
    has_phone = as.numeric(have_phone == "Yes")
  ) %>%
  filter(observed_table_a, census_first_responder)

balance_vars <- c(
  age.census = "Age",
  female_num = "Female",
  has_phone = "Has phone",
  mu_d = "Expected distance to PoT"
)

first_responder_balance_fits <- feols(
  as.formula(
    paste0(
      "c(", paste(names(balance_vars), collapse = ", "), ") ~ ",
      "i(assigned_treatment, assigned_dist_group, ref = 'control') | county"
    )
  ),
  data = first_responder_balance_df,
  cluster = ~cluster.id
)

setFixest_dict(c(
  "assigned_treatment::ink:assigned_dist_group::close" = "Ink $\\times$ Close",
  "assigned_treatment::ink:assigned_dist_group::far" = "Ink $\\times$ Far",
  "assigned_treatment::calendar:assigned_dist_group::close" = "Calendar $\\times$ Close",
  "assigned_treatment::calendar:assigned_dist_group::far" = "Calendar $\\times$ Far",
  "assigned_treatment::bracelet:assigned_dist_group::close" = "Bracelet $\\times$ Close",
  "assigned_treatment::bracelet:assigned_dist_group::far" = "Bracelet $\\times$ Far",
  age.census = "Age",
  female_num = "Female",
  has_phone = "Has phone",
  mu_d = "Expected distance to PoT"
))

tex_postprocessing <- function(tex) {
  tex %>%
    str_remove("\\\\begin\\{table\\}\\[htbp\\]") %>%
    str_remove("\\\\end\\{table\\}")
}

etable(
  first_responder_balance_fits,
  headers = unname(balance_vars),
  depvar = FALSE,
  fitstat = c("n", "r2"),
  se.below = TRUE,
  tex = TRUE,
  title = "Covariate Balance in the First-Responder FOB Sample",
  label = "tab:first-responder-balance-diagnostics",
  notes = paste(
    "Sample is the Table 5 first-responder FOB sample: observed Table A respondents with the census non-backup flag.",
    "Standard errors clustered by community."
  ),
  file = file.path(params$table_output_path, "first-responder-balance-diagnostics.tex"),
  replace = TRUE,
  postprocess.tex = tex_postprocessing,
  digits = 3,
  digits.stats = 3,
  drop.section = "fixef",
  style.df = style.df(depvar.title = "", fixef.title = "", var.title = "", stats.title = "")
)

first_responder_balance_df %>%
  summarise(
    n = n(),
    across(all_of(names(balance_vars)), ~mean(.x, na.rm = TRUE))
  ) %>%
  write_csv(file.path(params$output_path, "first-responder-balance-diagnostics-means.csv"))
