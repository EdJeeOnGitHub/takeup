library(tidyverse)
library(broom)
library(data.table)
library(kableExtra)
library(knitr)
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
source(file.path("R/common/analysis.R"))
source(file.path("R/structural/legacy-utils.R"))
source(file.path("multilvlr", "multilvlr_util.R"))
source(file.path("scratch", "reduced-form-setup.R"))

dir.create("temp-data/tidy-rf-tes", showWarnings = FALSE, recursive = TRUE)
dir.create(params$table_output_path, showWarnings = FALSE, recursive = TRUE)

summ_know_A_df <- summ_endline_know_table %>%
  filter(fct_match(know.table.type, "table.A")) %>%
  select(KEY.individ, know.table.type, obs_know_person)

conservative_table_a_attrition_df <- endline_data %>%
  left_join(summ_know_A_df, by = "KEY.individ") %>%
  mutate(
    observed_table_a_main = !is.na(know.table.type) & obs_know_person > 0,
    missing_know_table = !in_know_table,
    conservative_table_a_frame = observed_table_a_main | missing_know_table,
    attrit_know_table = missing_know_table,
    assigned_treatment = assigned.treatment,
    assigned_dist_group = dist.pot.group,
    wt = 1
  ) %>%
  filter(conservative_table_a_frame)

conservative_table_a_attrition_counts <- conservative_table_a_attrition_df %>%
  count(
    assigned_treatment,
    assigned_dist_group,
    attrit_know_table,
    name = "n"
  ) %>%
  group_by(assigned_treatment, assigned_dist_group) %>%
  mutate(cell_n = sum(n), share = n / cell_n) %>%
  ungroup()

conservative_table_a_attrition_counts %>%
  write_csv("temp-data/tidy-rf-tes/conservative-table-A-attrition-counts.csv")

conservative_table_a_attrition_no_covs = function(data, weights) {
  feols(
    attrit_know_table ~ 0 + assigned_treatment + assigned_dist_group +
      i(assigned_treatment, assigned_dist_group, "control") | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

conservative_table_a_attrition_covs = function(data, weights) {
  feols(
    attrit_know_table ~ 0 + assigned_treatment + assigned_dist_group +
      i(assigned_treatment, assigned_dist_group, "control") +
      .[l_cov_vars] + mu_d | county,
    data = data,
    nthreads = 1,
    weights = ~wt
  )
}

conservative_table_a_attrition_no_covs_output <- wrapper_function(
  data = conservative_table_a_attrition_df,
  regression_spec = conservative_table_a_attrition_no_covs,
  tidy_summ_path = "temp-data/tidy-rf-tes/conservative-table-A-attrition-tidy-tes.csv",
  table_name = "conservative_table_A_attrition_spec_tbl",
  table_options = list(
    caption = "Conservative Table A Knowledge-Table Attrition",
    dependent_var = "Dependent variable: Not in Table A"
  )
)

conservative_table_a_attrition_covs_output <- wrapper_function(
  data = conservative_table_a_attrition_df,
  regression_spec = conservative_table_a_attrition_covs,
  tidy_summ_path = "temp-data/tidy-rf-tes/conservative-table-A-attrition-covs-tidy-tes.csv",
  table_name = "conservative_table_A_attrition_covs_spec_tbl",
  table_options = list(
    caption = "Conservative Table A Knowledge-Table Attrition",
    dependent_var = "Dependent variable: Not in Table A"
  )
)

conservative_table_a_attrition_no_covs_output$tidy_summary %>%
  select(assigned_treatment, assigned_dist_group, estimate, std_error, pval) %>%
  print(n = Inf)

conservative_table_a_attrition_covs_output$tidy_summary %>%
  select(assigned_treatment, assigned_dist_group, estimate, std_error, pval) %>%
  print(n = Inf)
