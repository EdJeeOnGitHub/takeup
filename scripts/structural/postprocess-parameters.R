#!/usr/bin/Rscript
script_options <- docopt::docopt(
  stringr::str_glue(
"Usage:
  scripts/structural/postprocess-parameters.R <fit-version> [options] [<chain>...]
  
Options:
  --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
  --output-path=<path>  Path to find results [default: temp-data]
  --model=<model>  Which model to postprocess
  --prior  Postprocess the prior predictive
  --save-error-draws  Save the entire posterior w/ each cluster's w^* draws

  "), 
  args = if (interactive()) "
  95
  --output-path=temp-data
  --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB
  1 
  " else commandArgs(trailingOnly = TRUE)
)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(tidybayes)

model_type = if_else(str_detect(script_options$model, "STRUCT"), "structural", "reduced form")
indiv_community_model = if_else(str_detect(script_options$model, "INDIV_DIST_COMMUNITY_FP"), TRUE, FALSE)
indiv_indiv_model = if_else(str_detect(script_options$model, "INDIV_DIST_INDIV_FP"), TRUE, FALSE)
indiv_model = indiv_community_model | indiv_indiv_model
source("R/structural/postprocess-functions.R")

# N.B. treat_idx (the second idx, is the mu (signalling) idx)
mu_idx_mapper = tibble(
  treat_idx = 1:4,
  mu_treatment = c("control", "ink", "calendar", "bracelet")
) %>%
  mutate(
    mu_treatment = factor(mu_treatment, levels = c("bracelet", "calendar", "ink", "control")) %>% fct_rev,
    mu_treatment = fct_relabel(mu_treatment, str_to_title)
    )
dist_idx_mapper = tibble(
  dist_treat_idx = 1:8,
  treatment = rep(c("control", "ink", "calendar", "bracelet"), 2),
  dist_group = rep(c("close", "far"), each = 4)
) %>%
  mutate(
    treatment = factor(treatment, levels = c("bracelet", "calendar", "ink", "control")) %>% fct_rev,
    treatment = fct_relabel(treatment, str_to_title)
    )
cluster_mapper = analysis_data %>%
  select(
    cluster_id,
    assigned_treatment,
    assigned_dist_group
  ) %>% unique()


fit_type_str = if_else(script_options$prior, "prior", "fit")
if (length(script_options$chain) > 1) {
  chain_str = str_glue("{min(script_options$chain)}-{max(script_options$chain)}")
} else {
  chain_str = script_options$chain
}


takeup_param_draws = load_param_draws(
  fit_version = script_options$fit_version,
  model = script_options$model,
  chain = script_options$chain,
  prior_predictive = script_options$prior,
  input_path = script_options$input_path,
  beta_intercept, beta_ink_effect, beta_calendar_effect, beta_bracelet_effect
)

takeup_param_draws = load_param_draws(
  fit_version = script_options$fit_version,
  model = script_options$model,
  chain = script_options$chain,
  prior_predictive = script_options$prior,
  input_path = script_options$input_path,
  # 1:4 are beta private value
  beta[i_beta], 
  #  this might index clusters but is homogeneous across clusters
  u_sd[i_u],
  # distance param
  dist_beta_v[i_dist_beta_v],
  base_mu_rep
)

clean_takeup_param_df = takeup_param_draws %>%
  filter(i_u == 1, i_dist_beta_v == 1)  %>%
  filter(i_beta <= 4) %>%
  select(-i_u, -i_dist_beta_v)  %>%
  pivot_wider(
    names_from = i_beta,
    values_from = beta,
    names_prefix = "beta_"
  ) %>%
  pivot_longer(
    cols = where(is_rvar),
  )  %>%
  mean_qi(value) %>%
  to_broom_names() %>%
  mutate(
    i_treat = str_extract(name, "\\d+") %>% as.numeric
  )

takeup_obs_draws = load_param_draws(
  fit_version = script_options$fit_version,
  model = script_options$model,
  chain = script_options$chain,
  prior_predictive = script_options$prior,
  input_path = script_options$input_path,
  centered_cluster_beta_beliefs[i_clust, i_treat],
  centered_cluster_dist_beta_beliefs[i_clust, i_treat]
)

clean_takeup_obs_df =  takeup_obs_draws %>%
  filter(i_clust == 1) %>%
  select(-i_clust) %>%
  pivot_longer(where(is_rvar)) %>%
  select(
    i_treat, name, value
  ) %>%
  mean_qi(value)  %>%
  to_broom_names() 

takeup_wtp_params = load_param_draws(
  fit_version = script_options$fit_version,
  model = script_options$model,
  chain = script_options$chain,
  prior_predictive = script_options$prior,
  input_path = script_options$input_path,
  hyper_wtp_mu,
  wtp_value_utility
)

clean_takeup_wtp_df = takeup_wtp_params %>%
  pivot_longer(where(is_rvar)) %>%
  mean_qi(value) %>%
  to_broom_names()


clean_param_df = bind_rows(
  clean_takeup_param_df %>% mutate(type = "param"),
  clean_takeup_obs_df %>% mutate(type = "obs"),
  clean_takeup_wtp_df %>% mutate(type = "wtp")
)

library(kableExtra)

summ_param_df = clean_param_df %>%
  select(
    name,
    i_treat,
    type, value, conf.low, conf.high
  ) %>%
  filter(name != "wtp_value_utility" ) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  bind_rows(
    clean_param_df %>%
      filter(name == "wtp_value_utility") %>%
      select(
        name,
        type, value, conf.low, conf.high
      ) %>%
      mutate(across(where(is.numeric), round, 5)) 
  ) %>%
  mutate(
    estim_value = linebreak(
      paste0(value, "\n", "(", conf.low, ", ", conf.high, ")"), 
      align = "c"
      )
  ) %>%
  mutate(
    treatment = case_when(
      i_treat == 1 ~ "Control",
      i_treat == 2 ~ "Ink",
      i_treat == 3 ~ "Calendar",
      i_treat == 4 ~ "Bracelet"
    ),
    clean_name = case_when(
      name == "u_sd" ~ "\\sigma_u",
      name == "dist_beta_v" ~ "\\delta",
      name == "base_mu_rep" ~ "\\lambda",
      str_detect(name, "beta_\\d") ~ str_glue("\\beta_{{{treatment}}}"),
      str_detect(name, "centered_cluster_beta_beliefs") ~ str_glue("\\beta^O_{{{treatment}}}"),
      str_detect(name, "centered_cluster_dist_beta_beliefs") ~ str_glue("\\gamma^O_{{{treatment}}}"),
      name == "hyper_wtp_mu" ~ "\\psi \\text{ (USD)}",
      name == "wtp_value_utility" ~ "\\rho"
    )
  ) %>%
  mutate(type = factor(type, levels = c("param", "obs", "wtp"))) %>%
  arrange(type, name, i_treat) %>%
  select(
    clean_name, estim_value, type
  )



param_tbl = summ_param_df %>%
  select(clean_name, estim_value) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    col.names = c("Parameter", "Estimate"),
    align = "lc"
  ) %>%
  pack_rows(
    "Takeup Parameters", 1, 7
  )  %>%
  pack_rows(
    "Visibility Parameters", 8, 15
  ) %>%
  pack_rows(
    "WTP Parameter", 16, 17
  )

custom_save_latex_table = function(table, table_filename){
  table_conn = file(
    table_filename
  )
  # Create space in latex due to rmd bug with \\ immediately followed by [
  # table = table  %>%
  #   str_replace_all(., "removeme12345", "  ")

  attr(table, "kable_meta")$contents = str_replace_all(attr(table, "kable_meta")$contents, "removeme12345", " ")
  table[1] = str_replace_all(table[1], "removeme12345", " ")


  clean_table = table %>%
    str_remove(
      ., 
      fixed("\\begin{table}")
    ) %>%
    str_remove(
      .,
      "\\\\caption\\{.*\\}"
    ) %>%
    str_remove(
      ., 
      "\\\\end\\{table\\}"
    ) 
    
    
    clean_table %>%
      writeLines(
        table_conn
      )
    close(table_conn)

    return(table)
}


custom_save_latex_table(
  param_tbl,
  str_glue("{script_options$output_path}/all-structural-parameters-table-{fit_type_str}-{chain_str}-{script_options$model}.tex")
)
