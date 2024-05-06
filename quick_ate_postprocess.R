#!/usr/bin/Rscript
script_options <- docopt::docopt(
  stringr::str_glue(
"Usage:
  quick_ate_postprocess.R <fit-version> [options] [<chain>...]
  
Options:
  --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
  --output-path=<path>  Path to find results [default: temp-data]
  --model=<model>  Which model to postprocess
  --prior  Postprocess the prior predictive
  --save-error-draws  Save the entire posterior w/ each cluster's w^* draws

  "), 
  args = if (interactive()) "
  104
  --output-path=temp-data
  --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_INDIV_DIST_INDIV_FP
  --save-error-draws
  1 2 
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
source("quick_postprocess_functions.R")
create_tes = function(.data, group_var, levels = FALSE) {
 if (levels == TRUE) {
  levels_mult = 0
 }  else {
  levels_mult = 1
 }

  tes =  bind_rows(
    .data %>%
      group_by(fit_type, model, fit_version, {{ group_var }}, dist_group) %>%
      summarise(
        pr_takeup = rvar_mean(pr_takeup),
        .groups = "drop"
      )   %>%
      group_by(dist_group) %>%
      mutate(
        pr_takeup = if_else({{ group_var }} == "Control", pr_takeup, pr_takeup - pr_takeup[{{ group_var }} == "Control"]*levels_mult)
      ),
    .data %>%
      group_by(fit_type, model, fit_version, {{ group_var }}) %>%
      summarise(
        pr_takeup = rvar_mean(pr_takeup),
        .groups = "drop"
      ) %>%
      mutate(
        pr_takeup = if_else({{ group_var }} == "Control", pr_takeup, pr_takeup - pr_takeup[{{ group_var }} == "Control"]*levels_mult)
      )  %>%
      mutate(dist_group = "combined")
  ) %>%
    mutate(
      dist_group = factor(dist_group, levels = c("far", "close", "combined"))
    ) 
  return(tes)
}


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


if (indiv_indiv_model) {

  cluster_error_draws_raw = load_param_draws(
    fit_version = script_options$fit_version,
    model = script_options$model,
    chain = script_options$chain,
    prior_predictive = script_options$prior,
    input_path = script_options$input_path,
    structural_cluster_takeup_prob[i]
  ) %>%
  bind_cols(
    analysis_data %>%
      mutate(indiv_id = row_number()) %>%
      select(
        cluster_id,
        assigned_treatment,
        assigned_dist_group
      ) 
  )

  cluster_error_draws = cluster_error_draws_raw %>%
    mutate(
      treatment = str_to_title(assigned_treatment), 
      dist_group = assigned_dist_group,
      mu_treatment = treatment
      ) %>%
    rename(pr_takeup = structural_cluster_takeup_prob) 


## Since output from indiv model much larger (9,800 rvars, instead of 144)
## calculate TEs here in same format as tables and save to RDS.

init_incentive_tes = cluster_error_draws %>%
      filter(dist_group == assigned_dist_group) %>%
      filter((mu_treatment == treatment) | model_type == "reduced form") %>%
      create_tes(group_var = treatment) %>%
      mutate(estimand = "overall") 

incentive_tes = bind_rows(
  init_incentive_tes,
  init_incentive_tes %>%
    filter(treatment %in% c("Bracelet", "Calendar")) %>%
    group_by(
      fit_type,
      model,
      fit_version,
      estimand,
      dist_group
    ) %>%
    summarise(
      pr_takeup = pr_takeup[treatment == "Bracelet"] - pr_takeup[treatment == "Calendar"]
    ) %>%
    mutate(
      treatment = "bracelet_minus_calendar"
    )
) %>%
  ungroup()


wide_incentive_tes = incentive_tes %>%
  pivot_wider(
    names_from = dist_group,
    values_from = pr_takeup
  )  %>%
  mutate(
    far_minus_close = far - close
  )  
  
wide_te_rvar_cols = wide_incentive_tes %>%
  select(where(is_rvar)) %>%
  colnames()
  
  
wide_incentive_tes  %>%
    pivot_longer(
      all_of(wide_te_rvar_cols),
      names_to = "variable"
    )   %>%
    saveRDS(
      file.path(
        script_options$output_path,
        str_glue(
          "rvar_processed_dist_{fit_type_str}{script_options$fit_version}_INDIV_incentive_tes_{script_options$model}_{chain_str}.rds"
        )
      ) 
    )

}


if (indiv_community_model) {

  cluster_error_draws_raw = load_param_draws(
    fit_version = script_options$fit_version,
    model = script_options$model,
    chain = script_options$chain,
    prior_predictive = script_options$prior,
    input_path = script_options$input_path,
    takeup_pr[i]
  ) %>%
  bind_cols(
    analysis_data %>%
      mutate(indiv_id = row_number()) %>%
      select(
        cluster_id,
        assigned_treatment,
        assigned_dist_group
      ) 
  )

  cluster_error_draws = cluster_error_draws_raw %>%
    mutate(
      treatment = str_to_title(assigned_treatment), 
      dist_group = assigned_dist_group,
      mu_treatment = treatment
      ) %>%
    rename(pr_takeup = takeup_pr) 


## Since output from indiv model much larger (9,800 rvars, instead of 144)
## calculate TEs here in same format as tables and save to RDS.

init_incentive_tes = cluster_error_draws %>%
      filter(dist_group == assigned_dist_group) %>%
      filter((mu_treatment == treatment) | model_type == "reduced form") %>%
      create_tes(group_var = treatment) %>%
      mutate(estimand = "overall") 

incentive_tes = bind_rows(
  init_incentive_tes,
  init_incentive_tes %>%
    filter(treatment %in% c("Bracelet", "Calendar")) %>%
    group_by(
      fit_type,
      model,
      fit_version,
      estimand,
      dist_group
    ) %>%
    summarise(
      pr_takeup = pr_takeup[treatment == "Bracelet"] - pr_takeup[treatment == "Calendar"]
    ) %>%
    mutate(
      treatment = "bracelet_minus_calendar"
    )
) %>%
  ungroup()


wide_incentive_tes = incentive_tes %>%
  pivot_wider(
    names_from = dist_group,
    values_from = pr_takeup
  )  %>%
  mutate(
    far_minus_close = far - close
  )  
  
wide_te_rvar_cols = wide_incentive_tes %>%
  select(where(is_rvar)) %>%
  colnames()
  
  
wide_incentive_tes  %>%
    pivot_longer(
      all_of(wide_te_rvar_cols),
      names_to = "variable"
    )   %>%
    saveRDS(
      file.path(
        script_options$output_path,
        str_glue(
          "rvar_processed_dist_{fit_type_str}{script_options$fit_version}_INDIV_incentive_tes_{script_options$model}_{chain_str}.rds"
        )
      ) 
    )


}

if (!indiv_model) {
  if (model_type == "structural") {
    cluster_error_draws_raw = load_param_draws(
      fit_version = script_options$fit_version,
      model = script_options$model,
      chain = script_options$chain,
      prior_predictive = script_options$prior,
      input_path = script_options$input_path,
      cluster_cf_cutoff[dist_treat_idx, treat_idx, cluster_idx],
      total_error_sd[treat_idx]
    )
  } else {
    cluster_error_draws_raw = load_param_draws(
      fit_version = script_options$fit_version,
      model = script_options$model,
      chain = script_options$chain,
      prior_predictive = script_options$prior,
      input_path = script_options$input_path,
      cluster_cf_cutoff[dist_treat_idx, treat_idx, cluster_idx]
    ) %>%
      mutate(total_error_sd = 1)
  }


rvar_pnorm = rfun(pnorm)
cluster_error_draws = cluster_error_draws_raw %>%
  mutate(
    pr_takeup = rvar_pnorm(-cluster_cf_cutoff, sd = total_error_sd)
  ) %>%
  left_join(
    dist_idx_mapper,
    by = c("dist_treat_idx")
  ) %>%
  left_join(
    mu_idx_mapper,
    by = "treat_idx"
 ) %>%
  left_join(
    cluster_mapper,
    by = c("cluster_idx" = "cluster_id")
  ) %>%
  mutate(model_type = model_type)

}



rvar_cols = cluster_error_draws %>%
  select(where(is_rvar)) %>%
  colnames()


if (script_options$save_error_draws) {
  cluster_error_draws %>%
    pivot_longer(
      all_of(rvar_cols),
      names_to = "variable"
    )  %>%
    saveRDS(
      file.path(
        script_options$output_path,
        str_glue(
          "rvar_processed_dist_{fit_type_str}{script_options$fit_version}_cluster_error_{script_options$model}_{chain_str}.rds"
        )
      ) 
    )
}




incentive_tes = cluster_error_draws %>%
      filter(dist_group == assigned_dist_group) %>%
      filter((mu_treatment == treatment) | model_type == "reduced form") %>%
      create_tes(group_var = treatment) %>%
      mutate(estimand = "overall")
if (model_type == "structural" & !indiv_model) {
  signal_tes = cluster_error_draws %>%
    filter(treatment == "Control") %>%
    create_tes(group_var = mu_treatment) %>%
    mutate(estimand = "signal")
  private_tes = cluster_error_draws %>%
    filter(dist_group == assigned_dist_group) %>%
    filter(mu_treatment == "Control") %>%
    create_tes(group_var = treatment) %>%
    filter(dist_group == "combined") %>%
    mutate(estimand = "private")
} else {
  signal_tes = NULL
  private_tes = NULL
}

all_tes = bind_rows(
  incentive_tes,
  signal_tes,
  private_tes
) %>%
  ungroup()

all_rvar_cols = all_tes %>%
  select(where(is_rvar)) %>%
  colnames()



all_tes %>%
  pivot_longer(all_of(all_rvar_cols), names_to = "variable") %>%
  saveRDS(
    file.path(
      script_options$output_path,
      str_glue(
        "rvar_processed_dist_{fit_type_str}{script_options$fit_version}_ates_{script_options$model}_{chain_str}.rds"
      )
    ) 
  )

incentive_levels = cluster_error_draws %>%
      filter(dist_group == assigned_dist_group) %>%
      filter((mu_treatment == treatment) | model_type == "reduced form") %>%
      create_tes(group_var = treatment, levels = TRUE) %>%
      mutate(estimand = "overall")


incentive_levels %>%
  pivot_longer(where(is_rvar), names_to = "variable") %>%
  saveRDS(
    file.path(
      script_options$output_path,
      str_glue(
        "rvar_processed_dist_{fit_type_str}{script_options$fit_version}_levels_{script_options$model}_{chain_str}.rds"
      )
    ) 
  )

  

## Munging stuff

# all_tes = read_rds(
#     file.path(
#       script_options$output_path,
#       "rvar_processed_dist_fit95_ates_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_STRATA_FOB_1-2.rds"
#     ) 
# )
create_cis = function(.data) {
  med_fun = function(x) {
    mean_x = mean(x) %>% round(3)
    conf.low = quantile(x, 0.025) %>% round(3)
    conf.high = quantile(x, 0.975) %>% round(3)
    return(paste0(
      mean_x, " (", conf.low, ", ", conf.high, ")"
    ))
  }
  .data %>%
    mutate(across(where(is_rvar), med_fun))
}

# wide_struct_tes = all_tes %>% 
#     filter(estimand == "overall") %>%
#     select(
#       dist_treatment,
#       dist_group,
#       pr_takeup
#     ) %>%
#     pivot_wider(
#       names_from = dist_group,
#       values_from = pr_takeup
#     ) %>%
#     select(dist_treatment, combined, close, far) %>%
#     arrange(dist_treatment) %>%
#     bind_rows(
#       # bracelet minus calendar row
#       all_tes %>%
#         filter(dist_treatment %in% c("bracelet", "calendar")) %>%
#         filter(estimand == "overall") %>%
#         pivot_wider(names_from = dist_treatment, values_from = pr_takeup) %>%
#         mutate(
#           bracelet_minus_calendar = bracelet - calendar
#         ) %>%
#         select(dist_group, bracelet_minus_calendar) %>%
#         pivot_wider(
#           names_from = dist_group,
#           values_from = bracelet_minus_calendar
#         ) %>%
#         mutate(dist_treatment = "bracelet_minus_calendar")
#     )

# clean_wide_tbl = wide_struct_tes %>%
#   mutate(dist_treatment = factor(dist_treatment, levels = c(
#     "bracelet",
#     "calendar",
#     "ink",
#     "bracelet_minus_calendar",
#     "control"
#   ))) %>%
#   mutate(far_minus_close = far - close) %>%
#   arrange(dist_treatment)   %>%
#   create_cis()





# all_tes %>%
#     filter(estimand == "overall") %>%
#     select(
#       treatment,
#       dist_group,
#       pr_takeup
#     ) %>%
#     pivot_wider(
#       names_from = dist_group,
#       values_from = pr_takeup
#     ) %>%
#     select(treatment, combined, close, far) %>%
#     arrange(treatment) %>%
#     bind_rows(
#       # bracelet minus calendar row
#       all_tes %>%
#         filter(treatment %in% c("Bracelet", "Calendar")) %>%
#         filter(estimand == "overall") %>%
#         pivot_wider(names_from = treatment, values_from = pr_takeup) %>%
#         mutate(
#           bracelet_minus_calendar = Bracelet - Calendar
#         ) %>%
#         select(dist_group, bracelet_minus_calendar) %>%
#         pivot_wider(
#           names_from = dist_group,
#           values_from = bracelet_minus_calendar
#         ) %>%
#         mutate(treatment = "bracelet_minus_calendar")
#     ) %>%
#     create_cis() %>%
#     mutate(
#       treatment = factor(treatment, levels = c(
#         "bracelet_minus_calendar",
#         "Control",
#         "Ink",
#         "Calendar",
#         "Bracelet"
#       )) %>% fct_rev
#     ) %>%
#     arrange(treatment) %>%
#     View()
