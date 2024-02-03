#!/usr/bin/Rscript
script_options <- docopt::docopt(
  stringr::str_glue(
"Usage:
  equilibrium-vstar-comparison.R <struct-fit-version> <rf-fit-version> [options] [<chain>...]
  
Options:
  --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
  --output-path=<path>  Path to find results [default: temp-data]
  --struct-model=<struct-model>  Which structural model to postprocess
  --rf-model=<rf-model>  Which reduced form model to postprocess
  --prior  Postprocess the prior predictive
  
  "), 
  args = if (interactive()) "
  101 98
  --output-path=temp-data
  --struct-model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIER_FOB
  --rf-model=brm_flat_dist_fit
  1
  " else commandArgs(trailingOnly = TRUE)
)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(tidybayes)


source("quick_postprocess_functions.R")

# N.B. treat_idx (the second idx, is the mu (signalling) idx)
treat_idx_mapper = tibble(
  treat_idx = 1:4,
  treatment = c("control", "ink", "calendar", "bracelet")
) %>%
  mutate(
    treatment = factor(treatment, levels = c("bracelet", "calendar", "ink", "control")) %>% fct_rev,
    treatment = fct_relabel(treatment, str_to_title)

  )
dist_idx_mapper = tibble(
  dist_treat_idx = 1:8,
  dist_treatment = rep(c("control", "ink", "calendar", "bracelet"), 2),
  dist_group = rep(c("close", "far"), each = 4)
) %>%
  mutate(dist_treatment = factor(dist_treatment, levels = c("bracelet", "calendar", "ink", "control")))
cluster_mapper = analysis_data %>%
  select(
    cluster_id,
    assigned_treatment,
    assigned_dist_group
  ) %>% unique()

roc_dist_idx_mapper = tibble(
    roc_distance = seq(0, 5000, 100)) %>% 
    mutate(roc_dist_idx = seq(n()))


        
fit_type_str = if_else(script_options$prior, "prior", "fit")
if (length(script_options$chain) > 1) {
  chain_str = str_glue("{min(script_options$chain)}-{max(script_options$chain)}")
} else {
  chain_str = script_options$chain
}


struct_vstar_draws_raw = load_param_draws(
    fit_version = script_options$struct_fit_version,
    model = script_options$struct_model,
    chain = script_options$chain,
    prior_predictive = script_options$prior,
    input_path = script_options$input_path,
    structural_cluster_obs_v[cluster_id],
    total_error_sd[cluster_id]
)




struct_vstar_draws_raw = struct_vstar_draws_raw %>%
  left_join(
    cluster_mapper,
    by =  "cluster_id"
  )  %>%
  mutate(
    structural_cluster_obs_v = -1*structural_cluster_obs_v / total_error_sd # to make comparable to linear predictor in brms/RF probit
  ) %>%
  select(-total_error_sd) %>%
  left_join(
    analysis_data %>%
        select(
            cluster_id,
            standard_cluster.dist.to.pot
        ) %>%
        unique()
  )

rf_dist_output =  readRDS(file.path(script_options$input_path, paste0(script_options$rf_model, ".rds")))
rf_nodist_output =  readRDS(file.path(script_options$input_path, paste0("brm_rf_fit_hier_cluster_fit", ".rds")))

rf_nodist_model = rf_nodist_output$brm_model
rf_nodist_data = rf_nodist_output$data %>%
    select(-cluster_id) %>%
    rename(cluster_id = cluster.id)
rf_dist_model = rf_dist_output$brm_model
rf_dist_data = rf_dist_output$data %>%
    select(-cluster_id) %>%
    rename(cluster_id = cluster.id)

rf_nodist_vstar_draws_raw = tidybayes::linpred_rvars(
    rf_nodist_model,
    rf_nodist_data %>% 
        mutate(cluster.id = cluster_id) %>%
        group_by(
            assigned_treatment, assigned_dist_group, 
            cluster_id, cluster.id, county
        ) %>%
        summarise(
            standard_cluster.dist.to.pot = mean(standard_cluster.dist.to.pot)
        ),
    ndraws = 400
) %>%
    ungroup() %>%
    select(
        assigned_treatment, assigned_dist_group, 
        cluster_id,
        structural_cluster_obs_v = .linpred
    ) %>%
    group_by(cluster_id) %>%
    mutate(
        cluster_id = cur_group_id()
    ) %>%
    ungroup()

rf_dist_vstar_draws_raw = tidybayes::linpred_rvars(
    rf_dist_model,
    rf_dist_data %>% 
        mutate(cluster.id = cluster_id) %>%
        group_by(
            assigned_treatment, assigned_dist_group, 
            cluster_id, cluster.id, county
        ) %>%
        summarise(
            standard_cluster.dist.to.pot = mean(standard_cluster.dist.to.pot)
        ),
    ndraws = 400
) %>%
    ungroup() %>%
    select(
        assigned_treatment, assigned_dist_group, 
        cluster_id,
        standard_cluster.dist.to.pot,
        structural_cluster_obs_v = .linpred
    ) %>%
    group_by(cluster_id) %>%
    mutate(
        cluster_id = cur_group_id()
    ) %>%
    ungroup()



rf_vstar_draws_raw = bind_rows(
    rf_nodist_vstar_draws_raw %>% mutate(model = "reduced form - no distance"),
    rf_dist_vstar_draws_raw %>% mutate(model = "reduced form - distance")
)

vstar_draws = bind_rows(
    struct_vstar_draws_raw %>% mutate(model = "structural"),
    rf_vstar_draws_raw
)


vstar_draws %>%
    mutate(
        e_v = mean(structural_cluster_obs_v)
    ) %>%
    filter(model != "reduced form - no distance") %>%
    ggplot(aes(
        x = e_v,
        y = interaction(assigned_treatment),
        fill = model
    ), alpha = 0.1) +
    ggridges::geom_density_ridges(alpha = 0.5) +
    ggridges::theme_ridges() +
    theme(legend.position = "bottom") +
    facet_wrap(~assigned_dist_group) +
    ggthemes::scale_fill_canva("", palette = "Primary colors with a vibrant twist") 

vstar_draws %>%
    group_by(
        model,
        assigned_treatment,
        assigned_dist_group
    ) %>%
    mutate(
        within_var = var(structural_cluster_obs_v),
        across_var = rvar_var(structural_cluster_obs_v)
    ) %>%
    select(
        model,
        cluster_id,
        assigned_treatment,
        assigned_dist_group,
        within_var,
        across_var
    )  %>%
    ggplot(aes(
        xdist = across_var,
        y = 
    ))
    
    
    
    
    
    


vstar_stats = vstar_draws %>%
    group_by(
        model,
        assigned_treatment,
        assigned_dist_group
    ) %>%
    mutate(
        within_var = var(structural_cluster_obs_v),
        across_var = rvar_var(structural_cluster_obs_v)
    ) %>%
    select(
        model,
        cluster_id,
        assigned_treatment,
        assigned_dist_group,
        within_var,
        across_var
    ) %>%
    group_by(
        model,
        assigned_treatment,
        assigned_dist_group
    )  %>%
    mutate(
        heterogeneity_frac = across_var/(within_var + across_var)
    ) %>%
    summarise(
        mean_within_var = mean(within_var),
        across_var = E(rvar_mean(across_var)),
        heterogeneity_frac = E(rvar_mean(heterogeneity_frac))
    ) %>% 
    pivot_wider(
        names_from = model,
        values_from = c(mean_within_var, across_var, heterogeneity_frac)
    ) %>%
    select(
        assigned_treatment,
        assigned_dist_group,
        contains("across"),
        contains("within"),
        contains("heterogeneity")
        )   

vstar_stats %>%
    select(-contains("no distance")) %>%
    View()

max_dist = max(vstar_draws$standard_cluster.dist.to.pot, na.rm = TRUE)*1.05

overall_mean_vstar = vstar_draws %>%
    group_by(
        model,
        assigned_treatment,
        assigned_dist_group
    ) %>%
    summarise(
        mean_v = E(rvar_mean(structural_cluster_obs_v)),
    ) %>%
    mutate(
        standard_cluster.dist.to.pot = if_else(
            assigned_dist_group == "close",
            0,
            max_dist
        )
    )


summ_vstar_draws =  vstar_draws %>%
    median_qi(
        structural_cluster_obs_v
    )  %>%
    to_broom_names()

stop()



summ_vstar_draws %>%
    filter(model != "reduced form - no distance")  %>%
    ggplot(aes(
        x = structural_cluster_obs_v,
        xmin = conf.low,
        xmax = conf.high,
        y = standard_cluster.dist.to.pot,
        colour = assigned_dist_group
    )) +
    facet_grid(model ~ assigned_treatment)  +
    geom_pointrange()  +
    geom_point(
        data = overall_mean_vstar %>% filter(model != "reduced form - no distance"),
        inherit.aes = FALSE,
        aes(
            x = mean_v,
            y = standard_cluster.dist.to.pot,
            group = assigned_dist_group,
            colour = assigned_dist_group
        ),
        shape = 4
    ) +
    theme_minimal() +
    theme(legend.position = "none")


vstar_udraws = vstar_draws %>%
    unnest_rvars()

vstar_udraws %>%
filter(model %in% c("reduced form - distance", "structural")) %>%
  mutate(
    assigned_dist_group = str_to_title(assigned_dist_group),
    assigned_treatment = str_to_title(assigned_treatment)
  ) %>%
  ggplot(aes(
    x = structural_cluster_obs_v
  )) +
  # geom_histogram(position = position_dodge(0.5), bins = 60) +
  geom_histogram(
    aes(
      structural_cluster_obs_v, 
      y = after_stat(ndensity), 
      fill = model)
      , 
        position = "dodge", alpha = 0.75) +
  facet_grid(assigned_treatment ~ assigned_dist_group)  +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ggthemes::scale_fill_canva("", palette = "Primary colors with a vibrant twist") +
  labs(
    y = "", x = latex2exp::TeX("$w*$")
  ) 
