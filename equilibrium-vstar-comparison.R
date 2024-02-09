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
  --output-path=temp-data/w_distribution
  --struct-model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_HIER_FOB
  --rf-model=brm_flat_dist_fit
  1
  " else commandArgs(trailingOnly = TRUE)
)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(brms)
library(ggthemes)
library(latex2exp)


dir.create(script_options$output_path, showWarnings = FALSE)

source("quick_postprocess_functions.R")

canva_pal = "Primary colors with a vibrant twist"

dist_sd = readRDS("temp-data/sd_of_dist.rds")
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
) %>%
    mutate(
        assigned_treatment = factor(assigned_treatment, levels = c("control", "ink", "calendar", "bracelet")) ,
        assigned_treatment = fct_relabel(assigned_treatment, str_to_title) %>% fct_rev 
    )


vstar_draws = vstar_draws %>%
    mutate(
        cluster_mean = mean(structural_cluster_obs_v),
        demean_w_cluster = structural_cluster_obs_v - cluster_mean 
    ) 


create_within_across_sigma_wstar = function(draws) {
    var_within = draws %>%
        group_by(model) %>%
        mutate(
            posterior_variance_within_cluster = var(structural_cluster_obs_v)
        )
    var_across = draws %>%
        group_by(model) %>%
        summarise(
            variance_across_clusters = rvar_var(structural_cluster_obs_v)
        ) %>%
        unnest_rvars()
    comp_w_draws = bind_rows(
        var_within %>%
            select(model,  sigma_w = posterior_variance_within_cluster)  %>%
            mutate(type = "Within"),
        var_across %>%
            select(model, sigma_w = variance_across_clusters)  %>%
            mutate(type = "Across")
    ) %>%
    ungroup()
    return(comp_w_draws)
}

comp_var_df = create_within_across_sigma_wstar(vstar_draws)


p_sigma_comp = comp_var_df %>%
    filter(model != "reduced form - no distance") %>%
    mutate(
        model = case_when(
            model == "reduced form - distance" ~ "Reduced Form",
            model == "structural" ~ "Structural"
        )
    ) %>%
    ggplot(aes(
        x = sigma_w,
        fill = model
    )) +
    geom_density(
        alpha = 0.5
    ) +
    facet_wrap(~type, scales = "free") + 
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_canva(
        "", 
        palette = canva_pal
    ) +
    labs(
        x = TeX('$\\sigma_{w^*}$')
    )



ggsave(
    plot = p_sigma_comp,
    file.path(
       script_options$output_path, 
       "sigma_comp.pdf"
    ),
    width = 8,
    height = 6
)


p_sigma_comp_hist = comp_var_df %>%
    filter(model != "reduced form - no distance") %>%
    mutate(
        model = case_when(
            model == "reduced form - distance" ~ "Reduced Form",
            model == "structural" ~ "Structural"
        )
    ) %>%
    ggplot(aes(
        x = sigma_w,
        fill = model
    )) +
    geom_histogram(
        alpha = 0.5,
        colour = "black",
        bins = 30,
        position = position_dodge()
    ) +
    facet_wrap(~type, scales = "free") + 
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_canva(
        "", 
        palette = canva_pal
    ) +
    labs(
        x = TeX('$\\sigma_{w^*}$')
    )

p_sigma_comp_hist

ggsave(
    plot = p_sigma_comp_hist,
    file.path(
       script_options$output_path, 
       "sigma_comp_hist.pdf"
    ),
    width = 8,
    height = 6
)

cluster_mean_draws = vstar_draws %>%
    filter(model != "reduced form - no distance") %>%
    select(
        model,
        assigned_treatment,
        assigned_dist_group,
        cluster_id,
        cluster_mean,
        standard_cluster.dist.to.pot
    ) 
within_cluster_draws = vstar_draws %>%
    select(
        model,
        assigned_treatment,
        assigned_dist_group,
        cluster_id,
        demean_w_cluster,
        standard_cluster.dist.to.pot
    ) %>%
    filter(model != "reduced form - no distance") %>%
    unnest_rvars() %>%
    ungroup()

comp_w_draws = bind_rows(
    within_cluster_draws %>% 
        rename(w_star = demean_w_cluster) %>%
        mutate(type = "Within Cluster"),
    cluster_mean_draws %>%
        rename(w_star = cluster_mean) %>%
        mutate(type = "Across Cluster")
)


comp_w_draws %>%
    ggplot(aes(
        x = w_star,
        fill = model
    )) +
    geom_density(alpha = 0.5)  +
    theme_minimal() +
    theme(legend.position = "bottom")  +
    facet_wrap(~type, scales = "free") +
    scale_fill_canva(
        "", 
        palette = canva_pal
        )  



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
    ) %>%
    mutate(
        model = case_when(
            model == "reduced form - distance" ~ "Reduced Form",
            model == "structural" ~ "Structural"
        ),
        dist = standard_cluster.dist.to.pot * dist_sd
    )


summ_vstar_draws =  vstar_draws %>%
    median_qi(
        structural_cluster_obs_v
    )  %>%
    to_broom_names() %>%
    mutate(
        model = case_when(
            model == "reduced form - distance" ~ "Reduced Form",
            model == "structural" ~ "Structural"
        )
    ) %>%
    mutate(
        dist = standard_cluster.dist.to.pot * dist_sd
    )




p_dist_plot = summ_vstar_draws %>%
    filter(model != "reduced form - no distance")  %>%
    ggplot(aes(
        x = structural_cluster_obs_v,
        y = dist/1000,
        colour = model
    )) +
    facet_grid(assigned_treatment~ model)  +
    geom_pointrange(aes(
        xmin = conf.low,
        xmax = conf.high
    ))  +
    geom_point(
        data = overall_mean_vstar %>% filter(model != "reduced form - no distance"),
        aes(
            x = mean_v,
            colour = model
        ),
        shape = 4
    ) +
    theme_minimal()  +
    ggthemes::scale_color_canva("", palette = "Primary colors with a vibrant twist")  +
    theme(legend.position = "bottom") +
    labs(
        x = TeX("$w*$"),
        y = "Distance (km)"
    )

p_dist_plot

ggsave(
    plot = p_dist_plot,
    file.path(
       script_options$output_path, 
       "vstar_dist_plot.pdf"
    ),
    width = 10,
    height = 10
)



p_dist_w = vstar_draws %>%
    select(
        model,
        fit_version,
        fit_type,
        assigned_treatment,
        assigned_dist_group,
        structural_cluster_obs_v
    ) %>%
    filter(
        model != "reduced form - no distance"
    ) %>%
    unnest_rvars() %>%
    ggplot(aes(
        x = structural_cluster_obs_v,
        fill = model
    )) +
    geom_histogram(
        alpha = 0.75,
        colour = "black",
        # bins = 30,
        position = position_dodge()
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_canva(
        "", 
        palette = canva_pal
    ) +
    labs(
        x = TeX('$w^*$')
    ) +
    facet_grid(
        assigned_treatment~assigned_dist_group
    ) 


ggsave(
    plot = p_dist_w,
    file.path(
       script_options$output_path, 
       "w_distribution.pdf"
    ),
    width = 10,
    height = 10
)
