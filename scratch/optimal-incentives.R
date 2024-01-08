#!/usr/bin/Rscript
print(commandArgs(trailingOnly = TRUE))
script_options = docopt::docopt(
    stringr::str_glue("Usage:
    optimal-incentives.R <fit-version> <private-benefit-z> <visibility-z> [options] 

    Takes posterior fits and calculates optimal subsidy (PoT distance). There 
    are multiple counterfactuals:
        1. posterior - All posterior draws across all treatments - find Bayes optimal distance.
        2. robust-externality - Find optimal distance varying the value of externalities for a given treatment.
        3. robust-lambda - Find optimal distance varying cost of public funds for a given treatment.



    Options:
        --from-csv  Load parameter draws from a csv file
        --to-csv   Load stanfit object and write draws to csv
        --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
        --output-path=<path>  Path to find results [default: {file.path('temp-data')}]
        --output-name=<output-name>  Prepended to output file
        --num-cores=<num-cores>  Number of cores to use [default: 8]
        --model=<model>  Which model to use [default: STRUCTURAL_LINEAR_U_SHOCKS]
        --data-input-path=<data-input-path>  Where village and PoT data is stored [default: {file.path('optim', 'data')}]
        --data-input-name=<data-input-name>  Filename of village and PoT data [default: full-experiment-4-extra-pots.rds]
        --suppress-reputation  Suppress reputational returns
        --static-signal-pm  Use static signalling value at distance cutoff for all treatments
        --static-signal-distance=<static-signal-distance>  Distance to use for static signalling value [default: 1000]
        --single-chain  Only use first chain for draws (useful for debugging) 
        --fit-type=<fit-type>  Which fit type to use - prior predictive or posterior draws [default: fit] 
        --lambda=<lambda>  Lambda to use for optimal incentives across all treats [default: 0]
        --externality=<externality>  Externality to use across all treats [default: 0]
        --num-post-draws=<num-post-draws>  Number of posterior draws to use [default: 200]
        --posterior  Estimate across entire posterior across all treatments.
        --robust-externality  Estimate across externality grid for posterior median
        --robust-lambda  Estimate across lambda grid for posterior median
        --check-derivative-positive  Verify that dw/dc is positive for parameter values we see        


    "),
    args = if (interactive()) "
                            86
                            bracelet
                            bracelet
                            --output-name=ramsey-b-bracelet-mu-bracelet-lambda-0-externality-0-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
                            --num-post-draws=500
                            --num-cores=12
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP

                            --robust-externality
                            --robust-lambda
                            --posterior

                            --static-signal-pm
                            --static-signal-distance=500

                              " 
           else commandArgs(trailingOnly = TRUE)
)



set.seed(19484)

library(posterior)
library(tidyverse)
library(tidybayes)
library(broom)
library(nleqslv)
library(cmdstanr)
library(furrr)

source(file.path("optim", "optim-functions.R"))
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path("dist_structural_util.R"))

fit_version = script_options$fit_version

script_options$lambda = as.numeric(script_options$lambda)
script_options$num_post_draws = as.numeric(script_options$num_post_draws)
script_options$static_signal_distance = as.numeric(script_options$static_signal_distance)

struct_param_draws = read_csv(
    file.path(
        script_options$input_path, 
        str_interp(
            "param_posterior_draws_dist_${script_options$fit_type}${fit_version}_${script_options$model}.csv"
        )
    )
)




if (script_options$static_signal_pm) {
    post_output_name = str_glue(
        "posterior-optimal-distance-static_vstar-lambda_0-externality_0-model_{script_options$model}-fit-version_{script_options$fit_version}.csv"
    )
} else {
    post_output_name = str_glue(
        "posterior-optimal-distance-dynamic_vstar-lambda_{script_options$lambda}-externality_{script_options$externality}-model_{script_options$model}-fit-version_{script_options$fit_version}.csv"
    )
}

posterior_mean_draw = struct_param_draws %>%
    group_by(.variable, k, j) %>%
    summarise(.value = mean(.value)) %>% 
    mutate(.draw = 1)


script_options$num_post_draws = as.numeric(script_options$num_post_draws)

sd_of_dist = read_rds("temp-data/sd_of_dist.rds")

mu_rep_type = switch(
    script_options$model,
    STRUCTURAL_LINEAR_U_SHOCKS = 0,
    STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP = 1,
    STRUCTURAL_LINEAR_U_SHOCKS_LINEAR_MU_REP = 2,
    STRUCTURAL_LINEAR_U_SHOCKS_NO_REP = 3,
    STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP = 4
)


recalc_takeup = function(distance, params, b_add = 0, mu_add = 0) {
    # additional net benefit to get dewormed
    params$beta = params$beta  +  b_add
    params$centered_cluster_beta_1ord = params$centered_cluster_beta_1ord + mu_add
    takeup_fun = find_pred_takeup(params)
    takeup_list = takeup_fun(distance)
    mu_rep_0 = calculate_mu_rep(
                dist = 0,
                base_mu_rep = params$base_mu_rep,
                mu_beliefs_effect = 1,
                beta = params$centered_cluster_beta_1ord,
                dist_beta = params$centered_cluster_dist_beta_1ord,
                beta_control = params$mu_beta_z_control,
                dist_beta_control = params$mu_beta_d_control,
                mu_rep_type = params$mu_rep_type, 
                control = params$visibility_treatment == "control")
    pr_vis_0 = mu_rep_0 / params$base_mu_rep
    takeup_list$pr_vis_0 = pr_vis_0
    return(takeup_list)
}


#' distance in meters
find_optimal_incentive = function(distance, lambda, params, b_add = 0, mu_add = 0, externality = 0) {
    takeup_list = recalc_takeup(distance, params, b_add, mu_add)
    if (params$suppress_reputation) {
        takeup_list$mu_rep = 0
        takeup_list$mu_rep_deriv = 0
        takeup_list$delta_v_star = 0
        takeup_list$delta_v_star_deriv = 0
    }
    delta = takeup_list$delta
    mu_rep_deriv = takeup_list$mu_rep_deriv
    delta_v_star = takeup_list$delta_v_star
    mu_rep = takeup_list$mu_rep
    delta_v_star_deriv = takeup_list$delta_v_star_deriv
    v_star = takeup_list$v_star
    net_b = takeup_list$b # b is net benefit so add back
    dist_norm = distance / sd_of_dist
    b = net_b + delta * dist_norm 
    total_error_sd = takeup_list$total_error_sd

    hazard = dnorm(v_star/total_error_sd) / (1 - pnorm(v_star/total_error_sd))

    lhs = (-1) * (delta - mu_rep_deriv * delta_v_star) * (v_star + b + externality - lambda * delta * dist_norm) / (1 + mu_rep * delta_v_star_deriv)

    rhs = lambda * delta / hazard

    diff = abs(lhs - rhs)
    return(diff)
}



#### Estimation Across Posterior 
treatments = c(
    "control",
    "ink",
    "calendar",
    "bracelet"
)

if (script_options$posterior) {
# generate grid of posterior draw IDs and treatments
unique_post_draw_ids = unique(struct_param_draws$.draw) 
full_posterior_df = expand.grid(
    draw_id = sample(unique_post_draw_ids, min(script_options$num_post_draws, length(unique_post_draw_ids))),
    treatment = treatments
) %>%
    as_tibble()


if (script_options$static_signal_pm == TRUE) {
    # generate dynamic (non-static) prediction functions to calculate static signal 
    # at required cutoff
    true_pred_functions = map(
        full_posterior_df$draw_id,
        ~extract_params(
            param_draws = struct_param_draws,
            private_benefit_treatment = script_options$private_benefit_z,
            visibility_treatment = script_options$visibility_z,
            draw_id = .x,
            dist_sd = sd_of_dist,
            j_id = 1,
            rep_cutoff = script_options$rep_cutoff,
            dist_cutoff = script_options$dist_cutoff, 
            bounds = script_options$bounds,
            mu_rep_type = mu_rep_type,
            suppress_reputation = script_options$suppress_reputation, 
            static_signal = NA,
            fix_mu_at_1 = FALSE,
            fix_mu_distance = script_options$fix_mu_distance,
            static_delta_v_star = NA
        ) %>% find_pred_takeup()
    )


    static_pred_outputs = map(
        true_pred_functions, 
        ~.x(script_options$static_signal_distance)
    )

    delta_v_stars = map_dbl(static_pred_outputs, "delta_v_star")

    full_posterior_df$static_delta_v_star = delta_v_stars
}  else {
    full_posterior_df$static_delta_v_star = NA
}



# extract params for each posterior draw
full_posterior_df = full_posterior_df %>%
    mutate(
        params = pmap(
            list(
                draw_id,
                treatment
            ),
            ~extract_params(
                param_draws = struct_param_draws,
                private_benefit_treatment = ..2,
                visibility_treatment = ..2,
                draw_id = ..1,
                dist_sd = sd_of_dist,
                j_id = 1,
                rep_cutoff = Inf,
                dist_cutoff = Inf,
                bounds = c(-Inf, Inf),
                mu_rep_type = mu_rep_type,
                suppress_reputation = FALSE,
                static_signal = NA,
                fix_mu_at_1 = FALSE,
                fix_mu_distance = NULL,
                static_delta_v_star = NA
            )
        )
    )

find_optimal_static_delta_vstar_incentive = function(distance, params, static_delta_v_star) {
    dist_renorm = distance/params$dist_sd
    mu_rep = calculate_mu_rep(
                dist = dist_renorm,
                base_mu_rep = params$base_mu_rep,
                mu_beliefs_effect = 1,
                beta = params$centered_cluster_beta_1ord,
                dist_beta = params$centered_cluster_dist_beta_1ord,
                beta_control = params$mu_beta_z_control,
                dist_beta_control = params$mu_beta_d_control,
                mu_rep_type = params$mu_rep_type, 
                control = params$visibility_treatment == "control")


    lhs = static_delta_v_star*mu_rep/params$dist_beta_v 
    rhs = dist_renorm

    diff = abs(lhs - rhs)
    return(diff)
}


if (script_options$static_signal_pm) {
    # create anon functions for each param draw
    full_posterior_df = full_posterior_df %>%
        mutate(
            fun_treatments = pmap(
                list(params, static_delta_v_star),
                ~function(x) find_optimal_static_delta_vstar_incentive(
                    distance = x, 
                    params = ..1, 
                    static_delta_v_star = ..2
                    )
            )
            )
} else {
    # create anon functions for each param draw
    full_posterior_df = full_posterior_df %>%
        mutate(
            fun_treatments = map(
                params,
                ~function(x) find_optimal_incentive(
                    distance = x, 
                    lambda = script_options$lambda, 
                    params = .x, 
                    b_add = 0,
                    mu_add = 0,
                    externality = as.numeric(script_options$externality)
                    )
            )
            )
}
# optimise funs to get optimal distance
full_posterior_df = full_posterior_df %>%
    mutate(
        fit_funs = map(
            fun_treatments,
            ~optim(2500, .x, method = "Brent", lower = 0, upper = 20000),
            .progress = TRUE
        )
    )
# extract results
full_posterior_df = full_posterior_df %>%
    mutate(
        optimal_distance = map_dbl(fit_funs, "par")
    )
# save to csv
full_posterior_df %>%
    select(treatment, draw_id,  optimal_distance)  %>%
    write_csv(
    file.path(
        script_options$output_path,
        post_output_name
    )
    )
} # End posterior estimation

# Estimate Optimal Distance holding visibility fixed, only vary private incentive
b_df = expand.grid(
    draw = 1,
    lambda = script_options$lambda,
    b_add = seq(from = -3, to = 3, length.out = 20)
) %>% as_tibble() %>%
    group_by(draw) %>%
    nest(data = c(lambda, b_add))
b_df = b_df %>%
    mutate(
        params_vis = map(
            draw,
            ~extract_params(
                param_draws = posterior_mean_draw,
                private_benefit_treatment = script_options$private_benefit_z,
                visibility_treatment = script_options$visibility_z,
                draw_id = .x,
                dist_sd = sd_of_dist,
                j_id = 1,
                rep_cutoff = Inf,
                dist_cutoff = Inf, 
                bounds = c(-Inf, Inf),
                mu_rep_type = mu_rep_type,
                suppress_reputation = FALSE, 
                static_signal = NA,
                fix_mu_at_1 = FALSE,
                fix_mu_distance = NULL,
                static_delta_v_star = NA
            )
            ),
        params_control_vis = map(
            draw,
            ~extract_params(
                param_draws = posterior_mean_draw,
                private_benefit_treatment = script_options$private_benefit_z,
                visibility_treatment = "control",
                draw_id = .x,
                dist_sd = sd_of_dist,
                j_id = 1,
                rep_cutoff = Inf,
                dist_cutoff = Inf, 
                bounds = c(-Inf, Inf),
                mu_rep_type = mu_rep_type,
                suppress_reputation = FALSE, 
                static_signal = NA,
                fix_mu_at_1 = FALSE,
                fix_mu_distance = NULL,
                static_delta_v_star = NA
            )
            )
        ) 
b_df = b_df %>%
    unnest(data)

params_check = b_df  %>%
    pull(params_vis) %>%
    first() 
params_check$centered_cluster_dist_beta_1ord
params_check$visibility_treatment

b_df = b_df %>%
    mutate(
        funs_vis = pmap(
            list(
                lambda,
                params_vis,
                b_add
            ),
            ~function(x) find_optimal_incentive(distance = x, lambda = ..1, params = ..2, b_add = ..3)
        ),
        funs_control_vis = pmap(
            list(
                lambda,
                params_control_vis,
                b_add
            ),
            ~function(x) find_optimal_incentive(distance = x, lambda = ..1, params = ..2, b_add = ..3)
        )
    ) %>%
    ungroup()

b_df = b_df %>%
    ungroup() %>%
    mutate(
        fit_vis = map(
            funs_vis,
            ~optim(2500, .x, method = "Brent", lower = 0, upper = 20000),
            .progress = TRUE
        ),
        fit_control_vis = map(
            funs_control_vis,
            ~optim(2500, .x, method = "Brent", lower = 0, upper = 20000),
            .progress = TRUE
        )
    )

b_df = b_df %>%
    mutate(
        res_vis = map_dbl(fit_vis, "par"),
        res_control_vis = map_dbl(fit_control_vis, "par"),
        takeup_list = pmap(
            list(
                res_vis,
                params_vis,
                b_add
            ),
            ~recalc_takeup(distance = ..1, params = ..2, b_add = ..3, mu_add = 0)
        )
        ) %>%
    mutate(
        pred_takeup = map_dbl(
            takeup_list, 
            "pred_takeup"
            )
    )

b_df %>%
    select(
        lambda, b_add, params_vis, res_vis
    ) %>%
    mutate(
        ed = map2_dbl(
            b_add,
            params_vis, 
            ~analytical_delta(-.x, .y$total_error_sd)
            )
    )

b_df %>%
    select(
        lambda, b_add, params_vis, res_vis 
    ) %>%
    slice(1) %>%
    pull(params_vis) %>%
    first()




find_optimal_incentive_first_best = function(distance, delta, mu) {
    takeup_list = recalc_takeup(distance, params, b_add, mu_add)
    if (params$suppress_reputation) {
        takeup_list$mu_rep = 0
        takeup_list$mu_rep_deriv = 0
        takeup_list$delta_v_star = 0
        takeup_list$delta_v_star_deriv = 0
    }
    delta = takeup_list$delta
    mu_rep_deriv = takeup_list$mu_rep_deriv
    delta_v_star = takeup_list$delta_v_star
    mu_rep = takeup_list$mu_rep
    delta_v_star_deriv = takeup_list$delta_v_star_deriv
    v_star = takeup_list$v_star
    net_b = takeup_list$b # b is net benefit so add back
    dist_norm = distance / sd_of_dist
    b = net_b + delta * dist_norm 
    total_error_sd = takeup_list$total_error_sd

    hazard = dnorm(v_star/total_error_sd) / (1 - pnorm(v_star/total_error_sd))

    lhs = (-1) * (delta - mu_rep_deriv * delta_v_star) * (v_star + b - lambda * delta * dist_norm) / (1 + mu_rep * delta_v_star_deriv)

    rhs = lambda * delta / hazard

    diff = abs(lhs - rhs)
    return(diff)
}


b_df %>% 
    select(
        draw,
        lambda,
        b_add,
        contains("res"),
        takeup_list
    ) %>%
    mutate(
        pred_takeup = map_dbl(
            takeup_list, 
            "pred_takeup"
            )
    )


long_b_df = b_df %>%
    select(draw, lambda, b_add, pred_takeup, contains("res")) %>%
    pivot_longer(
        c(contains('res'), pred_takeup)
    ) %>%
    mutate(name = case_when(
        name == "res_control_vis" ~ "B: Bracelet, Mu: Control",
        name == "res_vis" ~ "B: Bracelet, Mu: Bracelet",
        name == "pred_takeup" ~ "Predicted Takeup"
        ))  %>%
    mutate(name_type = if_else(
        name == "Predicted Takeup",
        "pred_var",
        "res_var"
    ))

p_private_incentive_only = long_b_df %>%
    filter(name_type == "res_var") %>%
    ggplot(aes(
        x = b_add,
        y = value,
        colour = name,
        group = interaction(draw, name)
    )) +
    geom_line(size = 2) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        x = "Additional Private Incentive to Get Dewormed",
        y = "Optimal PoT Distance (meters)"
    ) +
    ggthemes::scale_color_canva("", palette = "Primary colors with a vibrant twist")  

p_private_incentive_only


ggsave(
    plot = p_private_incentive_only,
    file = file.path(
        script_options$output_path,
        str_glue("{script_options$output_name}-private-incentive-only.pdf")
    ),
    width = 8,
    height = 6
)


## Adding Variable Externalities
externality_fun = function(dist_m, b, dist_sd) {
    x = dist_m / dist_sd
    diff = 3*1000/dist_sd

    c = 2*b + atan(b/diff)*diff/2
    y = -1*atan(b/diff)*x + c
    return(y)
}
# check this works ed
find_optimal_incentive_varying_externality = function(distance, lambda, params, b_add = 0, mu_add = 0) {
    takeup_list = recalc_takeup(distance, params, b_add, mu_add)
    if (params$suppress_reputation) {
        takeup_list$mu_rep = 0
        takeup_list$mu_rep_deriv = 0
        takeup_list$delta_v_star = 0
        takeup_list$delta_v_star_deriv = 0
    }
    delta = takeup_list$delta
    mu_rep_deriv = takeup_list$mu_rep_deriv
    delta_v_star = takeup_list$delta_v_star
    mu_rep = takeup_list$mu_rep
    delta_v_star_deriv = takeup_list$delta_v_star_deriv
    v_star = takeup_list$v_star
    net_b = takeup_list$b # b is net benefit so add back
    dist_norm = distance / sd_of_dist
    b = net_b + delta * dist_norm 
    total_error_sd = takeup_list$total_error_sd

    hazard = dnorm(v_star/total_error_sd) / (1 - pnorm(v_star/total_error_sd))

    externality = externality_fun(
        distance,
        abs(params$beta_b_control),
        dist_sd = sd_of_dist
    )

    lhs = (-1) * (delta - mu_rep_deriv * delta_v_star) * (v_star + b + externality - lambda * delta * dist_norm) / (1 + mu_rep * delta_v_star_deriv)

    rhs = lambda * delta / hazard

    diff = abs(lhs - rhs)
    return(diff)
}

#### Varying B and Mu
seq_size = 30
b_mu_df = expand.grid(
    draw = 1,
    lambda = script_options$lambda,
    b_add = c(seq(from = -2, to = 2, length.out = seq_size), 0) %>% unique(),
    mu_add = c(seq(from = -4, to = 4, length.out = seq_size), 0) %>% unique()
) %>% as_tibble() %>%
    group_by(draw) %>%
    nest(data = c(lambda, b_add, mu_add))

b_mu_df = b_mu_df %>%
    mutate(
        params_vis = map(
            draw,
            ~extract_params(
                param_draws = posterior_mean_draw,
                private_benefit_treatment = script_options$private_benefit_z,
                visibility_treatment = script_options$visibility_z,
                draw_id = .x,
                dist_sd = sd_of_dist,
                j_id = 1,
                rep_cutoff = Inf,
                dist_cutoff = Inf, 
                bounds = c(-Inf, Inf),
                mu_rep_type = mu_rep_type,
                suppress_reputation = FALSE, 
                static_signal = NA,
                fix_mu_at_1 = FALSE,
                fix_mu_distance = NULL,
                static_delta_v_star = NA
            )
            )
        ) 

b_mu_df = b_mu_df %>%
    unnest(data)

b_mu_df = b_mu_df %>%
    mutate(
        funs_vis = pmap(
            list(
                lambda,
                params_vis,
                b_add,
                mu_add
            ),
            ~function(x) { find_optimal_incentive(distance = x*1000, lambda = ..1, params = ..2, b_add = ..3, mu_add = ..4) }
        ),
        funs_vary_externality = pmap(
            list(
                lambda,
                params_vis,
                b_add,
                mu_add
            ),
            ~function(x) {find_optimal_incentive_varying_externality(
                distance = x*1000,
                lambda = ..1,
                params = ..2,
                b_add = ..3,
                mu_add = ..4

            )}
        )
        )

b_mu_df = b_mu_df %>%
    ungroup() %>%
    mutate(
        fit_vis = map(
            funs_vis,
            ~optim(2.5, .x, lower = 0, upper = 15, method = "Brent"),
            .progress = TRUE
        ),
        fit_vary_ext_vis = map(
            funs_vis,
            ~optim(2.5, .x, lower = 0, upper = 15, method = "Brent")
        )
    )

b_mu_df = b_mu_df %>%
    mutate(
        res_vis = map_dbl(fit_vis, "par"),
        res_vary_ext = map_dbl(fit_vary_ext_vis, "par")
        )

b_mu_df = b_mu_df %>%
    mutate(
        takeup_list = pmap(
            list(
                res_vis,
                params_vis,
                b_add,
                mu_add
            ),
            ~recalc_takeup(distance = ..1*1000, params = ..2, b_add = ..3, mu_add = ..4)
        ),
        ext_takeup_list = pmap(
            list(
                res_vary_ext,
                params_vis,
                b_add,
                mu_add
            ),
            ~recalc_takeup(distance = ..1*1000, params = ..2, b_add = ..3, mu_add = ..4)
        )
    )  %>%
    mutate(
        pr_obs = map_dbl(takeup_list, "pr_vis_0"),
        ext_pr_obs = map_dbl(ext_takeup_list, "pr_vis_0")
    )

p_contour = b_mu_df %>%
    select(
        draw,
        lambda,
        b_add,
        pr_obs,
        res_vis
    ) %>%
    ggplot(aes(
        x = b_add,
        y = pr_obs,
        z = res_vis
    )) +
    metR::geom_contour_fill() +
    theme_minimal() +
    scale_fill_viridis_c(
         option = "inferno"
    )  +
    labs(
        x = "Shift in Norms/Additional Private Incentive",
        y = "Visibility (%)",
        fill = "Optimal Distance (km)"
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme(
        legend.position = "bottom"
    )

ggsave(
    p_contour,
    filename = file.path(
        script_options$output_path,
        str_glue("{script_options$output_name}-contour-plot.pdf")
    ),
    width = 10,
    height = 10
)

p_ext_contour = b_mu_df %>%
    select(
        draw,
        lambda,
        b_add,
        ext_pr_obs,
        res_vary_ext
    ) %>%
    ggplot(aes(
        x = b_add,
        y = ext_pr_obs,
        z = res_vary_ext
    )) +
    metR::geom_contour_fill() +
    theme_minimal() +
    scale_fill_viridis_c(
         option = "inferno"
    )  +
    labs(
        x = "Shift in Norms/Additional Private Incentive",
        y = "Visibility (%)",
        fill = "Optimal Distance (km)"
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme(
        legend.position = "bottom"
    )

ggsave(
    p_ext_contour,
    filename = file.path(
        script_options$output_path,
        str_glue("{script_options$output_name}-externality-contour-plot.pdf")
    ),
    width = 10,
    height = 10
)


## Varying Lambda
if (script_options$robust_lambda) {
lambda_df = expand.grid(
    draw = 1,
    lambda = seq(from = 0, to = 0.15, length.out = seq_size),
    b_add = seq(from = -2, to = 2, length.out = seq_size),
    mu_add = 0
) %>% as_tibble() %>%
    group_by(draw) %>%
    nest(data = c(lambda, b_add, mu_add))

lambda_df = lambda_df %>%
    mutate(
        params_vis = map(
            draw,
            ~extract_params(
                param_draws = posterior_mean_draw,
                private_benefit_treatment = script_options$private_benefit_z,
                visibility_treatment = script_options$visibility_z,
                draw_id = .x,
                dist_sd = sd_of_dist,
                j_id = 1,
                rep_cutoff = Inf,
                dist_cutoff = Inf, 
                bounds = c(-Inf, Inf),
                mu_rep_type = mu_rep_type,
                suppress_reputation = FALSE, 
                static_signal = NA,
                fix_mu_at_1 = FALSE,
                fix_mu_distance = NULL,
                static_delta_v_star = NA
            )
            )
        ) 

lambda_df = lambda_df %>%
    unnest(data)

lambda_df = lambda_df %>%
    mutate(
        funs_vis = pmap(
            list(
                lambda,
                params_vis,
                b_add,
                mu_add
            ),
            ~function(x) { find_optimal_incentive(distance = x*1000, lambda = ..1, params = ..2, b_add = ..3, mu_add = ..4) }
        )
        )

lambda_df = lambda_df %>%
    ungroup() %>%
    mutate(
        fit_vis = map(
            funs_vis,
            ~optim(
                2.5, 
                .x, 
                lower = 0, 
                upper = 18, 
                method = "Brent",
                control = list(
                    abstol = 1e-12,
                    reltol = 1e-12
                    )
                ),
            .progress = TRUE
        )
    )

lambda_df = lambda_df %>%
    mutate(
        res_vis = map_dbl(fit_vis, "par")
        )


p_lambda = lambda_df %>%
    select(draw, lambda, b_add, contains("res")) %>%
    pivot_longer(
        contains('res')
    ) %>%
    # this seems to be doing v weird stuff:
    # maybe the multiple equilibria thing?
    filter(value > 1) %>%
    ggplot(aes(
        x = b_add,
        y = value,
        colour = lambda,
        group = lambda
    )) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        x = "Shift in Norms/Additional Private Incentive",
        y = "Optimal Distance (km)",
        colour = latex2exp::TeX("\\lambda")
    )


ggsave(
    plot = p_lambda,
    file = file.path(
        script_options$output_path,
        str_glue("{script_options$output_name}-lambda-control.pdf")
    ),
    width = 8,
    height = 6
)


}


## Varying Externality
if (script_options$robust_externality) {
## Externalities ... and repeated code
externality_df = expand.grid(
    draw = 1,
    externality = seq(from = 0, to = abs(params_check$beta_b_control), length.out = seq_size),
    b_add = seq(from = -2, to = 2, length.out = seq_size),
    mu_add = 0,
    lambda = 0
) %>% as_tibble() %>%
    group_by(draw) %>%
    nest(data = c(externality, lambda, b_add, mu_add))

params_check$beta_b_control

params_check
externality_df = externality_df %>%
    mutate(
        params_vis = map(
            draw,
            ~extract_params(
                param_draws = posterior_mean_draw,
                private_benefit_treatment = script_options$private_benefit_z,
                visibility_treatment = script_options$visibility_z,
                draw_id = .x,
                dist_sd = sd_of_dist,
                j_id = 1,
                rep_cutoff = Inf,
                dist_cutoff = Inf, 
                bounds = c(-Inf, Inf),
                mu_rep_type = mu_rep_type,
                suppress_reputation = FALSE, 
                static_signal = NA,
                fix_mu_at_1 = FALSE,
                fix_mu_distance = NULL,
                static_delta_v_star = NA
            )
            )
        ) 

externality_df = externality_df %>%
    unnest(data)

externality_df = externality_df %>%
    mutate(
        funs_vis = pmap(
            list(
                lambda,
                params_vis,
                b_add,
                mu_add,
                externality
            ),
            ~function(x) { 
                find_optimal_incentive(
                    distance = x*1000, 
                    lambda = ..1, 
                    params = ..2, 
                    b_add = ..3, 
                    mu_add = ..4,
                    externality = ..5
                    ) 
                }
        )
        )

externality_df = externality_df %>%
    ungroup() %>%
    mutate(
        fit_vis = map(
            funs_vis,
            ~optim(
                2.5, 
                .x, 
                lower = 0, 
                upper = 18, 
                method = "Brent",
                control = list(
                    abstol = 1e-12,
                    reltol = 1e-12
                    )
                ),
                .progress = TRUE
        )
    )

externality_df = externality_df %>%
    mutate(
        res_vis = map_dbl(fit_vis, "par")
        )


p_externality = externality_df %>%
    select(draw, externality, lambda, b_add, contains("res")) %>%
    pivot_longer(
        contains('res')
    ) %>%
    mutate(externality_pct = 100*externality/abs(params_check$beta_b_control)) %>%
    # this seems to be doing v weird stuff:
    # maybe the multiple equilibria thing?
    # filter(value > 1) %>%
    ggplot(aes(
        x = b_add,
        y = value,
        colour = externality_pct,
        group = externality_pct
    )) +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        x = "Shift in Norms/Additional Private Incentive",
        y = "Optimal Distance (km)",
        colour = "Value of Externality (Relative to Deworming Private Benefit, %)"
    )


ggsave(
    plot = p_externality,
    file = file.path(
        script_options$output_path,
        str_glue("{script_options$output_name}-externality-control.pdf")
    ),
    width = 8,
    height = 6
)
}



if (script_options$check_derivative_positive) {
anon_takeup_fun = function(distance, params) {
    takeup_fun = find_pred_takeup(params)
    takeup_output = takeup_fun(distance)
    return(takeup_output)
}



unique_post_draw_ids = unique(struct_param_draws$.draw) 
dist_posterior_df = expand.grid(
    draw_id = sample(unique_post_draw_ids, min(100, length(unique_post_draw_ids))),
    distance = c(
        seq(from = 0, to = 2500, length.out = 20),
        3000,
        4000,
        5000
        ),
    treatment = treatments
) %>%
    as_tibble()


dist_posterior_df = dist_posterior_df %>%
    mutate(
        params = pmap(
            list(
                draw_id,
                treatment
            ),
            ~extract_params(
                param_draws = struct_param_draws,
                private_benefit_treatment = ..2,
                visibility_treatment = ..2,
                draw_id = ..1,
                dist_sd = sd_of_dist,
                j_id = 1,
                rep_cutoff = Inf,
                dist_cutoff = Inf,
                bounds = c(-Inf, Inf),
                mu_rep_type = mu_rep_type,
                suppress_reputation = FALSE,
                static_signal = NA,
                fix_mu_at_1 = FALSE,
                fix_mu_distance = NULL,
                static_delta_v_star = NA
            ),
            .progress = TRUE
        )
    )

dist_posterior_df = dist_posterior_df %>%
    mutate(
        takeup_output = map2(
            distance,
            params,
            ~anon_takeup_fun(distance = .x, params = .y),
            .progress = TRUE
            ),
        takeup_output = map2(
            takeup_output,
            params,
            ~list_modify(.x, dist_delta = .y$dist_beta_v) 
        )
    )


calc_numerator = function(takeup_output) {
    dist_delta = takeup_output$dist_delta
    mu_rep_deriv = takeup_output$mu_rep_deriv
    mu_rep = takeup_output$mu_rep
    delta_v_star = takeup_output$delta_v_star
    delta_v_star_deriv = takeup_output$delta_v_star_deriv
    numerator = dist_delta - mu_rep_deriv*delta_v_star
    return(numerator)
}

calc_denominator = function(takeup_output) {
    dist_delta = takeup_output$dist_delta
    mu_rep_deriv = takeup_output$mu_rep_deriv
    mu_rep = takeup_output$mu_rep
    delta_v_star = takeup_output$delta_v_star
    delta_v_star_deriv = takeup_output$delta_v_star_deriv
    denominator = 1 + mu_rep*delta_v_star_deriv
    return(denominator)
}

dist_posterior_df = dist_posterior_df %>%
    mutate(
        numerator = map_dbl(takeup_output, calc_numerator),
        denominator = map_dbl(takeup_output, calc_denominator)
    )


p_deriv = dist_posterior_df %>%
    mutate(
        frac = numerator / denominator,
        treatment = factor(treatment, levels = c("control", "ink", "calendar", "bracelet")),
        treatment = fct_relabel(treatment, str_to_title)
    ) %>%
    ggplot(aes(
        x = frac,
        fill = treatment
    )) +
    geom_histogram() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    geom_vline(xintercept = 0, linetype = "longdash") +
    labs(
        y = "Count",
        x = "Derivative Value"
    ) +
    ggthemes::scale_fill_canva("", palette = "Primary colors with a vibrant twist")

ggsave(
    plot = p_deriv,
    filename = file.path(
        "temp-data",
        "dydc-derivative-positive-check.pdf"
    ),
    width = 8,
    height = 6
)

}


#### Ed Externalities


df = tibble(
    x = seq(0, 6000, length.out = 500)
) %>%
    mutate(y = map_dbl(x, ~externality_fun(.x, abs(params_check$beta_b_control), sd_of_dist)))
df

df %>%
    ggplot(aes(
        x = x,
        y = y
    )) +
    geom_point()


df %>%
    group_by(x >= 3) %>%
    summarise(
        mean_y = mean(y)
    )
