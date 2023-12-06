#!/usr/bin/Rscript
print(commandArgs(trailingOnly = TRUE))
script_options = docopt::docopt(
    stringr::str_glue("Usage:
    optimal-incentives.R <fit-version> <private-benefit-z> <visibility-z> [options] [ --from-csv | --to-csv ]

    Takes posterior fits and calculates estimated takeup for a given distance.


    Options:
        --from-csv  Load parameter draws from a csv file
        --to-csv   Load stanfit object and write draws to csv
        --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
        --output-path=<path>  Path to find results [default: {file.path('temp-data')}]
        --output-name=<output-name>  Prepended to output file
        --num-post-draws=<num-post-draws>  Number of posterior draws to use [default: 200]
        --num-cores=<num-cores>  Number of cores to use [default: 8]
        --model=<model>  Which model to use [default: STRUCTURAL_LINEAR_U_SHOCKS]
        --data-input-path=<data-input-path>  Where village and PoT data is stored [default: {file.path('optim', 'data')}]
        --data-input-name=<data-input-name>  Filename of village and PoT data [default: full-experiment-4-extra-pots.rds]
        --suppress-reputation  Suppress reputational returns
        --single-chain  Only use first chain for draws (useful for debugging) 
        --fit-type=<fit-type>  Which fit type to use - prior predictive or posterior draws [default: fit] 
        --lambda=<lambda>  Lambda to use for optimal incentives [default: 0]
    "),
    args = if (interactive()) "
                            86
                            control
                            control
                            --output-name=ramsey-control-lambda-0-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
                            --from-csv
                            --num-post-draws=1
                            --num-cores=12
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
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
library(econometr)
library(furrr)

source(file.path("optim", "optim-functions.R"))
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path("dist_structural_util.R"))

fit_version = script_options$fit_version

script_options$lambda = as.numeric(script_options$lambda)

struct_param_draws = read_csv(
    file.path(
        script_options$input_path, 
        str_interp(
            "param_posterior_draws_dist_${script_options$fit_type}${fit_version}_${script_options$model}.csv"
        )
    )
)


struct_param_draws = struct_param_draws %>%
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



max_draw = max(struct_param_draws$.draw)
draw_treat_grid = expand.grid(
    draw = 1,
    # lambda = seq(from = 0, to = 0.3, length.out = 3),
    lambda = script_options$lambda,
    b_add = seq(from = -3, to = 3, length.out = 20)
    # treatment = script_options$treatment
) %>% as_tibble() %>%
    group_by(draw) %>%
    nest(data = c(lambda, b_add))

draw_treat_grid = draw_treat_grid %>%
    mutate(
        params_vis = map(
            draw,
            ~extract_params(
                param_draws = struct_param_draws,
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
                param_draws = struct_param_draws,
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
draw_treat_grid = draw_treat_grid %>%
    unnest(data)

param_check = draw_treat_grid  %>%
    pull(params_vis) %>%
    first() 
param_check$centered_cluster_dist_beta_1ord
param_check$visibility_treatment

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




draw_treat_grid = draw_treat_grid %>%
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

draw_treat_grid = draw_treat_grid %>%
    ungroup() %>%
    mutate(
        fit_vis = map(
            funs_vis,
            ~optim(2500, .x, method = "Brent", lower = 0, upper = 20000)
        ),
        fit_control_vis = map(
            funs_control_vis,
            ~optim(2500, .x, method = "Brent", lower = 0, upper = 20000)
        )
    )

draw_treat_grid = draw_treat_grid %>%
    mutate(
        res_vis = map_dbl(fit_vis, "par"),
        res_control_vis = map_dbl(fit_control_vis, "par")
        )

long_draw_treat = draw_treat_grid %>%
    select(draw, lambda, b_add, contains("res")) %>%
    pivot_longer(
        contains('res')
    ) %>%
    mutate(name = case_when(
        name == "res_control_vis" ~ "B: Bracelet, Mu: Control",
        name == "res_vis" ~ "B: Bracelet, Mu: Bracelet"
        )) 

p_private_incentive_only = long_draw_treat %>%
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

ggsave(
    plot = p_private_incentive_only,
    file = file.path(
        script_options$output_path,
        str_glue("{script_options$output_name}-private-incentive-only.pdf")
    ),
    width = 8,
    height = 6
)



#### B and Mu

seq_size = 30
b_mu_draw_treat_grid = expand.grid(
    draw = 1:min(max_draw, script_options$num_post_draws),
    lambda = script_options$lambda,
    b_add = seq(from = -2, to = 2, length.out = seq_size),
    # b_add = 0,
    mu_add = seq(from = -4, to = 4, length.out = seq_size)
    # treatment = script_options$treatment
) %>% as_tibble() %>%
    group_by(draw) %>%
    nest(data = c(lambda, b_add, mu_add))
b_mu_draw_treat_grid = b_mu_draw_treat_grid %>%
    mutate(
        params_vis = map(
            draw,
            ~extract_params(
                param_draws = struct_param_draws,
                private_benefit_treatment = "bracelet",
                visibility_treatment = "bracelet",
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

b_mu_draw_treat_grid = b_mu_draw_treat_grid %>%
    unnest(data)

b_mu_draw_treat_grid = b_mu_draw_treat_grid %>%
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

b_mu_draw_treat_grid = b_mu_draw_treat_grid %>%
    ungroup() %>%
    mutate(
        fit_vis = map(
            funs_vis,
            ~optim(2.5, .x, lower = 0, upper = 15, method = "Brent")
        )
    )

b_mu_draw_treat_grid = b_mu_draw_treat_grid %>%
    mutate(
        res_vis = map_dbl(fit_vis, "par")
        )

b_mu_draw_treat_grid = b_mu_draw_treat_grid %>%
    mutate(
        takeup_list = pmap(
            list(
                res_vis,
                params_vis,
                b_add,
                mu_add
            ),
            ~recalc_takeup(distance = ..1*1000, params = ..2, b_add = ..3, mu_add = ..4)
        )
    )  %>%
    mutate(
        pr_obs = map_dbl(takeup_list, "pr_vis_0")
    )

b_mu_draw_treat_grid  %>%
    select(b_add, pr_obs, res_vis) %>%
    tail()

p_contour = b_mu_draw_treat_grid %>%
    select(
        draw,
        lambda,
        b_add,
        pr_obs,
        res_vis
    ) %>%
    # filter(abs(b_add) < 2) %>%
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


library(plotly)
b_mat = b_mu_draw_treat_grid %>%
    select(b_add, pr_obs, res_vis) %>%
    pivot_wider(
        names_from = b_add,
        values_from = res_vis
    ) 

b_mat_mu = b_mat$pr_obs
b_mat_b = b_mat %>%
    select(-pr_obs) %>%
    colnames() %>%
    as.numeric()
b_mat_z = b_mat %>%
    select(-pr_obs) %>%
    as.matrix()
b_mat_z = b_mat_z * 1000

fig = plot_ly(
    x = b_mat_mu, 
    y = b_mat_b, 
    z = b_mat_z,
    contours = list(
    z = list(show = TRUE, start = min(b_mat_z), end = max(b_mat_z), size = 500))
    ) %>% add_surface()
fig_3d = fig %>%
    hide_colorbar() %>%
layout(
    showlegend = FALSE,
    scene = list(
      camera=list(
        eye = list(x=1.9/2, y=-4/2, z=2.5/2)
        ),
        yaxis = list(title = "Private Incentive"),
        xaxis = list(title = "Visibility"),
        zaxis = list(title = "Optimal Distance")
      )
  )
fig_3d

# orca(
#     fig_3d, 
#     scale = 3,
#     file = file.path(
#         script_options$output_path,
#         str_glue("{script_options$output_name}-3d-plot.pdf")
#         )
#     )



lambda_treat_grid = expand.grid(
    draw = 1:min(max_draw, script_options$num_post_draws),
    lambda = seq(from = 0, to = 0.15, length.out = seq_size),
    b_add = seq(from = -2, to = 2, length.out = seq_size),
    mu_add = 0
) %>% as_tibble() %>%
    group_by(draw) %>%
    nest(data = c(lambda, b_add, mu_add))

lambda_treat_grid = lambda_treat_grid %>%
    mutate(
        params_vis = map(
            draw,
            ~extract_params(
                param_draws = struct_param_draws,
                private_benefit_treatment = "control",
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

lambda_treat_grid = lambda_treat_grid %>%
    unnest(data)

lambda_treat_grid = lambda_treat_grid %>%
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

lambda_treat_grid = lambda_treat_grid %>%
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
                )
        )
    )

lambda_treat_grid = lambda_treat_grid %>%
    mutate(
        res_vis = map_dbl(fit_vis, "par")
        )


p_lambda = lambda_treat_grid %>%
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


## Externalities ... and repeated code
externality_treat_grid = expand.grid(
    draw = 1:min(max_draw, script_options$num_post_draws),
    externality = seq(from = 0, to = abs(param_check$beta_b_control), length.out = seq_size),
    b_add = seq(from = -2, to = 2, length.out = seq_size),
    mu_add = 0,
    lambda = 0
) %>% as_tibble() %>%
    group_by(draw) %>%
    nest(data = c(externality, lambda, b_add, mu_add))

externality_treat_grid = externality_treat_grid %>%
    mutate(
        params_vis = map(
            draw,
            ~extract_params(
                param_draws = struct_param_draws,
                private_benefit_treatment = "control",
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

externality_treat_grid = externality_treat_grid %>%
    unnest(data)

externality_treat_grid = externality_treat_grid %>%
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

externality_treat_grid = externality_treat_grid %>%
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
                )
        )
    )

externality_treat_grid = externality_treat_grid %>%
    mutate(
        res_vis = map_dbl(fit_vis, "par")
        )

abs(param_check$beta_b_control)

p_externality = externality_treat_grid %>%
    select(draw, externality, lambda, b_add, contains("res")) %>%
    pivot_longer(
        contains('res')
    ) %>%
    mutate(externality_pct = 100*externality/abs(param_check$beta_b_control)) %>%
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
        colour = "Value of Externality (Relative to Control, %)"
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
