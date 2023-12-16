#!/usr/bin/Rscript
print(commandArgs(trailingOnly = TRUE))
script_options = docopt::docopt(
    stringr::str_glue("Usage:
    info-test.R <fit-version>  [options] [ --from-csv | --to-csv ]

    Takes posterior fits and calculates estimated takeup for a given distance.


    Options:
        --from-csv  Load parameter draws from a csv file
        --to-csv   Load stanfit object and write draws to csv
        --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
        --output-path=<path>  Path to find results [default: {file.path('optim', 'data')}]
        --output-name=<output-name>  Prepended to output file
        --num-post-draws=<num-post-draws>  Number of posterior draws to use [default: 200]
        --num-cores=<num-cores>  Number of cores to use [default: 8]
        --model=<model>  Which model to use [default: STRUCTURAL_LINEAR_U_SHOCKS]
        --data-input-path=<data-input-path>  Where village and PoT data is stored [default: {file.path('optim', 'data')}]
        --data-input-name=<data-input-name>  Filename of village and PoT data [default: full-experiment-4-extra-pots.rds]
        --suppress-reputation  Suppress reputational returns
        --single-chain  Only use first chain for draws (useful for debugging) 
        --fit-type=<fit-type>  Which fit type to use - prior predictive or posterior draws [default: fit] 

    "),
    args = if (interactive()) "
                            86
                            --output-name=optimal-incentives-TEST-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
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
struct_param_draws

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
    draw = 1
    # treatment = script_options$treatment
) %>% as_tibble() 

incentives = c("ink", "bracelet")
params = map(
    incentives,
    ~extract_params(
        param_draws = struct_param_draws,
        private_benefit_treatment = .x,
        visibility_treatment = .x,
        draw_id = 1,
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

takeup_df = tibble(
    incentive = incentives,
    params = params
) %>%
    mutate(
        takeup_fun = map(params, find_pred_takeup),
        distance = list(seq(100, 2500, length.out = 20))
    ) %>%
    unnest(distance)

takeup_df = takeup_df %>%
    mutate(
        fit_fun = map2(takeup_fun, distance, ~.x(.y))
    )
takeup_df = takeup_df %>%
    mutate(
        pred_takeup = map_dbl(fit_fun, "pred_takeup")
    )
takeup_df %>%
    ggplot(aes(
        x = distance,
        y = pred_takeup,
        colour = incentive
    )) +
    geom_line()


pct_ink = 0.144 + 0.666
pct_bra = pct_ink + 0.022
pct_s_df = tibble(
    incentive = c("ink", "bracelet"),
    p_s_given_d = c(pct_ink, pct_bra)
)

takeup_df = takeup_df %>%
    left_join(pct_s_df, by = "incentive")



calculate_mi = function(p_d, p_s_given_d){
    p_s = p_s_given_d * p_d + (1 - p_s_given_d) * (1 - p_d)

    p_s_d = p_s_given_d * p_d
    p_sc_d = (1 - p_s_given_d) * p_d
    p_s_dc = 0
    p_sc_dc = 1 - p_s_d - p_sc_d - p_s_dc

    h_d = - p_d*log2(p_d) - (1 - p_d)*log2(1 - p_d)
    h_d_s = -p_s_d*log2(p_s_d/p_s) - p_sc_d*log2(p_sc_d/p_s) - p_s_dc - p_sc_dc*log2(p_sc_dc/p_s)

    mi = h_d - h_d_s
    return(
        lst(
            p_s,

            p_s_d,
            p_sc_d,
            p_s_dc,
            p_sc_dc,

            h_d,
            h_d_s,
            mi
        )
    )
}

info_df = takeup_df %>%
    select(
        incentive,
        distance, 
        pred_takeup,
        p_s_given_d
    )  %>%
    mutate(
        information = map2(pred_takeup, p_s_given_d, calculate_mi)
    ) %>%
    unnest_wider(information) %>%
    mutate(
        incentive = str_to_title(incentive)
    )

info_df %>%
    ggplot(aes(
        x = distance/1000,
        y = mi,
        colour = incentive
    ))  +
    geom_line() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    ggthemes::scale_color_canva(
        "",
        palette = "Primary colors with a vibrant twist"
    ) +
    labs(
        x = "Distance (km)",
        y = "Mutual Information (bits)"
    )

ggsave(
    "temp-plots/info-mi.pdf",
    width = 8,
    height = 6
)
stop()

p_grid = expand.grid(
    p_d = seq(0.01, 0.99, 0.01),
    p_s_given_d = seq(0.01, 0.99, 0.01)
) %>%
    as_tibble()

p_grid = p_grid %>%
    mutate(
        information =  map2(p_d, p_s_given_d, calculate_mi)
    )

p_grid = p_grid %>%
    unnest_wider(information)

p_grid %>%
    select(
        mi
    ) %>%
    filter(is.na(mi))

p_grid %>%
    select(p_d, p_s, mi)


p_grid %>%
    ggplot(aes(
        x = p_d,
        y = p_s_given_d,
        z = mi
    )) +
    metR::geom_contour_fill()


p_grid %>%
    pull(p_s_given_d) %>%
    unique()


p_grid %>%
    filter(p_d < 0.5) %>%
    filter(p_s_given_d > 0.5) %>%
    ggplot(aes(
        x = p_d,
        y = mi,
        colour = p_s_given_d,
        group = p_s_given_d
    )) +
    geom_line() 



log(0.5)

x = 0.3

x < 0.5


log(x) < log(0.5)

log(0.6)
