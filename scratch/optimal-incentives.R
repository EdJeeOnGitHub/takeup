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
                            bracelet
                            bracelet
                            --output-name=optimal-incentives-TEST-b-bracelet-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
                            --from-csv
                            --num-post-draws=200
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


first_param_draw = extract_params(
    param_draws = struct_param_draws,
    private_benefit_treatment = script_options$private_benefit_z,
    visibility_treatment = script_options$visibility_z,
    draw_id = 2,
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


#' distance in meters
find_optimal_incentive = function(distance, lambda, params) {
    takeup_fun = find_pred_takeup(params)
    takeup_list = takeup_fun(distance)

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

    lhs = (delta - mu_rep_deriv * delta_v_star) * (v_star + b + lambda * delta * dist_norm) / (1 + mu_rep * delta_v_star_deriv)

    rhs = lambda * delta / hazard

    diff = abs(lhs - rhs)
    return(diff)
}

find_optimal_incentive(100, 0.1, first_param_draw)
anon_f = function(x, lambda) {
    find_optimal_incentive(x, lambda, first_param_draw)
}

lambda = seq(from = 0, to = 2, length.out = 20)

lambda_df = tibble(
    lambda = lambda
)

lambda_df = lambda_df %>%
    mutate(
        lambda_funs = map(
            lambda,
            ~function(x) anon_f(x, lambda = .x)
        ),
        lambda_res = map(
            lambda_funs,
            ~optim(1, .x, method = "Brent", lower = 0, upper = 5000)
        )
    )

lambda_df %>%
    mutate(
        res = map_dbl(lambda_res, "par")
    ) %>%
    select(lambda, res) %>%
    ggplot(aes(
        x = lambda,
        y = res
    )) +
    geom_point()


me = lambda_df %>%
    pull(lambda_funs) %>%
    first()
me(10)
result = optim(1, anon_f, method = "Brent", lower = 0, upper = 5000)


result






