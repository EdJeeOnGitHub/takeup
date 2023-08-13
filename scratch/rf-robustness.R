script_options = docopt::docopt(
    stringr::str_glue("Usage:
    rf-robustness.R  [options]

    Options:
        --load-fit
        --save-fit
        --fit-path=<fit-path>
        --fit-file=<fit-file>
        --include-reminder
        --output-path=<output-path>  Output path [default: temp-data]
    "),
    args = if (interactive()) "
    --save-fit
        --fit-path=data/stan_analysis_data
        --fit-file=REDUCED_FORM_ROBUSTNESS.rds
        --output-path=temp-data
    " else commandArgs(trailingOnly = TRUE)
    # args = if (interactive()) "takeup cv --models=REDUCED_FORM_NO_RESTRICT --cmdstanr --include-paths=stan_models --update --output-path=data/stan_analysis_data --outputname=test --folds=2 --sequential" else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(posterior)
library(tidybayes)
library(brms)
library(marginaleffects)
library(ggthemes)

options(mc.cores = 4)

canva_palette_vibrant <- "Primary colors with a vibrant twist"


source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

# Data --------------------------------------------------------------------

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)
# stick to monitored sms.treatment group
# remove sms.treatment.2
monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

monitored_sms_data <- analysis.data %>% 
  filter(mon_status == "monitored") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()



nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


analysis_data <- monitored_nosms_data %>% 
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group)


balance_data = read_rds(
  file.path(
    "temp-data",
    "saved_balance_data.rds"
  )
)

baseline_data = balance_data$baseline_balance_data

# - primary schooling
# - floor tile/cement
# - know everyone can be infected
# - maybe know biyearly

cluster_cov_df = baseline_data %>%
    group_by(cluster.id) %>%
    summarise(
        mean_primary = mean(completed_primary, na.rm = TRUE),
        mean_floor = mean(floor_tile_cement, na.rm = TRUE),
        mean_all_can_get_worms = mean(all_can_get_worms, na.rm = TRUE)
    )

analysis_data = analysis_data %>%
    left_join(
        cluster_cov_df,
        by = "cluster.id"
    ) %>%
    filter(!is.na(mean_primary)) %>%
    group_by(cluster.id) %>%
    mutate(
        std_cluster_dist_to_pot = mean(standard_cluster.dist.to.pot, na.rm = TRUE)
    )


default_priors = get_prior(
    data = analysis_data,
    dewormed ~ (assigned_treatment*assigned_dist_group | county + cluster.id) + 
        mean_primary + mean_floor + mean_all_can_get_worms + std_cluster_dist_to_pot, 
    family = bernoulli(link = "probit")
)


manual_priors = default_priors

# manual_priors$prior[9] = "normal(0, 0.1)"
# manual_priors$prior[10] = "normal(0, 1)"
# manual_priors$prior[18] = "normal(0, 0.1)"
# manual_priors$prior[19] = "normal(0, 1)"
manual_priors

if (script_options$load_fit) {
    rf_fit = read_rds(
        file.path(
            script_options$fit_path,
            script_options$fit_file
        )
    )
} else {

    rf_fit = brm(
        data = analysis_data,
        dewormed ~ (assigned_treatment*assigned_dist_group | county + cluster.id) + 
            mean_primary + mean_floor + mean_all_can_get_worms + std_cluster_dist_to_pot, 
        family = bernoulli(link = "probit"),
        prior = manual_priors,
        backend = "cmdstanr",
        silent = 0,
        refresh = 50,
        chains = 4,
        cores = 4,
        threads = threading(3),
    	control = list(adapt_delta = 0.99)
    )

}

if (script_options$save_fit) {
    saveRDS(
        rf_fit,
        file.path(
            script_options$fit_path,
            script_options$fit_file
        )
    )
}


# rf_priors = brm(
#     data = analysis_data,
#         dewormed ~ (assigned_treatment*assigned_dist_group | sms_treatment), 
#         family = bernoulli(link = "probit"),
#         sample_prior = "only", 
#         prior = manual_priors
# )

library(tidybayes)

treatments = c("control", "ink", "calendar", "bracelet")

pred_hat_dfs = map_dfr(
    treatments,
    ~analysis_data %>%
        select(
            cluster.id, 
            county, 
            assigned_treatment, 
            assigned_dist_group, 
            mean_primary, 
            mean_floor, 
            mean_all_can_get_worms,
            cluster_dist_to_pot
            ) %>%
        mutate(assigned_treatment = .x) %>%
        add_epred_rvars(rf_fit) %>%
        mutate(row_id = 1:n())
)

wide_pred_df = pred_hat_dfs %>%
    pivot_wider(
        values_from = `.epred`,
        id_cols = c(row_id, cluster.id, assigned_dist_group),
        names_from = assigned_treatment
    )

level_df = pred_hat_dfs %>%
    group_by(
        assigned_dist_group,
        assigned_treatment
    ) %>%
    summarise(
        level = rvar_mean(.epred)
    ) %>%
    bind_rows(
        .,
        pred_hat_dfs %>%
            group_by(assigned_treatment) %>%
            summarise(level = rvar_mean(.epred)) %>%
            mutate(assigned_dist_group = "combined")
    ) %>%
    mutate(assigned_dist_group = factor(assigned_dist_group, levels = c("far", "close", "combined")) %>%fct_rev)



canva_pal = "Primary colors with a vibrant twist"
level_df %>%
    median_qi(level, .width = 0.9) %>%
    to_broom_names() %>%
    mutate(assigned_treatment = factor(assigned_treatment, levels = treatments)) %>%
    ggplot(aes(
        x = level,
        xmin = conf.low,
        xmax = conf.high,
        colour = assigned_treatment,
        y = assigned_treatment
    )) +
    facet_wrap(~assigned_dist_group, ncol = 1) +
    geom_pointrange() +
    theme_minimal() +
    ggthemes::scale_color_canva(palette = canva_pal)  +
    guides(color = "none") +
    labs(
        x = "Takeup Level",
        y = "Treatment",
        title = "Levels - Controlling For Covariates"
    )

ggsave(
    "temp-plots/rf-robust-levels.pdf",
    width = 8,
    height = 6
)    


te_plot = level_df %>%
    mutate(
        te = level - level[assigned_treatment == "control"]
    ) %>%
    filter(assigned_treatment != "control") %>%
    median_qi(te, .width = 0.9) %>%
    to_broom_names() %>%
    mutate(assigned_treatment = factor(assigned_treatment, levels = treatments)) %>%
    ggplot(aes(
        x = te,
        xmin = conf.low,
        xmax = conf.high,
        y = assigned_treatment,
        colour = assigned_treatment
    )) +
    facet_wrap(~assigned_dist_group, ncol = 1) +
    geom_pointrange() +
    theme_minimal() +
    ggthemes::scale_color_canva(palette = canva_pal)  +
    guides(color = "none") +
    labs(
        x = "Takeup TE",
        y = "Treatment",
        title = "ATEs - Controlling For Covariates"
    ) +
    geom_vline(xintercept = 0, linetype = "longdash")

te_plot

ggsave(
    plot = te_plot,
    "temp-plots/rf-robust-te.pdf",
    width = 8,
    height = 6
)    