script_options = docopt::docopt(
    stringr::str_glue("Usage:
    analyse-reduced-form-sms.R  [options]

    Options:
        --load-fit
        --save-fit
        --fit-path=<fit-path>
        --fit-file=<fit-file>
        --include-reminder
        --output-path=<output-path>  Output path [default: temp-data]
    "),
    args = if (interactive()) "
        --load-fit
        --fit-path=data/stan_analysis_data
        --fit-file=SMS_BRMS_reminder_fit.rds
        --output-path=temp-data/sms-noreminder
        --include-reminder
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
monitored_sms_data <- analysis.data %>% 
  filter(mon_status == "monitored") %>% 
  filter(have_phone == "Yes") %>%
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_sms_data %>% 
    filter(have_phone == "Yes") %>%
    mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group, 
    sms_treatment = sms.treatment.2, 
    phone_owner = if_else(phone_owner == TRUE, "phone", "nophone"), 
    sms_treatment = str_replace_all(sms_treatment, "\\.", ""),
    female = gender == "female"
    ) %>%
    # reminder.only only present in control condition
    filter(phone_owner == "phone") %>%
    mutate(sms_treatment = factor(sms_treatment))

cluster_expected_dist_df = read_csv(here::here("data", "cluster_expected_dist.csv")) %>%
  rename(clust_expected_dist = dist) %>%
  mutate(cluster.id = cluster.id)

analysis_data = analysis_data %>%
    left_join(
        cluster_expected_dist_df,
        by = "cluster.id"
    )

l_cov_vars = c(
  "female",
  "age.census"
)

analysis_data %>%
    group_by(
        assigned_treatment, 
        sms_treatment
    ) %>%
    summarise(
        n = n()
    )



library(fixest)
rf_fit = feols(
    data = analysis_data,
    dewormed ~ assigned_treatment*assigned_dist_group*sms_treatment + .[l_cov_vars] + clust_expected_dist | county, 
    cluster =  ~cluster.id
    )



nobs(rf_fit)

analysis_data %>%
    count(
        sms_treatment
    )


breakdown_sms_comp = comparisons(
    rf_fit, 
    variables = "assigned_treatment",
    newdata = datagrid(
        assigned_dist_group = unique(analysis_data$assigned_dist_group), 
        sms_treatment = unique(analysis_data$sms_treatment))
) %>%
    as_tibble()

create_comp_dfs = function(fit, interval) {
    close_far_comp_df = avg_comparisons(
        fit,
        variable = "assigned_treatment",
        newdata = datagrid(
            assigned_dist_group = unique(analysis_data$assigned_dist_group),
            sms_treatment = unique(analysis_data$sms_treatment)), 
        conf.level = interval
    ) 

    new_low = paste0("conf.low_", (1 - interval)/2)
    new_high = paste0("conf.high_", 1 - (1 - interval)/2)
    comp_df = bind_rows(
        close_far_comp_df
    )  %>%
    mutate(interval = interval) %>%
    as_tibble()
    return(comp_df)

}


comp_df = avg_comparisons(
        rf_fit,
        variable = "assigned_treatment",
        newdata = datagrid(
            assigned_dist_group = unique(analysis_data$assigned_dist_group),
            assigned_treatment = unique(analysis_data$assigned_treatment),
            sms_treatment = unique(analysis_data$sms_treatment)
            ), 
        by = c("sms_treatment", "assigned_dist_group"),
        conf.level = 0.95
    ) 

comp_combined_df = avg_comparisons(
        rf_fit,
        variable = "assigned_treatment",
        newdata = datagrid(
            assigned_treatment = unique(analysis_data$assigned_treatment),
            sms_treatment = unique(analysis_data$sms_treatment)
            ), 
        by = "sms_treatment",
        conf.level = 0.95
    )

plot_df = bind_rows(
    comp_combined_df %>%
        tidy() %>%
        mutate(
            assigned_dist_group = "combined"
        ),
    comp_df %>%
        tidy() 
) %>%
   select(contrast, sms_treatment, assigned_dist_group, estimate, conf.low, conf.high) %>%
    mutate(
        lhs = str_extract(contrast, "(?<=\\()\\w+"),
        rhs = str_extract(contrast, "\\w+(?=\\)$)")
    ) %>%
    filter(rhs == "control") %>%
    filter(sms_treatment != "reminderonly")  %>%
    mutate(
        assigned_dist_group = factor(assigned_dist_group, levels = c("combined", "close", "far")),
        assigned_dist_group = fct_relabel(assigned_dist_group, str_to_title),
        lhs = factor(lhs, levels = c("ink", "calendar", "bracelet")),
        lhs = fct_relabel(lhs, str_to_title),
        sms_treatment = case_when(
            sms_treatment == "socialinfo" ~ "Social Info", 
            sms_treatment == "smscontrol" ~ "SMS Control"
        )
    )

p_sms_df = plot_df %>%
    ggplot(aes(
        x = estimate,
        xmin = conf.low,
        xmax = conf.high,
        y = lhs,
        colour = sms_treatment
    )) +
    facet_wrap(~assigned_dist_group, ncol = 1) +
    geom_pointrange(
        position = position_dodge(width = 0.3),
    )  +
    theme_minimal() +
    theme(legend.position = "bottom")  +
    labs(
        colour = "",
        y = "",
        x = ""
    ) +
    ggthemes::scale_color_canva("", palette = canva_palette_vibrant)   +
    geom_vline(xintercept = 0, linetype = "dotted")

ggsave(plot = p_sms_df, filename = file.path(script_options$output_path, "sms-TE-by-dist-incentive.pdf"), width = 8, height = 6)