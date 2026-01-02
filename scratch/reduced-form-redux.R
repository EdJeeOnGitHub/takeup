library(tidyverse)
library(fixest)
library(ggthemes)
library(broom)



treat_levels_c = c("control", "ink", "calendar", "bracelet")
treat_levels = c("ink", "calendar", "bracelet")
dist_levels = c("close", "far")
model_level_order = c("reduced form", "structural")

quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)



canva_palette_vibrant <- "Primary colors with a vibrant twist"

theme_set(theme_minimal() +
            theme(legend.position = "bottom"))

source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))
source(file.path("multilvlr", "multilvlr_util.R"))

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

rct.schools.data <- read_rds(file.path("data", "takeup_rct_schools.rds"))
rct.cluster.selection <- read_rds(file.path("data", "rct_cluster_selection_2.0.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))
load(file.path("data", "takeup_village_pot_dist.RData"))

# Data -------------------------------------------------------------------------
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

create_community_variables = function(data) {
    data %>%
        group_by(cluster.id) %>%
        mutate(
            ave_floor = mean(floor == "Cement", na.rm = TRUE),
            have_primary = str_detect(school, "Secondary|College|University") | school == "Primary 8",
            ave_primary = mean(have_primary, na.rm = TRUE),
            n_indiv = replace_na(number_individuals, 99),
            n_indiv = factor(n_indiv),
            n_hhs = n()
        ) %>%
        ungroup() 
}



analysis_sms_data <- monitored_sms_data %>% 
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
  create_community_variables()


analysis_data <- monitored_nosms_data %>% 
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
  create_community_variables()



analysis_sms_data = analysis_sms_data %>%
    mutate(
        county = factor(county),
        cluster.id = factor(cluster.id),
        assigned_treatment = assigned.treatment,
        assigned_dist_group = dist.pot.group,
        signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
        signal = factor(signal, levels = c("no signal", "signal")),
        sms.treatment = replace_na(sms.treatment, "not in sms arm"),
        dist.to.pot = dist.to.pot / 1000,
        cluster.dist.to.pot = cluster.dist.to.pot / 1000
    )

analysis_data = analysis_data %>%
    mutate(
        county = factor(county),
        cluster.id = factor(cluster.id),
        assigned_treatment = assigned.treatment,
        assigned_dist_group = dist.pot.group,
        signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
        signal = factor(signal, levels = c("no signal", "signal")),
        dist.to.pot = dist.to.pot / 1000,
        cluster.dist.to.pot = cluster.dist.to.pot / 1000
    )

controls = c(
    "phone_owner", 
    "age.census", 
    "gender",
    # "cluster.dist.to.pot",
    # "dist.to.pot",
    "n_hhs",
    "ave_floor",
    "ave_primary",
    "assigned_dist_group"
    )



generate_fits = function(data, sms_data, controls, distance_var) {
    data = data %>%
        mutate(distance := {{ distance_var }})
    sms_data = sms_data %>%
        mutate(distance := {{ distance_var }})

    ols_fit = data %>%
        feols(
            data = .,
            fml = 
            dewormed ~ 0 + assigned_treatment +  distance + i(assigned_treatment, distance, ref = "control") | county,
            cluster = ~cluster.id
        ) 

    ols_sms_fit = sms_data %>%
        feols(
            data = .,
            fml = 
            dewormed ~ 0 +  assigned_treatment  + sms.treatment + distance + i(assigned_treatment, distance, ref = "control") | county,
            cluster = ~cluster.id
        )

    ols_sms_control_fit = sms_data %>%
        feols(
            data = .,
            fml = 
            dewormed ~ 0 + .[controls] + sms.treatment +  assigned_treatment + distance + i(assigned_treatment, distance, ref = "control") | county,
            cluster = ~cluster.id
        )

    ols_control_fit = data %>%
        feols(
            data = .,
            fml = 
            dewormed ~ 0 + .[controls] + distance + assigned_treatment +  i(assigned_treatment, distance, ref = "control") | county,
            cluster = ~cluster.id
        )

    tidy_ols_fit = ols_fit %>%
        tidy() %>%
        mutate(p.value = round(p.value, 4))
    tidy_ols_control_fit = ols_control_fit %>%
        tidy() %>%
        mutate(p.value = round(p.value, 4))
    tidy_ols_sms_fit = ols_sms_fit %>%
        tidy() %>%
        mutate(p.value = round(p.value, 4))
    tidy_ols_sms_control_fit = ols_sms_control_fit %>%
        tidy() %>%
        mutate(p.value = round(p.value, 4))

    return(lst(
        tidy_ols_fit,
        tidy_ols_control_fit,
        tidy_ols_sms_fit,
        tidy_ols_sms_control_fit
    ))
}

fits = generate_fits(analysis_data, analysis_sms_data, controls, cluster.dist.to.pot)
tidy_ols_fit = fits$tidy_ols_fit
tidy_ols_control_fit = fits$tidy_ols_control_fit


tidy_ols_control_fit %>%
    filter(str_detect(term, "treat|dist"))

tidy_ols_sms_fit = fits$tidy_ols_sms_fit
tidy_ols_sms_control_fit = fits$tidy_ols_sms_control_fit

tidy_ols_sms_control_fit %>%
    filter(str_detect(term, "treat|dist"))    



hh_fits = generate_fits(analysis_data, analysis_sms_data, controls, dist.to.pot)

hh_tidy_ols_fit = hh_fits$tidy_ols_fit
hh_tidy_ols_control_fit = hh_fits$tidy_ols_control_fit


hh_tidy_ols_sms_fit = hh_fits$tidy_ols_sms_fit
hh_tidy_ols_sms_control_fit = hh_fits$tidy_ols_sms_control_fit

hh_tidy_ols_sms_control_fit %>%
    filter(str_detect(term, "treat|dist"))    

analysis_data %>%
    colnames()



