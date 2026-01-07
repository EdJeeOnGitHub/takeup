library(magrittr)
library(tidyverse)
library(broom)
library(ggrepel)
library(ggmap)
library(ggstance)
library(gridExtra)
library(cowplot)
library(rgeos)
library(sp)
library(knitr)
library(modelr)
library(car)
library(rstan)
library(latex2exp)
library(ggthemes)
library(cmdstanr)
library(posterior)
library(tidybayes)

library(econometr)
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))


source(file.path("multilvlr", "multilvlr_util.R"))

fit_version <- 66

# 66 ed fit
# 60 Karim fit
# 62 also Karim fit


model_fit_by = if_else(fit_version %in% c(60, 62), "Karim", "Ed")


models_we_want = c(
  "STRUCTURAL_LINEAR_U_SHOCKS"
)


quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

output_basepath = str_glue("temp-data/output_dist_fit{fit_version}")

canva_palette_vibrant <- "Primary colors with a vibrant twist"

theme_set(theme_minimal() +
            theme(legend.position = "bottom"))

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

rct.schools.data <- read_rds(file.path("data", "takeup_rct_schools.rds"))
rct.cluster.selection <- read_rds(file.path("data", "rct_cluster_selection_2.0.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))

load(file.path("data", "takeup_village_pot_dist.RData"))

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data


## Fit Loading

load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.RData")))

dist_fit_data %<>%
  left_join(tribble(
    ~ fit_type,        ~ model_color,
      "fit",           "black", 
      "prior-predict", "darkgrey",
  ), by = "fit_type")

delta <- function(v, ...) dnorm(v, ...) / ((pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE)))




dist_fit_data


rf_files = fs::dir_ls("data/stan_analysis_data/", regex = str_glue("dist_fit{fit_version}_REDUCED_FORM_NO_RESTRICT.*\\.csv$"))

rf_fit = as_cmdstan_fit(rf_files)



rf_fit

rf_draws = as_draws_df(rf_fit)


colnames(rf_draws)


cluster_cf_cutoff_draws = rf_fit %>%
    gather_draws(cluster_cf_cutoff[a,b,c])

    # est_takeup_level = list(cluster_cf_cutoff, FALSE) %>% # map_lgl(model_type, fct_match, "structural")) %>%  
    #   pmap(organize_by_treatment, mu_assigned_treatment, assigned_treatment, assigned_dist_group) %>% 
    #   map2(cluster_cf_cutoff,
    #        ~ filter(.y, assigned_dist_group_obs == assigned_dist_group) %>%
    #          organize_by_treatment(condition_on_dist = TRUE, mu_assigned_treatment, assigned_treatment) %>%
    #          bind_rows(.x, .)),