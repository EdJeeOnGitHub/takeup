
library(tidyverse)
library(fixest)

load(file.path("scratch/cal-pref/analysis.RData"))


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

pref_df = analysis_data %>%
  filter(!is.na(gift_choice), monitored, monitor.consent, !hh.baseline.sample.pool, !is.na(sms.treatment)) %>% 
  filter(
    assigned.treatment %in% c("control", "ink")
  ) %>%
  select(
    gift_choice, 
    assigned_treatment = assigned.treatment, 
    dist_group = dist.pot.group,
    cluster.id,
    county,
    dewormed,
    standard_cluster.dist.to.pot
    ) %>%
    mutate(pref_cal = gift_choice == "calendar") 



pref_fit = pref_df %>%
    filter(assigned_treatment == "control") %>%
    feols(
      pref_cal ~ dist_group | county,
      cluster = ~cluster.id
    )

pref_df %>%
  filter(assigned_treatment == "control" & dist_group == "close") %>%
  summarise(
    mean_pref_cal = mean(pref_cal)
  )

pref_fit %>%
    etable(
        tex = TRUE,
        depvar = TRUE,
        # title = "Preference for Calendar Gift Across Distance",
        style.tex = fixest::style.tex("aer"),
        file = "archive/code/scratch/cal-pref/pref-cal-gift.tex",
        replace = TRUE
    )

colnames(analysis_data)  = str_replace_all(colnames(analysis_data), "\\.", "_")
colnames(analysis_data)  = str_replace_all(colnames(analysis_data), "-", "_")

analysis_data %>%
  haven::write_dta("scratch/cal-pref/analysis-data.dta")

colnames(pref_df)  = str_replace_all(colnames(pref_df), "\\.", "_")
colnames(pref_df)  = str_replace_all(colnames(pref_df), "-", "_")

pref_df %>%
  haven::write_dta("scratch/cal-pref/pref-cal-data.dta")



