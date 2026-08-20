
library(tidyverse)

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


monitored_sms_data = analysis.data %>% 
  filter(mon_status == "monitored") %>% 
  filter(have_phone == "Yes") %>%
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()
monitored_nosms_data

monitored_sms_data %>%
    count(sms.treatment)

monitored_data = analysis.data %>%
    filter(mon_status == "monitored") 

monitored_data %>%
    select(contains("consent")) %>%
    skimr::skim()

    summarise(
        n_distinct(KEY.individ)
    )
stop()


nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


analysis_data <- monitored_nosms_data %>% 
  ungroup() %>%
  mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group,
    standard_hh.dist.to.pot = standardize(dist.to.pot),
    indiv_id = row_number()
    )


cluster_dispersion_df = analysis_data %>%
  group_by(
    assigned_treatment,
    cluster_id
  ) %>%
  summarise(
    mse_dist_to_cluster = mean(((dist.to.pot - cluster.dist.to.pot)/1000)^2)
  ) %>%
  mutate(
    dispersed_community = mse_dist_to_cluster > 0.5
  ) %>%
  ungroup()
analysis_data = analysis_data %>%
  left_join(cluster_dispersion_df %>% select(-assigned_treatment), by = "cluster_id")  

if (str_detect(script_options$models, "INDIV_DIST_INDIV_FP")) {
  analysis_data = analysis_data %>%
    mutate(
      old_cluster_id = cluster_id,
      cluster.id = indiv_id, 
      cluster_id = indiv_id,
      standard_cluster.dist.to.pot = standard_hh.dist.to.pot
      )
}

if (str_detect(script_options$models, "NO_OUTLIERS")) {
  # recalculate cluster id
  analysis_data = analysis_data %>%
    filter(!dispersed_community) %>%
    group_by(cluster.id) %>% 
    mutate(cluster_id = cur_group_id()) %>% 
    ungroup()
}

# Models ------------------------------------------------------------------


total = 9935
baseline = 2250
endline = 3750
sms = 3935



baseline + endline



baseline + endline + sms


baseline.data

endline.know.table.data %>%
    filter(fct_match(know.table.type, "table.A"))  %>%
    group_by(KEY.individ) %>%
    mutate(num.recognized = sum(num.recognized))  %>%
    filter(num.recognized > 0) %>%
    ungroup() %>%
    summarise(
      n_indiv = n_distinct(KEY.individ)
    )


know_a_table = endline.know.table.data %>%
  filter(fct_match(know.table.type, "table.A")) 

sms_treat = monitored_sms_data %>%
  filter(sms.treatment %in% c("social.info", "reminder.only")) 

# out of total endline sample, we randomly sampled 1627 people to implement beliefs.
# Amongst them, XYZ didn't recognize anyone, and remaining ABC in the SMS sample, 
# who are not used in the main analysis.o


# out of endline sample, we randomly sampled 1627 people to implement beliefs.
# In the non-SMS sample this leads to ABC total beliefs observations, of which 
# XYZ didn't recognize anyone.


know_a_clusters = know_a_table %>%
  pull(cluster.id)

ana_clusters = analysis_data %>%
  pull(cluster.id)

setdiff(ana_clusters, know_a_clusters)


know_a_ids = know_a_table %>%
  pull(KEY.individ)

ana_ids = analysis_data %>%
  pull(KEY.individ)

setdiff( know_a_ids, ana_ids)

not_in_monitored = setdiff(know_a_ids, monitored_data$KEY.individ)

stop()



know_a_table %>%
  filter(KEY.individ %in% not_in_monitored)  %>%
  colnames()



know_a_table = know_a_table %>%
  group_by(KEY.individ) %>%
  mutate(num.recognized = sum(num.recognized)) %>%
  mutate(any_recognized = num.recognized > 0) %>%
  ungroup() %>%
  mutate(
    id_in_analysis = if_else(KEY.individ %in% analysis_data$KEY.individ, 1, 0),
    id_in_sms = if_else(KEY.individ %in% sms_treat$KEY.individ, 1, 0),
    id_in_either = if_else(KEY.individ %in% analysis_data$KEY.individ | KEY.individ %in% sms_treat$KEY.individ, 1, 0),
    id_in_mon = if_else(KEY.individ %in% monitored_data$KEY.individ, 1, 0)
  )

know_a_table %>%
  summarise(
    pr_in_ana = mean(id_in_analysis, na.rm = TRUE),
    pr_in_sms = mean(id_in_sms, na.rm = TRUE),
    pr_in_either = mean(id_in_either, na.rm = TRUE)
  )

know_a_table %>%
  group_by(id_in_analysis, id_in_sms, id_in_mon) %>%
  summarise(
    n_indiv = n_distinct(KEY.individ)
  )


endline.data %>%
  filter(
    KEY.individ %in% in_know_not_in_mon
  ) 
endline.data %>%
  filter(
    KEY.individ %in% in_know_not_in_mon
  ) %>%
  group_by(assigned.treatment) %>%
  summarise(n = n()) %>%
  count()

all.endline.data %>%
  filter(
    KEY.individ %in% in_know_not_in_mon
  )

in_know_not_in_mon = know_a_table %>%
  filter(id_in_mon == 0) %>%
  pull(KEY.individ)




all.endline.data %>%
  filter(KEY.individ %in% in_know_not_in_mon)  %>%
  skimr::skim()


all.endline.data %>%
  filter(!(KEY.individ %in% in_know_not_in_mon))  %>%
  skimr::skim()


all.endline.data %>%
  filter(!(KEY.individ %in% in_know_not_in_mon))  %>%
  count(return)


know_a_table %>%
  filter(id_in_mon == 0)  %>%
  count(dewormed)


know_a_table %>%
  filter(id_in_mon == 1)  %>%
  count(dewormed)

know_a_table %>%
  filter(id_in_analysis == 1) %>%
  group_by(KEY.individ) %>%
  summarise(
    n_recognized = sum(num.recognized)
  ) %>%
  count(n_recognized) 

know_a_table %>%
  group_by(id_in_analysis, id_in_sms, any_recognized) %>%
  summarise(
    n_indiv = n_distinct(KEY.individ)
  )

analysis_data %>%
    count(obs_know_person > 0)
