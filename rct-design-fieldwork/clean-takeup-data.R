library(tidyverse)




# Village data
load(file.path("data", "takeup_village_pot_dist.RData"))
# Analysis data
load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

# Only take sms.control HHs
nosms_data = analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


# Only take monitored, no sms HHs
monitored_nosms_data = analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  mutate(standard_dist.to.pot = standardize(dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


nosms_data %>%
  write_csv("data/clean-data/nosms-takeup-data.csv")

monitored_nosms_data %>%
  write_csv("data/clean-data/monitored-nosms-takeup-data.csv")


nosms_data %>%
  write_rds("data/clean-data/nosms-takeup-data.rds")

monitored_nosms_data %>%
  write_rds("data/clean-data/monitored-nosms-takeup-data.rds")

analysis.data %>%
  write_rds("data/clean-data/full-takeup-data.rds")