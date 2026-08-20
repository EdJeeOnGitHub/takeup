library(tidyverse)

# ---- Load census data (contains monitored / true.monitored) ----
census_data_env = new.env()
load(file.path("data", "takeup_census.RData"), envir = census_data_env)
census_data = census_data_env$census.data %>%
  rename(census.consent = consent)

# ---- 1. Load the dict and raw wave 1 survey data ----
takeup_dict = read_csv("data/takeup_survey_dict_wave1.csv") %>%
  transmute(person_key, KEY.individ, cluster_key)

# Read each file separately (as the field notebook does) because:
# - Wave 1 has numeric person_key (needs dict lookup to get KEY.individ)
# - Kakamega has string person_key (already IS the KEY.individ)
raw_wave1 = read_csv("data/raw-data/Wave 1 Point of Treatment Form-ind.csv")
raw_kakamega = read_csv("data/raw-data/Kakamega Point of Treatment Form-ind.csv")

# ---- 2. Wave 1: How many surveyed individuals don't match the dict? ----
# Only wave 1 uses the dict (numeric person_key)
# Filter out empty/garbage rows: NA cluster_key, NA/0 person_key
valid_wave1 = raw_wave1 %>%
  filter(!is.na(cluster_key), !is.na(person_key), person_key != 0)

unmatched_wave1 = valid_wave1 %>%
  anti_join(takeup_dict, by = c("cluster_key", "person_key"))

cat("=== Wave 1 (numeric person_key, needs dict) ===\n")
cat("Total wave 1 survey rows:", nrow(raw_wave1), "\n")
cat("Valid rows (non-NA cluster_key & person_key > 0):", nrow(valid_wave1), "\n")
cat("Valid but not in takeup.dict:", nrow(unmatched_wave1), "\n")
cat("Unique clusters with unmatched:", n_distinct(unmatched_wave1$cluster_key), "\n")

# Breakdown of excluded rows
cat("\nExcluded row breakdown:\n")
raw_wave1 %>%
  mutate(
    reason = case_when(
      is.na(cluster_key) ~ "NA cluster_key",
      is.na(person_key) ~ "NA person_key",
      person_key == 0 ~ "person_key = 0",
      TRUE ~ "valid"
    )
  ) %>%
  count(reason) %>%
  print()

# ---- 3. Kakamega: person_key is already KEY.individ (string) ----
cat("\n=== Kakamega (string person_key = KEY.individ) ===\n")
cat("Total Kakamega survey rows:", nrow(raw_kakamega), "\n")
cat("person_key type:", class(raw_kakamega$person_key), "\n")

# Check if any Kakamega person_keys don't match census
kakamega_unmatched = raw_kakamega %>%
  filter(person_key != "0", !is.na(person_key)) %>%
  anti_join(census_data, by = c("person_key" = "KEY.individ"))

cat("Kakamega surveyed but not in census:", nrow(kakamega_unmatched), "\n")

# ---- 4. How many individuals are monitored but not true.monitored? ----
mon_not_truemon = census_data %>%
  filter(monitored & !true.monitored)

cat("\n=== Monitored vs True.Monitored ===\n")
cat("monitored=TRUE, true.monitored=FALSE:", nrow(mon_not_truemon), "\n")

census_data %>%
  count(wave, monitored, true.monitored) %>%
  print()

# ---- 5. Do any unmatched wave 1 entries fall in clusters with mon!=truemon? ----
mon_not_truemon_clusters = mon_not_truemon %>%
  distinct(cluster.id)

unmatched_in_affected_clusters = unmatched_wave1 %>%
  filter(cluster_key %in% mon_not_truemon_clusters$cluster.id)

cat("\nUnmatched wave 1 rows in clusters with mon!=truemon:", nrow(unmatched_in_affected_clusters), "\n")

# ---- 6. Inspect unmatched wave 1 entries ----
if (nrow(unmatched_wave1) > 0) {
  cat("\nUnmatched wave 1 survey entries:\n")
  unmatched_wave1 %>%
    select(cluster_key, person_key, any_of(c("person", "dewormed", "deworming_day"))) %>%
    print(n = 50)
}

# ---- 7. Inspect Kakamega unmatched entries ----
if (nrow(kakamega_unmatched) > 0) {
  cat("\nUnmatched Kakamega survey entries:\n")
  kakamega_unmatched %>%
    select(any_of(c("cluster_key", "person_key", "person", "dewormed", "deworming_day"))) %>%
    print(n = 50)
}

# ---- 8. Did name-matched individuals end up in analysis.data? ----
load(file.path("data", "analysis.RData"))

cat("\n=== Name-matched individuals in analysis.data ===\n")
analysis.data %>%
  filter(sms.treatment.2 == "sms.control") %>%
  count(name_matched, monitored) %>%
  print()

cat("\nDeworming status of name-matched individuals:\n")
analysis.data %>%
  filter(sms.treatment.2 == "sms.control") %>%
  group_by(name_matched) %>%
  summarise(
    n = n(),
    dewormed_any_true = sum(dewormed.any, na.rm = TRUE),
    dewormed_any_false = sum(!dewormed.any, na.rm = TRUE),
    dewormed_any_na = sum(is.na(dewormed.any)),
    dewormed_matched_true = sum(dewormed.matched, na.rm = TRUE),
    dewormed_matched_false = sum(!dewormed.matched, na.rm = TRUE),
    dewormed_matched_na = sum(is.na(dewormed.matched)),
    .groups = "drop"
  ) %>%
  print()

cat("\nWave breakdown of name-matched:\n")
analysis.data %>%
  filter(sms.treatment.2 == "sms.control") %>%
  count(wave, name_matched) %>%
  print()

    
