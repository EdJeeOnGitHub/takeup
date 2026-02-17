library(magrittr)
library(tidyverse)
library(lubridate)
library(haven)
library(broom)
library(scales)
library(sf)
library(here)


source("clean-analysis-util.R")
source("rct-design-fieldwork/takeup_rct_assign_clusters.R")
dir.create("data/clean-data", showWarnings = FALSE)

cluster_wave_county_data = read_rds(here("data", "takeup_cluster_wave_county_5.0.rds"))
cluster_strat_data = read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))
cluster_survey_data = read_rds(here("data", "takeup_cluster_survey.rds"))

# Load analysis data for deworming status in endline knowledge table
load(file.path("data", "analysis.RData"))

census.data
# Loading census data
# multiple datasets inside it with common df names
census_data_env = new.env()
with_env = function(f, e = parent.frame()) {
    stopifnot(is.function(f))
    environment(f) = e
    f
}
load_census_function = function(){
  load(file.path("data", "takeup_census.RData"))
  return(census.data)
}
census_data = with_env(load_census_function, census_data_env)() %>%
  rename(census.consent = consent) # Rename this to reduce chance of error


pot_info = read_rds(file.path("data", "pot_info.rds")) # Point of treatment information
# Geographic data for clusters in the study
pot_geo_info = pot_info %>% 
  filter(!is.na(wave)) %>% 
  transmute(cluster.id, 
            pot.lon = lon.post.rct.update,
            pot.lat = lat.post.rct.update) 

# Add distance to pot info to census data
census_data = census_data %>%
  left_join(pot_geo_info, by = "cluster.id") %>%
  group_by(cluster.id) %>%
  group_modify(~ {
    if (is.na(first(.y$cluster.id))) {
      individ.dist.to.pot <- NA
    } else {
      cluster.pot.sf <- st_as_sf(slice(.x, 1), coords = c("pot.lon", "pot.lat"), crs = wgs.84) %>%
        st_transform(kenya.proj4)

      individ.sf <- st_as_sf(.x, coords = c("lon", "lat"), crs = wgs.84) %>%
        st_transform(kenya.proj4)

      individ.dist.to.pot <- st_distance(individ.sf, cluster.pot.sf) %>%
        units::drop_units() %>%
        c()
    }

    mutate(.x, dist.to.pot = individ.dist.to.pot)
  }) %>%
  ungroup()

# Village location info
rct.villages <- read_rds(here("data", "rct_target_villages_2.0.rds")) %>%
  mutate(new.village = FALSE) %>%
  bind_rows(read_rds("data/rct_target_villages_2.0-4.rds") %>%
              mutate(new.village = TRUE))

all.villages <- rct.villages %>%
  bind_rows(anti_join(cluster_survey_data, ., c("cluster.id", "target.village.id"))) %>%
  mutate(targeted.village = !is.na(new.village),
         cluster.id = as.integer(cluster.id))

rm(rct.villages)

known.village.locations <- all.villages %>%
  filter(targeted.village, !is.na(target.lon), !is.na(target.lat)) %>%
  distinct(cluster.id, target.village.id, .keep_all = TRUE)


raw_endline_df = read_csv(raw.data.path("Endline Survey.csv"), guess_max = 10000) %>%
    mutate(
        across(c(SubmissionDate, starttime, endtime),
        ~ parse_datetime(., datetime.format, locale = locale(tz = "America/New_York")) %>%
            format(tz = "Africa/Nairobi") %>%
            parse_datetime(locale = locale(tz = "Africa/Nairobi"))
        )
    )


all_endline_data = raw_endline_df %>%
    filter(deviceid != "(web)", present == 1, SubmissionDate >= "2016-10-18") %>% 
    rename(cluster.id = cluster_id,
        KEY.individ = person,
        survey.type = survey_type,
        sms.treatment = sms,
        endline.type = endline_type,
        dewormed.reported = treated,
        any.sms.reported = receive_text,
        num.sms.reported = number_text,
        sms.content.reported = text_content, 
        lon = `gps-Longitude`,
        lat = `gps-Latitude`) %>%
  mutate(monitor.consent = coalesce(consent | reconsent, NA),
         sms.content.reported = map(sms.content.reported, ~ str_split(., " "))) %>%
  filter(!is.na(lon), !is.na(lat)) %>% 
  select(-county) %>% 
  left_join(cluster_wave_county_data, "cluster.id") %>% 
    # Data fabricated by enumerators 111
  filter(cluster.id != 1163 | SubmissionDate >= "2016-11-14", enumerator != 111) %>%
  filter(date(SubmissionDate) != "2016-11-7" | test == 1) %>% 
  left_join(identify_closest_cluster(.), "KEY")  %>%
  left_join(select(census_data, KEY.individ, lon, lat), "KEY.individ", suffix = c(".survey", ".census")) %>% 
  rename(lon = lon.survey, lat = lat.survey)  %>%
  mutate(
    dist.to.census = st_distance(st_as_sf(., coords = c("lon", "lat"), crs = wgs.84) %>% st_transform(kenya.proj4),
                                 st_as_sf(., coords = c("lon.census", "lat.census"), crs = wgs.84) %>% st_transform(kenya.proj4),
                                 by_element = TRUE) %>% as.numeric,
    far.from.census = dist.to.census > 50,
    sms_status = case_when(
        sms.treatment == "sms.control" ~ "control",
        sms.treatment %in% c("reminder.only", "social.info") ~ "treatment",
        TRUE ~ NA_character_
    )
  )


all_endline_data %>%
    filter(consent == 1) %>%
    group_by(sms_status) %>%
    summarise(
        n = n_distinct(KEY.individ)
    )

n_reach_df = all_endline_data %>%
    group_by(sms_status) %>%
    summarise(
        n = n_distinct(KEY.individ)
    ) %>%
    bind_rows(
        summarise(., n = sum(n)) %>%
            mutate(sms_status = "total")
    )


cat("Reached total: ", n_reach_df %>% filter(sms_status == "total") %>% pull(n), " individuals in endline survey\n")
cat("Reached SMS control: ", n_reach_df %>% filter(sms_status == "control") %>% pull(n), " individuals in endline survey\n")
cat("Reached SMS treatment: ", n_reach_df %>% filter(sms_status == "treatment") %>% pull(n), " individuals in endline survey\n")


endline_data = all_endline_data %>%
    filter(consent == 1)  %>%
    prepare_endline_data(census_data, cluster_strat_data)


n_complete_df = endline_data %>%
    group_by(sms_status) %>%
    summarise(
        n = n_distinct(KEY.individ)
    ) %>%
    bind_rows(
        summarise(., n = sum(n)) %>%
            mutate(sms_status = "total")
    )
cat("Completed total: ", n_complete_df %>% filter(sms_status == "total") %>% pull(n), " individuals in endline survey\n")
cat("Completed SMS control: ", n_complete_df %>% filter(sms_status == "control") %>% pull(n), " individuals in endline survey\n")
cat("Completed SMS treatment: ", n_complete_df %>% filter(sms_status == "treatment") %>% pull(n), " individuals in endline survey\n")





#### Create endline knowledge tables
raw_endline_know_table_data = bind_rows(table.A = read_csv(file.path("data", "raw-data", "Endline Survey-survey-sec_D-tableA.csv")) %>% 
                                       transmute(num.recognized = recogniseA,
                                                 dewormed = dewormedA,
                                                 second.order = order_2ndA,
                                                 second.order.reason = order_2nd_reason,
                                                 relationship = relationshipA,
                                                 relationship.other = relationshipA_other,
                                                 times.seen = times_seenA,
                                                 visited = visitedA, 
                                                 know.other.index = instanceA,
                                                 PARENT_KEY, KEY),
                                     table.B = read_csv(file.path("data", "raw-data", "Endline Survey-survey-sec_D-tableB.csv")) %>% 
                                       transmute(num.recognized = recogniseB,
                                                 dewormed = dewormedB,
                                                 dewormed.know.only = dewormedBB,
                                                 know.other.index = instanceB + 0:9, 
                                                 know.other.index.2 = know.other.index + 1, 
                                                 PARENT_KEY, KEY),
                                     .id = "know.table.type") %>% 
  mutate(recognized = num.recognized > 0)

raw_endline_know_table_data %>%
  summarise(
    n = n(),
    n_indiv = n_distinct(PARENT_KEY)
  )

survey_know_list = file.path("instruments", "SurveyCTO Forms", "Endline Survey", "Deployed Form Version", 
                              c("knowledge_list.csv", "kak_knowledge_list.csv")) %>% 
  map_dfr(read_csv, .id = "wave") %>% 
  mutate(wave = as.integer(wave))  %>%
  filter(wave == 1 | !is.na(survey.type))

endline_data %>%
    select(contains('treatment'), contains('dist'))

endline_know_table_data = raw_endline_know_table_data %>%
  # Join on actual individuals from the endline
  left_join(select(endline_data, KEY.individ, wave, cluster.id, assigned.treatment, sms.treatment, dist.pot.group, KEY), 
             by = c("PARENT_KEY" = "KEY")) %>% 
  # Use the knowledge list so we know who knows each other
  left_join(select(survey_know_list, person_key, know.other.index, KEY.individ.other),
             by = c("KEY.individ" = "person_key", "know.other.index")) %>% 
  left_join(select(survey_know_list, person_key, know.other.index, KEY.individ.other), # Get the second person (table B)
            by = c("KEY.individ" = "person_key", "know.other.index.2" = "know.other.index"),
            suffix = c(".1", ".2")) %>% 
  mutate(across(c(recognized, dewormed, dewormed.know.only), factor, levels = c(0:1, 98), labels = c("no", "yes", "don't know"))) %>% 
  mutate(across(second.order, factor, levels = c(1:2, 97:98), labels = c("yes", "no", "prefer not say", "don't know"))) %>% 
  mutate(across(relationship, factor, levels = c(1:5, 99), labels = c("hh member", "extended family", "friend", "neighbor", "church", "other"))) %>%
  mutate(relationship = if_else(fct_match(relationship, "other") & fct_match(str_to_lower(relationship.other), c("village member", "village mate", "same village", "village elder", "village mates", "villager")),
                                "village member", as.character(relationship)) %>% factor,
         dewormed = if_else(fct_match(know.table.type, "table.B") & is.na(dewormed), dewormed.know.only, dewormed)) %>% 
  left_join(filter(analysis.data, monitored) %>% transmute(KEY.individ, respondent.dewormed.any = dewormed.any),
            by = "KEY.individ") %>% 
  left_join(filter(analysis.data, monitored) %>% transmute(KEY.individ, actual.other.dewormed.any.1 = dewormed.any),
            by = c("KEY.individ.other.1" = "KEY.individ")) %>% 
  left_join(filter(analysis.data, monitored) %>% transmute(KEY.individ, actual.other.dewormed.any.2 = dewormed.any),
            by = c("KEY.individ.other.2" = "KEY.individ")) %>% 
  mutate(actual.other.dewormed.any.either = actual.other.dewormed.any.1 | (fct_match(know.table.type, "table.B") & actual.other.dewormed.any.2)) 


# Missing some people in here - let's find from where
endline_know_table_data %>%
  filter(is.na(KEY.individ)) %>%
  summarise(
    n_parents = n_distinct(PARENT_KEY)
  )



# HHs in endline, not in know table
in_endline_not_know_table = endline_data %>% 
  select(KEY.individ, KEY) %>%
  anti_join(raw_endline_know_table_data, by = c("KEY" = "PARENT_KEY")) %>%
  pull(KEY) %>%
  unique()

in_know_table_not_endline = raw_endline_know_table_data %>% 
  select(PARENT_KEY) %>%
  anti_join(endline_data, by = c("PARENT_KEY" = "KEY")) %>%
  pull(PARENT_KEY) %>%
  unique()

raw_endline_know_table_data %>%
  select(PARENT_KEY) %>%
  anti_join(all.endline.data %>% select(KEY), by = c("PARENT_KEY" = "KEY")) %>%
  pull(PARENT_KEY) %>%
  unique()  %>%
  length()

n_in_know_table = raw_endline_know_table_data %>%
  summarise(n = n_distinct(PARENT_KEY)) %>%
  pull(n)
n_in_endline = endline_data %>%
  summarise(n = n_distinct(KEY)) %>%
  pull(n)

n_in_know_table
n_in_endline
in_endline_not_know_table %>% length()
in_know_table_not_endline %>% length()


endline_know_table_data %>%
  summarise(
    n_parent_keys = n_distinct(PARENT_KEY),
    n_indiv_keys = n_distinct(KEY.individ),
    n_indiv_keys_know_treat = n_distinct(KEY.individ[!is.na(assigned.treatment)]),
    n_indiv_keys_know_dewormed = n_distinct(KEY.individ[!is.na(dewormed)])
  )

endline_data = endline_data %>%
  mutate(
    in_know_table = KEY %in% raw_endline_know_table_data$PARENT_KEY
  )


summ_endline_know_table = endline_know_table_data %>%
  filter(!is.na(KEY.individ)) %>%
  mutate(
    # convert table A first person (as only 1 in table A) to yes/no
    actual_other_dewormed_1_lgl = case_when(
      actual.other.dewormed.any.1 == TRUE ~ "yes",
      actual.other.dewormed.any.1 == FALSE ~ "no",
      TRUE ~ NA_character_
    ),
    # if "don't know", set to NA, otherwise check if classification matches actual deworming status of other person
    correct_classification_yesno = case_when(
      dewormed == "don't know" ~ NA,
      dewormed != "don't know" ~ (actual_other_dewormed_1_lgl == dewormed),
      TRUE ~ NA
    ),
    # if "don't know", set to FALSE, otherwise check if classification matches actual deworming status of other person
    correct_classification_yesnodk = case_when(
      dewormed == "don't know" ~ FALSE,
      dewormed != "don't know" ~ (actual_other_dewormed_1_lgl == dewormed),
      TRUE ~ NA
    )
    ) %>%
  group_by(KEY.individ, know.table.type) %>%
  summarise(
    obs_know_person = sum(num.recognized),
    obs_know_person_prop = mean(num.recognized),
    knows_other_dewormed = sum(fct_match(dewormed, c("yes", "no")), na.rm = TRUE),
    knows_other_dewormed_yes = sum(fct_match(dewormed, "yes"), na.rm = TRUE),
    knows_other_dewormed_no = sum(fct_match(dewormed, "no"), na.rm = TRUE),
    thinks_other_knows = sum(fct_match(second.order, c("yes", "no")), na.rm = TRUE),
    thinks_other_knows_yes = sum(fct_match(second.order, "yes"), na.rm = TRUE),
    thinks_other_knows_no = sum(fct_match(second.order, "no"), na.rm = TRUE),
    pct_correct_classification_yesno = mean(correct_classification_yesno, na.rm = TRUE),
    pct_correct_classification_yesnodk = mean(correct_classification_yesnodk, na.rm = TRUE)
  ) %>%
  ungroup()


summ_know_A_df = summ_endline_know_table %>%
  filter(know.table.type == "table.A") 

summ_know_B_df = summ_endline_know_table %>%
  filter(know.table.type == "table.B")

summ_endline_know_table %>%
  summarise(
    n_indiv = n_distinct(KEY.individ)
  )

endline_know_table_data %>%
  summarise(
    n_indiv = n_distinct(KEY.individ)
  )


# Save clean endline data
endline_data %>%
    write_csv("data/clean-data/clean-endline-data.csv")

summ_endline_know_table %>%
    write_csv("data/clean-data/clean-endline-know-table-data.csv")


endline_know_table_data %>%
    write_csv("data/clean-data/clean-endline-know-table-data-long.csv")

endline_data %>%
    write_rds("data/clean-data/clean-endline-data.rds")

summ_endline_know_table %>%
    write_rds("data/clean-data/clean-endline-know-table-data.rds")



endline_know_table_data %>%
    write_rds("data/clean-data/clean-endline-know-table-data-long.rds")
