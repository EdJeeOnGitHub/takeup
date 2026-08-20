library(tidyverse)
library(here)

library(magrittr)
library(tidyverse)
library(lubridate)
library(haven)
library(broom)
library(scales)
library(knitr)
library(sf)

# Setup
source("rct-design-fieldwork/takeup_rct_assign_clusters.R")

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

datetime.format <- "%b %d, %Y %I:%M:%S %p"
takeup.datetime.type <- col_datetime(datetime.format)
takeup.date.type <- col_date(datetime.format)

raw.data.path <- . %>% here("data", "raw-data", .)

# Cluster Survey

## Point-of-Treatment

clusters.to.drop <- c(277, # Too close other cluster (I think 503)
                      491, 492, # Problematic urban clusters
                      1, # Village dispute about PoT
                      678, # Hostile community member
                      737) # Data fabrication and medication theft

cluster.wave.county.data <- read_rds(here("data", "takeup_cluster_wave_county_5.0.rds"))

cluster.info.to.drop <- c("uuid:3a491628-a0e8-4be0-9afc-c2cca08fc450", # Bad entries
                          "uuid:fcd3137e-7ddb-4f2e-9140-680e7638e42c",
                          "uuid:00fe16c6-5de1-4ff8-a9c8-92101ef4bd01",
                          "uuid:e3a47c54-f7f2-4af0-bd5b-0de12d997cb0")

cluster.info <- read_csv(raw.data.path("Cluster Survey V3 July 04.csv"),
                         col_types = list(SubmissionDate = takeup.date.type)) %>% 
  bind_rows(v3 = ., v1 = read_dta(here("data", "Cluster Survey V1.dta")) %>%
                   transmute(cluster.id = clusterid,
                             SubmissionDate = submissiondate,
                             KEY = key,
                             `gps2-Longitude` = gps2longitude,
                             `gps2-Latitude` = gps2latitude,
                             SubmissionDate = parse_date(SubmissionDate, "%d/%m/%Y %T"),
                             deviceid = as.character(deviceid),
                             manual_long2, manual_lat2, location_type, alt_name, comments, enumerator), 
            .id = "data.source") %>%
  filter(!is.na(cluster.id), deviceid != "(web)", !KEY %in% cluster.info.to.drop) %>% 
  arrange(cluster.id, SubmissionDate) %>% 
  group_by(cluster.id) %>%
  slice(n()) %>%
  ungroup %>%
  rename(alt.pot.lon = `gps2-Longitude`,
         alt.pot.lat = `gps2-Latitude`) %>%
  mutate(alt.pot.lon = ifelse(is.na(alt.pot.lon), read_lon_ant, alt.pot.lon),
         alt.pot.lat = ifelse(is.na(alt.pot.lat), read_lat_ant, alt.pot.lat),
         alt.pot.lon = ifelse(is.na(alt.pot.lon), manual_long2, alt.pot.lon),
         alt.pot.lat = ifelse(is.na(alt.pot.lat), manual_lat2, alt.pot.lat),
         location_type = factor(location_type, levels = 1:5, labels = c("Clinic", "Church", "Market", "Home", "Other")),
         cluster.id = as.integer(cluster.id)) %>% 
  rename(old.county.code = county) %>%
  mutate(old.county.code = factor(old.county.code, levels = 1:3, labels = c("Busia", "Kakamega", "Siaya"))) %>% 
  left_join(cluster.wave.county.data, "cluster.id") 

pot.info <- cluster.info %>%
  select(enumerator, enumerator_other, wave, county, cluster.id, alt.pot.lon, alt.pot.lat, alt_name, location_type, comments, SubmissionDate, data.source)

pot.verify.data <- read_csv(raw.data.path("POT verification.csv"), 
                            col_types = list(SubmissionDate = takeup.datetime.type)) %>%
  rename(cluster.id = cluster_id,
         lon.verify = `gps-Longitude`,
         lat.verify = `gps-Latitude`) %>% 
  select(-county) %>% 
  left_join(cluster.wave.county.data, "cluster.id") %>% 
  filter(wave == 1 | SubmissionDate >= "2016-10-16", !is.na(lon.verify), !is.na(lat.verify)) %>% 
  group_by(cluster.id) %>% 
  mutate(num.entries = n()) %>% 
  filter(min_rank(SubmissionDate) == n()) %>% 
  ungroup 

pot.info %<>% 
  left_join(#filter(pot.verify.data, !is.na(lon.verify), !is.na(lat.verify)),
            pot.verify.data,
            c("wave", "county", "cluster.id"), 
            suffix = c(".original", ".verify")) %>% 
  set_names(str_replace(names(.), "\\.original$", "")) 


## Villages

# `village.cluster.id` is not be trusted. Geographic analysis (below) of the closest cluster shows that `cluster.id` is more reliable. 

rct.village.codes <- read_csv(here("data", "village_codes_2.csv"), skip = 1, col_names = c("village.cluster.id", "village_name", "village")) %>% 
  mutate(village.cluster.id = as.integer(village.cluster.id))

# We're already reading this above (cluster-survey-data), but not in exactly the same format...
cluster.survey.data <- read_rds(here("data", "takeup_cluster_survey.rds"))

rct.villages <- read_rds(here("data", "rct_target_villages_2.0.rds")) %>% 
  mutate(new.village = FALSE) %>% 
  bind_rows(read_rds("data/rct_target_villages_2.0-4.rds") %>% 
              mutate(new.village = TRUE))
  
all.villages <- rct.villages %>% 
  bind_rows(anti_join(cluster.survey.data, ., c("cluster.id", "target.village.id"))) %>% 
  mutate(targeted.village = !is.na(new.village),
         cluster.id = as.integer(cluster.id),
         target.village_name = str_trim(target.village_name) %>% str_replace_all("\\s+", " ")) %>% 
  left_join(rct.village.codes, c("cluster.id" = "village.cluster.id", "target.village_name" = "village_name")) %>% 
  left_join(select(pot.info, cluster.id, matches("alt.pot.(lon|lat)"), location_type), "cluster.id") %>% 
  group_by(is.na(village)) %>% 
  do({ # Adding some village IDs for village that weren't in the code file shared by Arthur
    if(is.na(first(.$village))) {
      mutate(., village = 1000 + seq_len(nrow(.)))
    } else {
      (.)
    }
  }) %>% 
  ungroup %>% 
  # EJ: No longer runs without rgeos which is deprecated
  # mutate(dist.group = convert.to.sp(., ~ target.lon + target.lat, wgs.84) %>% 
  #          spTransform(kenya.proj4) %>% 
  #          gDistance(byid = TRUE) %>% 
  #          as.dist %>% 
  #          hclust %>% 
  #          cutree(h = 1000)) %>% 
  # group_by(cluster.id) %>% 
  # mutate(num.dist.groups = unique(dist.group) %>% length) %>% 
  ungroup %>% 
  left_join(cluster.wave.county.data, "cluster.id") %>% 
  mutate(village.name.group = target.village_name %>% 
           str_trim %>% 
           str_replace("\\s*[A-Z]\\d?$", "") %>% 
           str_replace_all("\\s+", " ") %>% 
           str_replace_all("'", "") %>% 
           str_to_upper %>% 
           str_replace("\\s+[A-Z]\\d?$", "") %>% 
           str_replace(regex("\\s+(village|upper|lower|group|east|west|north|south|rural|urban|township|central)$", ignore_case = TRUE), "") %>% 
           str_replace(regex("\\s+(MWILUECHINA)$", ignore_case = TRUE), "") %>% 
           str_replace(regex("\\s*(\\(.+\\)|estate|township)$", ignore_case = TRUE), "") %>% 
           str_replace("\\s+[A-Z]\\d?$", ""),
         vill.name.dist.group = adist(village.name.group) %>% 
           as.dist %>% 
           hclust %>% 
           cutree(h = 1))

rm(rct.villages)

rct.village.codes %<>%
  left_join(select(all.villages, village, targeted.village), "village")

pre.census.data <- read_rds("data/pre.census.processed.rds")

all.villages %>% count(vill.name.dist.group) %>% filter(n > 1) %>% left_join(all.villages, "vill.name.dist.group") %>% select(village.name.group, vill.name.dist.group) %>% print(n = 300)

# Census

## Data

rct.counties <- c("Busia", "Siaya", "Kakamega")
busia.subcounties <- c("butula", "nambale", "teso south", "teso north") 
siaya.subcounties <- c("gem", "ugenya", "ugunja")

ke.lvl2.adm.data <- read_rds(here("data", "adm", "KEN_adm2.rds"))

counties.adm.data <- ke.lvl2.adm.data[ke.lvl2.adm.data$NAME_1 %in% rct.counties, ] #, "Vihiga"), ]
# subcounties.adm.data <- counties.adm.data[!counties.adm.data$NAME_1 %in% c("Busia", "Siaya") | counties.adm.data$NAME_2 %in% str_to_title(c(busia.subcounties, siaya.subcounties)), ]

county.bbox <- counties.adm.data@bbox

day1.wave1 <- as_date("2016-10-03")
day12.wave1 <- day1.wave1 + days(11)
day1.wave2 <- as_date("2016-10-24")
day12.wave2 <- day1.wave2 + days(11)

validate.coords <- . %>% 
  mutate(invalid.coord = 
           (!is.na(lon) & (lon > county.bbox["x", "max"] | lon < county.bbox["x", "min"])) |
           (!is.na(lat) & (lat > county.bbox["y", "max"] | lat < county.bbox["y", "min"])))

tu.data.reader <- function(file.name, submit.datetime.type = NULL, .other.types = NULL) { # =  "%b %d, %Y %I:%M:%S %p") {
  col.types <- list(SubmissionDate = if (is.null(submit.datetime.type)) takeup.datetime.type else submit.datetime.type,
                                       manual_long = col_number(),
                                       manual_lat = col_number()) %>% 
    c(.other.types)
  
  read_csv(file.name, col_types = col.types, guess_max = 100000) %>% 
    mutate(isValidated = isValidated == "true") %>% 
    rename(lat = `gps-Latitude`,
           lon = `gps-Longitude`,
           cluster.id = cluster_id) %>% 
    filter(deviceid != "(web)") %>% 
    mutate(lon = ifelse(is.na(lon), manual_long, lon),
           lat = ifelse(is.na(lat), manual_lat, lat)) %>% 
    validate.coords()
}

hh.census.data.pre24th <- tu.data.reader(here("data", "census_gps_util23.csv"), 
                                         submit.datetime.type = col_datetime("%d/%m/%Y %H:%M")) %>% 
  filter(SubmissionDate >= "2016-08-22") %>% 
  mutate_at(vars(present, return, number_individuals, instance, check), as.numeric)

factor.have.phone <- . %>% factor(levels = 0:2, labels = c("No", "Yes", "Don't know number"))

raw_census_df = read_csv(raw.data.path("Census.csv"), guess_max = 100000)

raw_census_df %>%
    count(cluster_id)

hh.census.data <- tu.data.reader(raw.data.path("Census.csv"), 
                                 .other.types = list(`gps-Latitude` = col_double(),
                                                     `gps-Longitude` = col_double(),
                                                     `gps-Altitude` = col_double(),
                                                     `gps-Accuracy` = col_double(),
                                                     deviceid = col_character(),
                                                     gps_work = col_integer())) %>% 
  filter(SubmissionDate >= "2016-08-25") %>%
  select(!c(starts_with("hhh_two_digits"))) |> 
  bind_rows(select(hh.census.data.pre24th, !c(hhh_name1st_check, directions, starts_with("hhh_two_digits")))) %>% 
  select(-c(RI, matches("^IC\\d{1,2}$"))) %>% 
  left_join(rct.village.codes, by = "village") %>%
  mutate(village_name = factor(village_name),
         present = ifelse(is.na(present) & SubmissionDate <= "2016-08-22", NA, present),
         return = as.numeric(return),
         hhh_have_phone = factor.have.phone(hhh_have_phone))


missing.gps.data <- read_csv(here("data", "missing_gps.csv"), guess_max = 1000000) %>% 
  select(KEY, `gps-Latitude`, `gps-Longitude`, gps_work, manual_long, manual_lat) %>% 
  mutate(lon = ifelse(is.na(`gps-Longitude`), manual_long, `gps-Longitude`),
         lat = ifelse(is.na(`gps-Latitude`), manual_lat, `gps-Latitude`),
         data.source = "missing_gps",
         across(c("lon", "lat"), as.double)) %>% 
  select(-c(starts_with("gps-"), starts_with("manual_"))) %>% 
  filter(!is.na(KEY), !is.na(lon), !is.na(lat))  

missing.gps.149.data <- read_csv(here("data", "missing_gps_149.csv")) %>% 
  transmute(lon = gpslongitude,
            lat = gpslatitude, 
            KEY = `key of form with missing GPS data`,
            data.source = "missing_gps_149")

missing.gps.360.data <- read_csv(here("data", "missing_gps_360.csv")) %>% 
  transmute(lon = gpslongitude,
            lat = gpslatitude, 
            KEY = `key of form missing GPS data`,
            data.source = "missing_gps_149")

missing.gps.844.data <- read_csv(here("data", "missing_gps_844.csv")) %>%
  rename(key.gps = `Key with GPS`,
         key.missing.gps = `Key missing GPS`) %>%
  inner_join(hh.census.data, c("key.gps" = "KEY")) %>% 
  transmute(lon, lat, KEY = key.missing.gps, data.source = "missing_gps_844")

hh.census.data %<>% {
  mask <- .$KEY %in% missing.gps.844.data$KEY

  mutate_at(., vars(lon, lat), funs(ifelse(mask, NA, .)))  
}

missing.gps.data %<>% 
  bind_rows(missing.gps.149.data, missing.gps.360.data, missing.gps.844.data) %>% 
  validate.coords()

hh.census.data %<>% 
  left_join(missing.gps.data, by = "KEY") %>% 
  filter(is.na(data.source) | data.source != "missing_gps_149" | cluster.id == 149) %>%
  mutate(lon = ifelse(is.na(lon.x), lon.y, lon.x),
         lat = ifelse(is.na(lat.x), lat.y, lat.x),
         gps_work = ifelse(is.na(gps_work.y), gps_work.x, as.numeric(gps_work.y)),
         invalid.coord = ifelse(is.na(invalid.coord.y), invalid.coord.x, invalid.coord.y)) %>% 
  select(-matches("\\.[xy]$"))


to_remove_data = map_dfr(paste0("data/", c("takeup_returns_toremove.csv", 
                    "takeup_moved_notvalid_toremove.csv",
                    "takeup_noconsent_toremove.csv")), read_csv) %>% 
  transmute(key = key)

hh.census.data %<>% 
  anti_join(to_remove_data,
            by = c("KEY" = "key"))

hh.id.dict <- read_rds("data/takeup_hh_id_dict.rds")

hh.census.data %<>% 
  select(-county) %>% 
  left_join(cluster.wave.county.data, "cluster.id") %>% 
  left_join(hh.id.dict, "KEY")


person.id.dict <- read_rds("data/takeup_person_id_dict.rds")

clean.names <- . %>% 
  mutate_at(vars(matches("1st|mid|2nd|clan")), str_to_upper) %>% 
  mutate_at(vars(matches("1st|mid|2nd|clan")), 
            funs(ifelse(. %in% c("NONE", "N/A", "NA", "NO"), NA, .))) %>% 
  mutate_at(vars(matches("1st|2nd|clan")), 
            funs(ifelse(nchar(.) < 2, NA, .))) %>% 
  mutate(last_name = ifelse(is.na(name2nd) & !is.na(name_mid), name_mid, name2nd),
         name_mid = ifelse(is.na(name2nd) & !is.na(name_mid), NA, name_mid)) 

census.data <- hh.census.data %>% {
  individ.data <- read_csv(raw.data.path("Census-survey-individual.csv"))
  anti_join(., individ.data, c("KEY" = "PARENT_KEY")) %>% 
    count(SubmissionDate) %>% 
    rename(num.hhs.without.individuals = n) %>% 
    knitr::kable(col.names = c("Submission Date", "Number HH w/o individuals")) %>% 
    print
  
  inner_join(., individ.data, c("KEY" = "PARENT_KEY"), suffix = c(".hh", ".individ"))
} %>% 
  mutate(two.digit.match = ifelse(have_phone == 1, two_digits == two_digits_check, NA),
         have_phone = factor.have.phone(have_phone)) %>% 
  group_by(KEY) %>% 
  mutate(num.individuals = n(),
         hh.has.phone = any(have_phone == "Yes"),
         hh.has.non.phone = any(have_phone != "Yes")) %>% 
  ungroup %>% 
  clean.names %>% 
  filter(!is.na(name1st) & (!is.na(name_mid) | !is.na(name2nd))) %>% 
  left_join(person.id.dict, "KEY.individ")

#This step would also remove all households with no individuals (inner join)
hh.census.data %<>%
  inner_join(distinct(census.data, KEY, num.individuals, hh.has.phone, hh.has.non.phone), "KEY")

# cluster.wave.county.data <- hh.census.data %>% 
#   filter(!is.na(cluster.id)) %>% 
#   distinct(cluster.id, wave, county)

cluster.wave.county.data <- cluster.info %>% 
  distinct(wave, county, cluster.id) %>% 
  mutate(wave = if_else(.$cluster.id %in% clusters.to.drop, 0, wave) %>% na_if(0))


load("data/takeup_cluster_review.RData")


hh.census.data %>% 
  group_by(SubmissionDate) %>% 
  summarize(num.hhs = n(),
            num.individuals = sum(num.individuals),
            num.clusters = unique(cluster.id) %>% { length(.) - any(is.na(.)) },
            num.villages = unique(village) %>% { length(.) - any(is.na(.)) }) %>% 
  ungroup %>% 
  knitr::kable(col.names = c("Submission Date", "Households", "Individuals", "Clusters", "Villages"))

# <!-- Testing the number of individuals in the data against reported number: -->

# This isn't very important; we are expecting some difference (see Arthur's response on Slack)
hh.census.data %>%
  select(number_individuals, num.individuals) %>%
  summarise(across(c(number_individuals, num.individuals), sum, na.rm = TRUE) )

hh.census.data %>% 
  validate::check_that(number_individuals == num.individuals) %>% 
  summary %>% 
  knitr::kable()



pot.ranges <- pot.info %>% 
  convert.to.sp(~ alt.pot.lon + alt.pot.lat, wgs.84) %>% 
  buffer.clusters(.width = 2500) %>% 
  spTransform(wgs.84) %>% {
    merge(tidy(.), 
          select(.@data, cluster.id), 
          by.x = "id", by.y = "row.names", all.x = TRUE, all.y = FALSE) 
  } %>% 
  mutate(boundary.type = "pot") 

get.survey.boundary <- function(.data, .boundary.type, .group.by = c("cluster.id")) {
  .data %>% 
    filter(!is.na(lon), !is.na(lat)) %>% 
    group_by_(.dots = .group.by) %>% 
    do((function(cluster.pc.hh) {
      cluster.pc.hh %>% 
        convert.to.sp(~ lon + lat, wgs.84) %>% 
        gConvexHull %>% 
        spTransform(kenya.proj4) %>% 
        gBuffer %>% 
        spTransform(wgs.84) %>% 
        tidy %>%
        mutate(group = paste(group, cluster.pc.hh$cluster.id[1], sep = "-")) 
    })(.)) %>% 
    ungroup %>% 
    mutate(boundary.type = .boundary.type)
}

pre.census.hhs.boundary <- pre.census.data %>% 
  mutate(cluster.id = as.integer(act.cluster.id)) %>% 
  get.survey.boundary("pre.census.convex.hull")

boundaries <- bind_rows(pot.ranges, pre.census.hhs.boundary) %>% 
  mutate(group = paste(group, boundary.type, sep = "-")) %>% 
  rename(lon = long)  

known.pot.locations <- pot.info %>%
  filter(!is.na(alt.pot.lon) | !is.na(lon.verify), !is.na(alt.pot.lat) | !is.na(lat.verify)) %>%
  mutate(alt.pot.lon = ifelse(is.na(alt.pot.lon), lon.verify, alt.pot.lon),
         alt.pot.lat = ifelse(is.na(alt.pot.lat), lat.verify, alt.pot.lat)) %>% 
  distinct(cluster.id, .keep_all = TRUE)

known.village.locations <- all.villages %>%  #rct.villages %>%
  filter(targeted.village, !is.na(target.lon), !is.na(target.lat)) %>%
  distinct(cluster.id, target.village.id, .keep_all = TRUE)

identify.closest.cluster <- function(.data, data.coords.formula = ~ lon + lat,  key.variable = "KEY") {
  # min.dist.df <- .data %>% 
  .data %>%
    filter(rowSums(is.na(model.frame(data.coords.formula, data = ., na.action = NULL))) == 0) %>% 
    # filter(!is.na(lon), !is.na(lat)) %>%
    (function (.located.data) {
      convert.to.sp(.located.data, data.coords.formula, wgs.84) %>% 
        spTransform(kenya.proj4) %>% 
        gDistance(known.village.locations %>%
                    convert.to.sp(~ target.lon + target.lat, wgs.84) %>%
                    spTransform(kenya.proj4),
                  byid = TRUE) %>%
        plyr::adply(2, function(dist.col) {
          in.range.dist <- dist.col[dist.col <= 2000]
      
          if (length(in.range.dist) == 0) {
            tibble(min.dist = NA, closest.cluster = NA)
          } else {
            tibble(min.dist = min(in.range.dist),
                   closest.village = known.village.locations$target.village.id[which(dist.col <= min.dist)], 
                   closest.cluster = known.village.locations$cluster.id[which(dist.col <= min.dist)]) %>% 
              distinct(closest.cluster, .keep_all = TRUE)
          }
          
        }, .parallel = TRUE) %>% 
        bind_cols(.located.data, .) %>% 
        select_(.dots = c(key.variable, "min.dist", "closest.cluster")) 
    })
}

hh.census.data %<>%
  filter(!is.na(lon), !is.na(lat), !invalid.coord, cluster.id %in% known.village.locations$cluster.id) %>%
  identify.closest.cluster %>% 
  right_join(hh.census.data, by = "KEY") %>% 
  mutate(misidentified.cluster = cluster.id != closest.cluster)

hh.census.data %>% 
  count(misidentified.cluster)

# Aside from some corrections of cluster ID due to incorrect coding, 

# * we're dropping households from cluster 954 that were found in the census that are very different from the pre-census

cluster.954.to.drop <- read_rds("data/cluster_954_to_drop.rds")

get.outside.keys <- function(.data, .cluster.id) 
  .data %>% 
    filter(cluster.id == .cluster.id, !invalid.coord, !is.na(lon), !is.na(lat)) %>% 
    convert.to.sp(~ lon + lat, wgs.84) %>% 
    magrittr::extract(!gWithin(., boundaries %>% 
                                filter(cluster.id == .cluster.id, boundary.type == "pre.census.convex.hull") %>% 
                                convert.to.sp(~ lon + lat, wgs.84) %>% 
                                gConvexHull %>% 
                                buffer.clusters(50) %>% 
                                spTransform(wgs.84), byid = TRUE) %>% 
                        c, ) %>% `@`(data) %$% KEY
  
corrected.hhs <- hh.census.data %>%
  filter(cluster.id %in% c(492, 254, 149, 360, 834, 587, 951, 1272, 735, 293, 1313), 
         misidentified.cluster) %>% 
  mutate(cluster.id = closest.cluster,
         misidentified.cluster = FALSE) %>% 
  select(-c(wave, county)) %>% 
  left_join(cluster.wave.county.data, "cluster.id")

cluster.hh.to.drop <- c(954, 197, 1442, 1242, 378, 431, 1426) %>% 
  map(~ get.outside.keys(hh.census.data, .)) %>% 
  c(c(147, 587, 735, 293, 844) %>% # Drop out of range households 
      map(~ filter(hh.census.data, cluster.id == ., is.na(closest.cluster)) %$% KEY)) %>% 
  unlist %>% 
  c("uuid:13ca71aa-d966-4956-8d8b-83673fd3f74a",
    "uuid:04d0d59b-e33b-49c6-b55f-ccf2eda14213",
    "uuid:fbedfc11-9e02-46c6-bf04-dfb817facefc",
    "uuid:611bf6c7-a84c-4210-809a-1d2b072dc682",
    "uuid:937590a7-b024-43c0-9bff-53d87defc663")


hh.census.data <- corrected.hhs %>% 
  bind_rows(anti_join(hh.census.data, ., "KEY")) %>% 
  mutate(cluster.id = ifelse(KEY %in% cluster.hh.to.drop, NA, cluster.id)) #%>% 
         # county = if_else(cluster.id == 103, "Siaya", as.character(county)) %>% factor,
         # wave = if_else(cluster.id == 103, 1, wave))

census.data %<>% { 
  bind_rows(anti_join(., select(corrected.hhs, KEY), "KEY"), 
            left_join(select(corrected.hhs, KEY, wave, county, cluster.id), select(., -c(wave, county, cluster.id)), "KEY"))
} %>% 
  mutate(cluster.id = ifelse(KEY %in% cluster.954.to.drop, NA, cluster.id))
         # county = if_else(cluster.id == 103, "Siaya", as.character(county)) %>% factor,
         # wave = if_else(cluster.id == 103, 1, wave))

### Households in Clusters

cluster.census.data <- hh.census.data %>% 
  group_by(cluster.id) %>% 
  summarize(num.hhs = n(), 
            num.invalid.gps.coord = sum(invalid.coord),
            num.hhs.gps.work = sum(gps_work, na.rm = TRUE),
            num.hhs.gps.work.na = sum(is.na(gps_work)),
            num.hhs.no.coord = sum(is.na(lon) | is.na(lat)),
            num.hhs.misidentified = sum(misidentified.cluster, na.rm = TRUE),
            num.hhs.out.of.range = sum(is.na(closest.cluster)),
            num.hhs.present = sum(present, na.rm = TRUE),
            num.hhs.present.na = sum(is.na(present)),
            num.hhs.return = sum(return, na.rm = TRUE),
            num.hhs.return.na = sum(is.na(return)),
            last.submission = max(SubmissionDate)) %>% 
  ungroup %>% 
  arrange(last.submission)

knitr::kable(cluster.census.data, 
             col.names = c("Cluster", "Total Count", "Invalid Coords", "GPS Worked", "GPS Worked NA", "No Coords", "Misidentified", "Out of Range", "Present", "Present NA", "Return", "Return NA", "Last Submission"))


# Remove all households that were not present, moved (_and_ there are no individuals record in that household), or did not consent.

hh.census.data %<>% filter(!invalid.coord, 
                           is.na(present) | (present == 1),
                           is.na(moved) | (moved == 0) | (num.individuals > 0),
                           is.na(consent) | (consent == 1),
                           !is.na(cluster.id),
                           !is.na(lon), !is.na(lat)) %>% 
  mutate(old.cluster.id = ifelse(cluster.id %in% clusters.to.drop, cluster.id, NA),
         cluster.id = ifelse(cluster.id %in% clusters.to.drop, NA, cluster.id)) %>% 
                           # !cluster.id %in% clusters.to.drop) %>% 
  group_by(cluster.id) %>% 
  filter(n() >= 50) %>% 
  ungroup

old.census.villages <- census.data %>% 
  select(KEY.individ, cluster.id, village)

census.data %<>% select(-cluster.id) %>% right_join(distinct(hh.census.data, cluster.id, old.cluster.id, KEY), "KEY") 

census.hhs.boundaries <- hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  get.survey.boundary("census.convex.hull") %>% 
  mutate(group = paste(group, boundary.type, sep = "-")) %>% 
  rename(lon = long)  

village.census.hh.boundaries <- hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  get.survey.boundary("census.convex.hull", .group.by = c("cluster.id", "village_name")) %>% 
  mutate(group = paste(group, boundary.type, sep = "-")) %>% 
  rename(lon = long)  

boundaries %<>% bind_rows(census.hhs.boundaries)  

# Below is a table comparing the number of households in the census to the pre-census:

num.hh.data <- hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  count(cluster.id) %>% 
  left_join(pre.census.data %>%
              mutate(act.cluster.id = as.integer(act.cluster.id)) %>% 
              count(act.cluster.id),
            by = c("cluster.id" = "act.cluster.id"),
            suffix = c(".census", ".pre.census"))

num.hh.data

# It appears that quite a few of the houses that were found during the pre-census are not actually standalone households. Some indviduals live in a separate structure, but eat their means under a different roof.

# Distribution of number of households:
# summary(cluster.census.data$num.hhs)



num.hh.data %>% 
  gather(survey, num.hhs, -cluster.id) %>% 
  ggplot() +
  geom_bar(aes(factor(cluster.id), y = num.hhs, fill = survey), stat = "identity", position = "dodge") +
  scale_fill_discrete("Survey", labels = c("Census", "Pre-census")) +
  ylab("Number of Households") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  count(cluster.id) %>% 
  ggplot() +
  geom_histogram(aes(n), binwidth = 50, boundary = 50) +
  scale_x_continuous("Individuals/Cluster", breaks = seq(0, 600, 50))

census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  distinct(cluster.id, KEY) %>% 
  count(cluster.id) %>% { 
    ggplot(.) +
      geom_histogram(aes(n), binwidth = 50, boundary = 50) +
      scale_x_continuous("Households/Cluster", breaks = seq(0, 600, 50)) 
  } 

cluster.census.data %>% 
  ggplot() +
  geom_histogram(aes(num.hhs), binwidth = 50) +
  scale_y_continuous(breaks = seq(0, 40, 5)) 


hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  ggplot() +
  geom_bar(aes(factor(num.individuals))) +
  xlab("Individuals/HH")  



### Number of Individuals in Clusters and Villages

census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  count(cluster.id)

census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  count(cluster.id, village) # village_name)

### Cluster Distribution by County

hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  distinct(cluster.id, county) %>% 
  ggplot() +
  geom_bar(aes(county), alpha = 0.5, color = "black") +
  labs(x = "", y = "")


## Cluster Household Maps

# Below is a separate map per cluster. Small filled circles indicate the point-of-treatments.  Triangles indicate targeted village centers. 

plot.census.hhs <- function(.data, .rct.villages.data, .cluster.id, .zoom = 13, 
                            .label = TRUE, 
                            .hh.types = c("all", "focal.only"), 
                            .hh.color = c("cluster.relation", "reported.village"),
                            use.google.maps = TRUE) {
  factor.cluster.type <- . %>% 
    factor(levels = c("focal", "not.focal", "unkown.cluster", "misidentified.cluster", "out.of.range"), 
           labels = c("Focal", "Not Focal", "Unknown Cluster", "Incorrect Cluster", "Out of Range"))
  
  .data %<>%
    mutate(cluster.type = case_when(is.na(.$cluster.id) ~ "unknown.cluster",
                                    is.na(.$closest.cluster) ~ "out.of.range",
                                    .$closest.cluster != .$cluster.id ~ "misidentified.cluster",
                                    .$cluster.id %in% .cluster.id ~ "focal",
                                    TRUE ~ "not.focal") %>% 
             factor.cluster.type)
  
  .rct.villages.data %<>% 
    mutate(cluster.type = ifelse(is.na(cluster.id), 
                                 "unknown.cluster",
                                 ifelse(cluster.id %in% .cluster.id,
                                        "focal",
                                        "not.focal")) %>% 
             factor.cluster.type)
  
  focal.hhs <- .data %>% 
    filter(cluster.id == .cluster.id, cluster.type == "Focal")
  
  cluster.boundaries <- boundaries %>%
    filter(cluster.id %in% .cluster.id)
  
  if (empty(focal.hhs)) {
    hhs.center <- .rct.villages.data[1, c("alt.pot.lon", "alt.pot.lat"), drop = TRUE] %>% 
      unlist
  } else {
    hhs.center <- focal.hhs %>% 
      select(lon, lat) %>% 
      bind_rows(cluster.boundaries %>% 
                  filter(boundary.type == "pot") %>% 
                  select(lon, lat)) %>% 
      make_bbox(lon, lat, data = .) 
    
  }
  
  hh.color.var <- switch(match.arg(.hh.color), 
                         cluster.relation = "cluster.type",
                         reported.village = "villaage")
  
  if (use.google.maps) {
    hhs.center %<>% 
      matrix(nrow = 2) %>% 
      rowMeans()
    
    plot.obj <- ggmap(get_googlemap(center = hhs.center,
                                    maptype = "terrain", # "hybrid", 
                                    zoom = .zoom, 
                                    scale = 2,
                                    style = "element:labels|visibility:off", # Drop all labels from map
                                    key = config$google_api_key)) 
  } else {
    plot.obj <- ggplot() +
      coord_fixed(xlim = hhs.center[c("left", "right")], ylim = hhs.center[c("bottom", "top")]) 
  }
  
  plot.obj <- plot.obj +
    geom_point(aes_string("lon", "lat", fill = hh.color.var), shape = 22, alpha = 0.25, size = 2, stroke = 0, 
               data = switch(match.arg(.hh.types), all = .data, focal.only = focal.hhs)) 
  
  if (.label) {
    plot.obj <- plot.obj +
      geom_label_repel(aes(alt.pot.lon, alt.pot.lat, label = cluster.id), segment.color = "black",  
                       data = .rct.villages.data %>% 
                         filter(cluster.type == "Not Focal") %>% 
                         distinct(cluster.id, .keep_all = TRUE)) 
  }
  
  
  plot.obj <- plot.obj +
    geom_point(aes(alt.pot.lon, alt.pot.lat, fill = cluster.type), shape = 21, color = "black", size = 2, stroke = 1, 
               data = distinct(.rct.villages.data, cluster.id, .keep_all = TRUE)) +
    geom_point(aes(target.lon, target.lat, fill = cluster.type), shape = 24, color = "black", size = 2, stroke = 1, 
               data = .rct.villages.data) +
    ggtitle(sprintf("Cluster %s", paste(.cluster.id, collapse = ", "))) +
    scale_color_discrete("Boundaries", labels = c("Census", "PoT", "Pre-census")) +  
    scale_fill_discrete(if (hh.color.var == "cluster.type") "Cluster Type" else "Village", 
                            drop = hh.color.var != "cluster.type") +
    labs(x = "", y = "") +
    theme(axis.text = element_blank())
  

  if (!is.null(cluster.boundaries)) {
    plot.obj <- plot.obj +
      geom_polygon(aes(lon, lat, group = group, color = boundary.type), linetype = "dotted", alpha = 0, data = cluster.boundaries)
  }
  
  return(plot.obj)
}

generate.hh.plots <- function(.data, cluster.ids = NULL) { 
  .data %>% 
    filter(!invalid.coord, !is.na(cluster.id)) %>% {
      if (!is.null(cluster.ids)) filter(., cluster.id %in% cluster.ids) else return(.)
    } %>% 
    filter(!is.na(lon), !is.na(lat)) %>% { # BUGBUG
      plot.fun <- partial(plot.census.hhs, 
                          .data = ., 
                          .rct.villages.data = all.villages %>% filter(targeted.village),
                          use.google.maps = FALSE) 
      
      
      map(unique(.$cluster.id), ~ plot.fun(.cluster.id = .))
    }
}


problem.clusters <- unmatched.hh.data %>% 
  filter(group.dist == 50, unmatched.census.prop > 0.5 | unmatched.pre.census.prop > 0.5) %$% 
  prox.cluster.id %>% 
  unique

grouping.polygons.data <- grouping.data %>% 
  filter(prox.cluster.id %in% problem.clusters, group.dist == 50) %>% {
  left_join(plyr::ddply(., 
                  .(prox.cluster.id, hh.prox.group, group.dist), 
                  . %>% get.group.polygons(convert.to.wgs.84 = FALSE)),
            distinct(., hh.prox.group, group.dist, .keep_all = TRUE) %>% 
              select(-c(lon, lat)), 
            by = c("prox.cluster.id", "hh.prox.group", "group.dist"))
}

plot.simple.prox.groups <- function(cluster.group.polygons) {
  cluster.bbox <- make_bbox(lon, lat, data = filter(cluster.group.polygons, !ignore.centering))
  mat.cluster.bbox <- matrix(cluster.bbox, nrow = 2)
  
  cluster.group.polygons %>% 
    mutate(mix.cat = factor(mix.cat, 
                            levels = c("pure.pre.census", "mixed", "pure.census"),
                            labels = c("Pre-census Only", "Both", "Census Only"))) %>% 
  ggplot() +
    coord_fixed(xlim = mat.cluster.bbox[1, ], ylim = mat.cluster.bbox[2, ]) +
    geom_polygon(aes(lon, lat, group = hh.prox.group, color = mix.cat, fill = mix.cat), alpha = 0.5) +
    labs(x = "", y = "") +
    scale_color_discrete("") +
    scale_fill_discrete("") +
    ggtitle(sprintf("Cluster %d", cluster.group.polygons$prox.cluster.id[1])) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_wrap(~ group.dist)
}

plot.prox.groups <- function(cluster.hhs, .group.polygons, .zoom = 15) {
  hhs.center <- cluster.hhs %>% 
    filter(!is.na(cluster.id)) %>% 
    make_bbox(lon, lat, data = .) %>% 
    matrix(nrow = 2) %>% 
    rowMeans()
  
  plot.obj <- ggmap(get_googlemap(center = hhs.center,
                                  maptype = "hybrid",
                                  zoom = .zoom, 
                                  scale = 2,
                                  style = "element:labels|visibility:off", # Drop all labels from map
                                  key = config$google_api_key)) +
    geom_point(aes(lon, lat, fill = survey.type, shape = is.na(cluster.id)), color = "white", alpha = 1, size = 2, stroke = 1, data = cluster.hhs) +
    ggtitle(sprintf("Cluster %d", cluster.hhs$cluster.id[1])) +
    scale_fill_discrete("Survey Type") +
    scale_shape_manual("Correct Cluster ID", values = c(21, 24)) +
    labs(x = "", y = "") +
    theme(axis.text = element_blank())
  
  if (!is.null(.group.polygons)) {
    plot.obj <- plot.obj +
      geom_polygon(aes(lon, lat, group = factor(hh.prox.group)), color = "orange", alpha = 0.1, size = 1, 
                   data = filter(.group.polygons, prox.cluster.id == cluster.hhs$prox.cluster.id[1])) 
  }
  
  plot(plot.obj)
}

grouping.plots <- grouping.polygons.data %>% 
  plyr::dlply(.(group.dist, prox.cluster.id), plot.simple.prox.groups) 


## Demographics

### Phones

hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  ggplot() +
  geom_bar(aes(hhh_have_phone)) +
  labs(x = "Head of household has phone")

census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  ggplot() +
  geom_bar(aes(have_phone)) + 
  labs(x = "household member has phone")

mean.phone.own.data <- census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  group_by(cluster.id) %>% 
  summarize(mean.phone.own = mean(have_phone == "Yes")) 

mean.phone.own.data %>% {
    ggplot(.) +
      geom_histogram(aes(mean.phone.own), binwidth = 0.1, alpha = 0.5, color = "black") +
      geom_vline(xintercept = median(.$mean.phone.own), color = "red", linetype = "dashed")
  }


### Gender

census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  mutate(gender = factor(gender, levels = 1:2, labels = c("Male", "Female"))) %>% 
  ggplot() +
  geom_bar(aes(gender)) 

### Age

census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  mutate(know_age = factor(know_age, levels = 1:4, labels = c("Exact", "ID", "Roughly", "Doesn't Know/Guess"))) %>% 
  ggplot() +
  geom_bar(aes(know_age)) 

census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  ggplot() +
  geom_histogram(aes(age), binwidth = 5) +
  scale_x_continuous(breaks = seq(10, 100, 10)) 


## Geographic Characteristics

### Village-PoT Distances
factor.dist.pot <- . %>% { ifelse(. <= 2500/2, "close", "far") %>% factor }

village.centers <- hh.census.data %>%
  filter(!is.na(cluster.id), !is.na(wave), !is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  group_by(wave, county, cluster.id) %>%
  summarise(geometry = st_centroid(st_combine(geometry)), .groups = "drop") %>%
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry() %>%
  left_join(select(pot.info, wave, county, cluster.id, starts_with("alt.pot."), ends_with(".verify")),
            c("wave", "county", "cluster.id"))


get.vill.pot.dist <- function(filtered.village.centers,
                              pot.coords = c("alt.pot.lon", "alt.pot.lat")) {
  village.sf <- st_as_sf(filtered.village.centers, coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(kenya.proj4)

  pot.sf <- st_as_sf(filtered.village.centers, coords = pot.coords, crs = 4326) %>%
    st_transform(kenya.proj4)

  dist.mat <- st_distance(village.sf, pot.sf) %>%
    units::drop_units()

  filtered.village.centers %>%
    mutate(
      dist.to.own.pot = diag(dist.mat),
      dist.pot.group = factor.dist.pot(dist.to.own.pot),
      closest.other = `diag<-`(dist.mat, NA) %>%
        apply(1, min, na.rm = TRUE)
    ) %>%
    bind_cols(seq(2500/2, 7500, 2500/2) %>% {
      setNames(map(., ~ rowSums(dist.mat < ., na.rm = TRUE)), paste0("num.within.", .))
    })
}

vill.pot.dist <- village.centers %>%
  filter(cluster.id %in% known.pot.locations$cluster.id) %>%
  get.vill.pot.dist()



verify.vill.pot.dist <- village.centers %>%
  filter(cluster.id %in% known.pot.locations$cluster.id, !is.na(lon.verify), !is.na(lat.verify)) %>%
  get.vill.pot.dist(pot.coords = c("lon.verify", "lat.verify"))

# save(vill.pot.dist, verify.vill.pot.dist, village.centers, file = file.path("..", "data", "takeup_village_pot_dist.RData"))

vill.pot.dist %>% 
  right_join(verify.vill.pot.dist, c("wave", "county", "cluster.id"), suffix = c(".origin", ".verify")) %>% 
  select(wave, county, cluster.id, matches("dist.pot.group"), matches("dist.to.own")) %>% 
  left_join(select(pot.info, cluster.id, verify.dist), "cluster.id") %>% 
  filter(wave == 2, dist.pot.group.origin != dist.pot.group.verify)

verify.vill.pot.dist %>% { 
  ggplot(.) +
    geom_histogram(aes(dist.to.own.pot), binwidth = 2500/4, boundary = 0, alpha = 0.5, color = "black") + 
    scale_x_continuous("Distance to Own PoT", breaks = seq(0, 5000, 2500/4)) + 
    scale_y_continuous(breaks = seq(0, 100, 2)) 
} 

vill.pot.dist %>% 
  select(starts_with("num.within"), cluster.id, lon, lat, dist.to.own.pot) %>% 
  gather(within.m, num.pot, -c(cluster.id, lon, lat, dist.to.own.pot)) %>% 
  separate(within.m, into = c("temp", "within.m"), sep = "num.within.") %>% 
  filter(within.m <= 5000) %>% 
  select(-temp) %>% {
  ggplot(.) +
    geom_freqpoly(aes(num.pot, color = within.m), binwidth = 1) +
    scale_x_continuous("Number of PoT", breaks = 0:10) +
    scale_y_continuous(breaks = seq(0, 200, 20)) +
    scale_color_discrete("Within (meters)")
  }

plot.marginal.dist.to.other <- function(.data, .color = NULL) {
  .data %>% 
    mutate(marginal.dist = closest.other - dist.to.own.pot,
           dist.to.own.pot.cat = cut(dist.to.own.pot, 
                                     breaks = c(0,  2500/2, 10000),
                                     labels = c("Close", "Far"))) %>% {
    plot.obj <- ggplot(.) +
      geom_histogram(aes_string("marginal.dist", group = .color, fill = .color), alpha = 0.5, color = "black", binwidth = 500, boundary = 0, position = "dodge") +
      scale_y_continuous(breaks = seq(0, 50, 2)) +
      scale_x_continuous("Marginal Distance to Other PoT", breaks = seq(-1000, 10000, 1000))
                                     
    plot(plot.obj + 
           geom_vline(xintercept = quantile(.$marginal.dist, probs = c(0.25, .5, 0.75)), 
                      color = "red", linetype = "dashed"))
    
    plot(plot.obj + facet_wrap(~ dist.to.own.pot.cat))
  }
}

vill.pot.dist %>% plot.marginal.dist.to.other


### Intra-cluster Household Distance

get.hh.dist <- function(village.hhs) {
  st_as_sf(village.hhs, coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(kenya.proj4) %>%
    st_distance() %>%
    units::drop_units() %>% {
      upper.map <- upper.tri(.)
      tibble(hh.dist = magrittr::extract(., upper.map),
             from = village.hhs$KEY %>%
               magrittr::extract(-length(.)) %>%
               rep(times = seq(length(.), 1)),
             to = village.hhs$KEY %>%
               matrix(ncol = length(.), nrow = length(.), byrow = TRUE) %>%
               magrittr::extract(upper.map))
    }
}

hh.pair.distance <- hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  group_by(cluster.id) %>% 
  filter(!is.na(lon), !is.na(lat)) %>%
  do(get.hh.dist(.)) %>% 
  ungroup

hh.med.dist <- hh.pair.distance %>% 
  mutate(temp = from, from = to, to = temp) %>% 
  select(-temp) %>% 
  bind_rows(hh.pair.distance) %>% 
  group_by(cluster.id, from) %>% 
  summarize(median.dist = median(hh.dist))


hh.pair.distance %>% 
  ggplot() +
  geom_density(aes(x = hh.dist)) +
  coord_cartesian(xlim = c(0, 5000)) + 
  theme(legend.position = "none")

### Density

cluster.density <- hh.pair.distance %>% 
  mutate(temp = from, from = to, to = temp) %>% 
  select(-temp) %>% 
  bind_rows(hh.pair.distance) %>% 
  arrange(from, hh.dist) %>% 
  group_by(cluster.id, from) %>% 
  filter(row_number() <= 50) %>% 
  summarize(med.dist = median(hh.dist)) %>% 
  group_by(cluster.id) %>% 
  summarize(med.dist = median(med.dist)) %>% 
  ungroup 

ggplot(cluster.density) +
  geom_histogram(aes(med.dist), binwidth = 25, boundary = 0, alpha = 0.5, color = "black") +
  geom_vline(xintercept = quantile(cluster.density$med.dist, prob = c(0.25, 0.5, 0.75)), color = "red", linetype = "dashed") +
  scale_x_continuous("", breaks = seq(0, 700, 50))

# Study Randomization and Sampling

wave.1.hh.data <- hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  filter(wave == 1) %>% 
  select(cluster.id, KEY, tu.hh.id) 

wave.1.individuals <- wave.1.hh.data %>% 
  distinct(KEY) %>% 
  left_join(census.data, "KEY")  %>% 
  select(cluster.id, KEY, KEY.individ, tu.person.id, name1st, check1st, name_mid, check_mid, name2nd, check2nd, clan, clan_check, two_digits, two_digits_check, age, have_phone) 

hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  count(cluster.id) %>% {
    summary(.$n) %>% print
    mutate(., cluster.size = cut(n, breaks = seq(0, 300, 25), include.lowest = TRUE)) %>% 
      count(cluster.size) %>% 
      print
    ggplot(.) +
      geom_histogram(aes(n), color = "black", alpha = 0.5, binwidth = 25, boundary = 25) +
      scale_y_continuous(breaks = seq(0, 30, 2), limits = c(0, 30)) +
      scale_x_continuous("", breaks = seq(0, 300, 25)) 
  }

## Stratified Cluster Randomization

# http://blogs.worldbank.org/impactevaluations/tools-of-the-trade-doing-stratified-randomization-with-uneven-numbers-in-some-strata

arms <- c("control", "ink", "bracelet", "calendar")



cluster.strat.data <- hh.census.data %>% 


  filter(!is.na(cluster.id)) %>% 
  group_by(cluster.id, county) %>% 
  summarize(num.hh = n()) %>% 
  ungroup %>%
  left_join(vill.pot.dist %>% select(-county), "cluster.id") %>% 
  left_join(cluster.density, "cluster.id") %>% 
  left_join(mean.phone.own.data, "cluster.id") %>% 
  mutate(county.group = ifelse(county == "Kakamega", "Kakamega", "Busia-Siaya") %>% factor,
         density.group = ifelse(med.dist <= median(med.dist), "high", "low") %>% factor) %>% 
  group_by(county, dist.pot.group) %>% 
  mutate(assigned.treatment = sample(c(rep(arms, times = n() %/% length(arms)), sample(arms, n() %% length(arms))))) %>% 
  ungroup


## Re-creating treatment assignment
#### Setup
wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"
raw.data.path <- . %>% here("data", "raw-data", .)

datetime.format <- "%b %d, %Y %I:%M:%S %p"
takeup.datetime.type <- col_datetime(datetime.format)
takeup.date.type <- col_date(datetime.format)
rct.counties <- c("Busia", "Siaya", "Kakamega")
busia.subcounties <- c("butula", "nambale", "teso south", "teso north") 
siaya.subcounties <- c("gem", "ugenya", "ugunja")

ke.lvl2.adm.data <- read_rds(here("data", "adm", "KEN_adm2.rds"))

counties.adm.data <- ke.lvl2.adm.data[ke.lvl2.adm.data$NAME_1 %in% rct.counties, ] #, "Vihiga"), ]
subcounties.adm.data <- counties.adm.data[!counties.adm.data$NAME_1 %in% c("Busia", "Siaya") | counties.adm.data$NAME_2 %in% str_to_title(c(busia.subcounties, siaya.subcounties)), ]

county.bbox <- counties.adm.data@bbox
subcounty.bbox <- subcounties.adm.data@bbox

#### Loading data

# First, generate hh.census.data
validate.coords <- . %>% 
  mutate(invalid.coord = 
           (!is.na(lon) & (lon > county.bbox["x", "max"] | lon < county.bbox["x", "min"])) |
           (!is.na(lat) & (lat > county.bbox["y", "max"] | lat < county.bbox["y", "min"])))

tu.data.reader <- function(file.name, submit.datetime.type = NULL, .other.types = NULL) { # =  "%b %d, %Y %I:%M:%S %p") {
  col.types <- list(SubmissionDate = if (is.null(submit.datetime.type)) takeup.datetime.type else submit.datetime.type,
                                       manual_long = col_number(),
                                       manual_lat = col_number()) %>% 
    c(.other.types)
  
  read_csv(file.name, col_types = col.types) %>% 
    mutate(isValidated = isValidated == "true") %>% 
    rename(lat = `gps-Latitude`,
           lon = `gps-Longitude`,
           cluster.id = cluster_id) %>% 
    filter(deviceid != "(web)") %>% 
    mutate(lon = ifelse(is.na(lon), manual_long, lon),
           lat = ifelse(is.na(lat), manual_lat, lat)) %>% 
    validate.coords()
}

hh.census.data.pre24th <- tu.data.reader(here("data", "census_gps_util23.csv"), 
                                         submit.datetime.type = col_datetime("%d/%m/%Y %H:%M")) %>% 
  filter(SubmissionDate >= "2016-08-22") %>% 
  mutate_at(vars(present, return, number_individuals, instance, check), as.numeric)

factor.have.phone <- . %>% factor(levels = 0:2, labels = c("No", "Yes", "Don't know number"))

   
hh.census.data <- tu.data.reader(raw.data.path("Census.csv"), 
                                 .other.types = list(`gps-Latitude` = col_double(),
                                                     `gps-Longitude` = col_double(),
                                                     `gps-Altitude` = col_double(),
                                                     `gps-Accuracy` = col_double(),
                                                     deviceid = col_character(),
                                                     gps_work = col_integer())) %>% 
  filter(SubmissionDate >= "2016-08-25") %>%
  select(!c(starts_with("hhh_two_digits"))) |> 
  bind_rows(select(hh.census.data.pre24th, !c(hhh_name1st_check, directions, starts_with("hhh_two_digits")))) %>% 
  select(-c(RI, matches("^IC\\d{1,2}$"))) %>% 
  left_join(rct.village.codes, by = "village") %>%
  mutate(village_name = factor(village_name),
         present = ifelse(is.na(present) & SubmissionDate <= "2016-08-22", NA, present),
         return = as.numeric(return),
         hhh_have_phone = factor.have.phone(hhh_have_phone))

missing.gps.data <- read_csv(here("data", "missing_gps.csv")) %>% 
  select(KEY, `gps-Latitude`, `gps-Longitude`, gps_work, manual_long, manual_lat) %>% 
  mutate(lon = ifelse(is.na(`gps-Longitude`), manual_long, `gps-Longitude`),
         lat = ifelse(is.na(`gps-Latitude`), manual_lat, `gps-Latitude`),
         data.source = "missing_gps",
         across(c("lon", "lat"), as.double)) %>% 
  select(-c(starts_with("gps-"), starts_with("manual_"))) %>% 
  filter(!is.na(KEY), !is.na(lon), !is.na(lat))  

missing.gps.149.data <- read_csv(here("data", "missing_gps_149.csv")) %>% 
  transmute(lon = gpslongitude,
            lat = gpslatitude, 
            KEY = `key of form with missing GPS data`,
            data.source = "missing_gps_149")

missing.gps.360.data <- read_csv(here("data", "missing_gps_360.csv")) %>% 
  transmute(lon = gpslongitude,
            lat = gpslatitude, 
            KEY = `key of form missing GPS data`,
            data.source = "missing_gps_149")

missing.gps.844.data <- read_csv(here("data", "missing_gps_844.csv")) %>%
  rename(key.gps = `Key with GPS`,
         key.missing.gps = `Key missing GPS`) %>%
  inner_join(hh.census.data, c("key.gps" = "KEY")) %>% 
  transmute(lon, lat, KEY = key.missing.gps, data.source = "missing_gps_844")

hh.census.data %<>% {
  mask <- .$KEY %in% missing.gps.844.data$KEY

  mutate_at(., vars(lon, lat), funs(ifelse(mask, NA, .)))  
}

missing.gps.data %<>% 
  bind_rows(missing.gps.149.data, missing.gps.360.data, missing.gps.844.data) %>% 
  validate.coords()

hh.census.data %<>% 
  left_join(missing.gps.data, by = "KEY") %>% 
  filter(is.na(data.source) | data.source != "missing_gps_149" | cluster.id == 149) %>%
  mutate(lon = ifelse(is.na(lon.x), lon.y, lon.x),
         lat = ifelse(is.na(lat.x), lat.y, lat.x),
         gps_work = ifelse(is.na(gps_work.y), gps_work.x, as.numeric(gps_work.y)),
         invalid.coord = ifelse(is.na(invalid.coord.y), invalid.coord.x, invalid.coord.y)) %>% 
  select(-matches("\\.[xy]$"))

hh.census.data %<>% 
  anti_join(plyr::ldply(paste0("data/", c("takeup_returns_toremove.csv", 
                    "takeup_moved_notvalid_toremove.csv",
                    "takeup_noconsent_toremove.csv")), read_csv),
            by = c("KEY" = "key"))

hh.id.dict <- read_rds("data/takeup_hh_id_dict.rds")

hh.census.data %<>% 
  select(-county) %>% 
  left_join(cluster.wave.county.data, "cluster.id") %>% 
  left_join(hh.id.dict, "KEY")

person.id.dict <- read_rds("data/takeup_person_id_dict.rds")

clean.names <- . %>% 
  mutate_at(vars(matches("1st|mid|2nd|clan")), str_to_upper) %>% 
  mutate_at(vars(matches("1st|mid|2nd|clan")), 
            funs(ifelse(. %in% c("NONE", "N/A", "NA", "NO"), NA, .))) %>% 
  mutate_at(vars(matches("1st|2nd|clan")), 
            funs(ifelse(nchar(.) < 2, NA, .))) %>% 
  mutate(last_name = ifelse(is.na(name2nd) & !is.na(name_mid), name_mid, name2nd),
         name_mid = ifelse(is.na(name2nd) & !is.na(name_mid), NA, name_mid)) 

census.data <- hh.census.data %>% {
  individ.data <- read_csv(raw.data.path("Census-survey-individual.csv"))
  anti_join(., individ.data, c("KEY" = "PARENT_KEY")) %>% 
    count(SubmissionDate) %>% 
    rename(num.hhs.without.individuals = n) %>% 
    knitr::kable(col.names = c("Submission Date", "Number HH w/o individuals")) %>% 
    print
  
  inner_join(., individ.data, c("KEY" = "PARENT_KEY"), suffix = c(".hh", ".individ"))
} %>% 
  mutate(two.digit.match = ifelse(have_phone == 1, two_digits == two_digits_check, NA),
         have_phone = factor.have.phone(have_phone)) %>% 
  group_by(KEY) %>% 
  mutate(num.individuals = n(),
         hh.has.phone = any(have_phone == "Yes"),
         hh.has.non.phone = any(have_phone != "Yes")) %>% 
  ungroup %>% 
  clean.names %>% 
  filter(!is.na(name1st) & (!is.na(name_mid) | !is.na(name2nd))) %>% 
  left_join(person.id.dict, "KEY.individ")

#This step would also remove all households with no individuals (inner join)
hh.census.data %<>%
  inner_join(distinct(census.data, KEY, num.individuals, hh.has.phone, hh.has.non.phone), "KEY")















### Old code from Karim's script to recreate treatment assignment
census_data_env = new.env()
load("data/takeup_census.RData", envir = census_data_env)

vill_pot_env = new.env()
load("data/takeup_village_pot_dist.RData", envir = vill_pot_env)

census_data_env$hh.census.data
vill_pot_env$vill.pot.dist

census_data_env$hh.census.data %>%
  count(cluster.id)

# cluster.strat.data <- hh.census.data %>% 
#   filter(!is.na(cluster.id)) %>% 
#   group_by(cluster.id, county) %>% 
#   summarize(num.hh = n()) %>% 
#   ungroup %>%
#   left_join(vill.pot.dist, "cluster.id") %>% 
#   left_join(cluster.density, "cluster.id") %>% 
#   left_join(mean.phone.own.data, "cluster.id") %>% 
#   mutate(county.group = ifelse(county == "Kakamega", "Kakamega", "Busia-Siaya") %>% factor,
#          density.group = ifelse(med.dist <= median(med.dist), "high", "low") %>% factor) %>% 
#   group_by(county, dist.pot.group) %>% 
#   mutate(assigned.treatment = sample(c(rep(arms, times = n() %/% length(arms)), sample(arms, n() %% length(arms))))) %>% 
#   ungroup



#### Saved treatment assignment from Karim

cluster_strat_data = read_rds("data/takeup_cluster_randomization_1.0.rds")

join_cluster_strat_data = read_rds("data/takeup_processed_cluster_strat.rds")

join_cluster_strat_data %>%
  colnames()

cluster_strat_data %>%
  select(cluster.id, assigned.treatment, dist.pot.group)

join_cluster_strat_data %>%
  select(cluster.id, assigned.treatment, dist.pot.group)

