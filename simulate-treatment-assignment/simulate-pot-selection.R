library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(stringr)
library(haven)
library(purrr)
library(ggplot2)
library(broom)
library(parallel)
library(foreach)
library(sf)
library(here)

doParallel::registerDoParallel(cores = 12) 
# Constants ---------------------------------------------------------------

kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"
wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

rct.counties <- c("Busia", "Siaya", "Kakamega")
busia.subcounties <- c("butula", "nambale", "teso south", "teso north") 
siaya.subcounties <- c("gem", "ugenya", "ugunja")


# Geographic admin boundaries ---------------------------------------------
ke.lvl2.adm.data <- read_rds(here::here("data", "adm", "KEN_adm2.rds")) %>%
  st_as_sf() 
ke.lvl3.adm.data <- read_rds(here::here("data", "adm", "KEN_adm3.rds")) %>%
  st_as_sf()

counties.adm.data <- ke.lvl2.adm.data[ke.lvl2.adm.data$NAME_1 %in% rct.counties, ] #, "Vihiga"), ]
subcounties.adm.data <- counties.adm.data[!counties.adm.data$NAME_1 %in% c("Busia", "Siaya") | counties.adm.data$NAME_2 %in% str_to_title(c(busia.subcounties, siaya.subcounties)), ]

ke.ward.lvl.adm.data <- ke.lvl3.adm.data  

ke.ward.lvl.adm.data =  ke.ward.lvl.adm.data %>%
  st_as_sf() %>%
  # select(NAME_1, NAME_2, NAME_3) %>% 
  rename(county = NAME_1,
         constituency = NAME_2,
         ward = NAME_3)

rct.ward.lvl.adm.data <- ke.ward.lvl.adm.data %>% 
  magrittr::extract(.$county %in% rct.counties, ) %>% 
  magrittr::extract(!.$county %in% c("Busia", "Siaya") | .$constituency %in% str_to_title(c(busia.subcounties, siaya.subcounties)), ) 

rct.ward.lvl.adm.data$ward.id <- rownames(rct.ward.lvl.adm.data)
  

kakamega.wards <- ke.ward.lvl.adm.data %>% 
  magrittr::extract(.$county == "Kakamega", ) 


rct.schools.data = read_rds(here::here("data", "takeup_rct_schools.rds")) %>%
  st_as_sf(wgs.84) %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2])
# Locate schools based on ADM boundary data -------------------------------


# 1) Convert ward polygons to sf, project to kenya.proj4
wards_sf <- st_as_sf(rct.ward.lvl.adm.data) %>%
  st_transform(kenya.proj4)

# 2) Convert schools to sf, project to kenya.proj4
schools_sf <- st_as_sf(rct.schools.data) %>%
  st_transform(kenya.proj4)

# 3) For each ward, find the schools it contains:
#    - st_contains(ward, school) returns which points (schools) fall inside the polygon (ward).
#    - By using st_join(wards_sf, schools_sf, join = st_contains),
#      we attach each matching school's attributes to the ward's row.
ward_school_join <- st_join(
  x    = wards_sf,
  y    = schools_sf,
  join = st_contains
)

# 4) Keep just the columns you care about. Suppose wards have 'ward.id'
#    and schools have 'cluster.id' (like the original code).
school.ward.dict <- ward_school_join %>%
  select(ward.id, cluster.id) %>%
  # Remove any rows where a ward didn't match a school (NA in cluster.id).
  filter(!is.na(cluster.id)) %>%
  # Drop geometry to get a regular data frame (not sf).
  st_drop_geometry() %>%
  distinct()



# School centered clusters ------------------------------------------------
# Move to separate file, for profiling using package lineprof
source(here::here("simulate-treatment-assignment", "simulate-functions.R"))

# ED NOTE: I can't find village.data and can't exclude the pilot schools here...
targetable.schools = read_rds(here::here("data", "takeup_rct_schools.rds")) %>%
  st_as_sf(wgs.84) %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
  mutate(county = County)




school.buffer.radius <- 1000
min.area.frac <- 0.6
school.area <- pi * (school.buffer.radius^2)

seeds = sample(1:1000, 100)
library(furrr)
plan(multicore)

simmed_rct_clusters = future_map(
  seeds,
  ~ get.rct.clusters(rct.schools.data, NULL, rct.schools.data,
                   num.clusters = c(control.ink = 79, bracelet.airtime = 79),
                   ca.outer.radius = c(control.ink = 3000, bracelet.airtime = 4000), 
                   ca.inner.radius = 2500, 
                   cluster.group.order = "random",
                   cluster.size.tester = generate.cluster.num.schools.tester2_sf(targetable.schools, 
                                                                              school.buffer.radius = school.buffer.radius, 
                                                                              min.num.schools = 1,
                                                                              min.area.frac = min.area.frac),
                   plot.rct.clusters = FALSE, 
                   seed = .x
                   ) %>%
    magrittr::extract(.$selected, ),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
)
names(simmed_rct_clusters) = seeds
saveRDS(simmed_rct_clusters, here::here("data", "simmed_rct_clusters.rds"))
