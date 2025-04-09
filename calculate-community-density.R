# Setup -------------------------------------------------------------------

library(magrittr)
library(tidyverse)
library(lubridate)
library(sf)
library(sp)

# library(rgeos)
# library(econometr)

source("analysis_util.R")

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

# Load the data ----------------------------------------------------------
load(file.path("data", "takeup_village_pot_dist.RData"))

## Ed note: this produces error "unknown input format" as of 2023/04/13 so 
## we use the RData file and some env munging instead
# This data was prepared in takeup_field_notebook.Rmd
# census.data <- read_rds(file.path("data", "takeup_census.rds")) %>% 
#   rename(census.consent = consent) # Rename this to reduce chance of error
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
census.data = with_env(load_census_function, census_data_env)() %>%
  rename(census.consent = consent) # Rename this to reduce chance of error


baseline.data <- read_rds(file.path("data", "takeup_baseline_data.rds"))
takeup.data <- read_rds(file.path("data", "takeup.rds"))
all.endline.data <- read_rds(file.path("data", "all_endline.rds"))
reconsent.data <- read_rds(file.path("data", "reconsent.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds")) # Study clusters metadata
sms.content.data <- read_rds(file.path("data", "takeup_sms_treatment.rds"))
pot.info <- read_rds(file.path("data", "pot_info.rds")) # Point of treatment information
wtp.data <- read_rds(file.path("data", "takeup_wtp.rds")) # Willingness-to-pay experiment data

# Geographic data for clusters in the study
pot_geo_info <- pot.info %>% 
  filter(!is.na(wave)) %>% 
  transmute(cluster.id, 
            pot.lon = lon.post.rct.update,
            pot.lat = lat.post.rct.update) 
  
# HH distance to PoT ------------------------------------------------------

# census.data %>% 
#   group_by(cluster.id) %>% 
#   group_modify(function(households, ...) {
#     hh_sf <- st_as_sf(households, coords = c("lon", "lat"), crs = wgs.84) %>% 
#       st_transform(kenya.proj4)
#     
#     # pot_sf <- st_as_sf(pot_geo_info, coords = c("pot.lon", "pot.lat"), crs = wgs.84)
#     
#     return(households)
#   }, pot_geo_info = pot_geo_info) %>% 
#   ungroup()
#   
#   
#   
#   nest_join(pot_geo_info, ., by = "cluster.id", name = "households") %>% 
#   # st_as_sf(coords = c("pot.lon", "pot.lat"), crs = wgs.84) %>% 
#   mutate(households = map(households, st_as_sf, coords = c("lon", "lat"), crs = wgs.84) %>% 
#            map2(geometry, ~ mutate(.x, 
#                                    dist.to.pot = st_distance(.x %>% st_transform(kenya.proj4),
#                                                              st_sfc(.y, crs = wgs.84) %>% st_transform(kenya.proj4)) %>% 
#                                      as.vector()))) 

# ── 2.  Point‑of‑treatment (PoT) geometry ------------------------------------
pot_sf <- pot_geo_info %>%                         # has cluster.id, pot.lon, pot.lat
  st_as_sf(coords = c("pot.lon", "pot.lat"),
           crs     = wgs.84, remove = FALSE) %>%    # keep lon/lat columns if you like
  st_transform(kenya.proj4) %>%                    # put into projected metres
  select(cluster.id, pot_geom = geometry)          # rename so we can keep two geoms

# ── 3.  Household geometry ----------------------------------------------------
hh_sf <- census.data %>%                           # has cluster.id, lon, lat
  st_as_sf(coords = c("lon", "lat"),
           crs     = wgs.84, remove = FALSE) %>% 
  st_transform(kenya.proj4)


pot_sf

hh_sf

# ── 4.  Join PoT geometry onto each household ---------------------------------
hh_sf <- hh_sf %>%
  left_join(pot_sf %>% st_drop_geometry(), by = "cluster.id")

  colnames(hh_sf)

  hh_sf %>%
    select(num.individuals) %>%
    skimr::skim()

# Calculate the convex hull for household points in each cluster

convex_hulls <- hh_sf %>%
  group_by(cluster.id) %>%                          # Group by cluster ID
  summarise(geometry = st_union(geometry), n_indiv = sum(num.individuals), n_hh = n()) %>%      # Combine all geometries in each cluster
  mutate(
    convex_hull = st_convex_hull(geometry),
    area = st_area(convex_hull),
    area_km2 = as.numeric(area) / 1e6,         # Convert area to km²
    density_km2 = n_indiv/ area_km2            # Calculate density
  )    # Calculate the convex hull

convex_hulls %>%
  ggplot(aes(
    x = density_km2
  )) +
  geom_histogram() 

sd(convex_hulls$density_km2, na.rm = TRUE)

dens_df = convex_hulls %>%
  filter(cluster.id %in% unique(takeup.data$cluster.id)) %>%
  summarise(
    mean_dens = mean(density_km2, na.rm = TRUE),
    sd_dens = sd(density_km2, na.rm = TRUE),
    prop_in_1sd = mean((density_km2 >= mean_dens - sd_dens) & (density_km2 <= mean_dens + sd_dens)),
    prop_in_2sd = mean((density_km2 >= mean_dens - 2*sd_dens) & (density_km2 <= mean_dens + 2*sd_dens)),
    n = n()
    ) 
library(kableExtra)

dens_tbl = dens_df %>%
  st_drop_geometry() %>%
  kbl(
    digits = 2,
    col.names = c("Mean density (km\\textsuperscript{2})", "SD density (km\\textsuperscript{2})", "Proportion in 1 SD", "Proportion in 2 SD", "N Communities"),
    format = "latex",
    caption = "Density of Communities",
    align = "c",
    escape = FALSE,
    booktabs = TRUE
  ) %>%
  kable_styling(
    latex_options = c("scale_down")
  )

dens_tbl %>%
  save_kable(
    file = "temp-data/density_table.tex",
    format = "latex"
  )

# View the resulting convex hulls
plot(convex_hulls)