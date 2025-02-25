
library(here)
library(tidyverse)
library(sf)
library(magrittr)

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"


# Identify valid target villages and clusters -----------------------------
rct.cluster.selection <- read_rds(here("data", "rct_cluster_selection_2.0.rds"))
rct.targetable.schools <- read_rds(here("data", "rct_targetable_schools_2.0.rds"))
cluster.survey.data <- read_rds(here("data", "takeup_cluster_survey.rds"))


rct.schools.data = read_rds(here("data", "takeup_rct_schools.rds")) %>%
  st_as_sf(wgs.84) %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2])



clusters.to.drop <- c("99", "691", "892", "1293", "239") %>% # Problems with clusters at cluster survey stage; dropped from study
  c("201", "853", "1402") # These are are study cluster that are no longer valid based on distance to other PoT

rct.cluster.selection %<>% magrittr::extract(!.$cluster.id %in% clusters.to.drop, )

rct.cluster.selection@data %<>% left_join(rct.schools.data %>% 
                                            select(cluster.id, lon, lat)) %>% 
  rename(pot.lon = lon, pot.lat = lat)

rct.cluster.selection %<>% 
  `@`(data) %>% 
  left_join(cluster.survey.data %>% 
              distinct(cluster.id, found, num.valid.villages, gpslatitude, gpslongitude, gps2latitude, gps2longitude) %>% 
              filter(!duplicated(cluster.id)),
            by = "cluster.id") %>% 
  mutate(alt.pot.lat = gps2latitude,
         alt.pot.lon = gps2longitude,
         pot.lat = ifelse(!is.na(gpslatitude), #& found %in% c("Found", "Name Changed"), 
                          gpslatitude, 
                          ifelse(!is.na(gps2latitude), #& found %in% c("Closed", "Not Found"), 
                                 alt.pot.lat, 
                                 pot.lat)),
         pot.lon = ifelse(!is.na(gpslongitude), #& found %in% c("Found", "Name Changed"), 
                          gpslongitude, 
                          ifelse(!is.na(gps2longitude), #& found %in% c("Closed", "Not Found"), 
                                 alt.pot.lon, 
                                 pot.lon))) %>%
  as.data.frame %>%
  st_as_sf(coords = c("pot.lon", "pot.lat"), crs = wgs.84) %>% 
  st_transform(kenya.proj4)

#' Create Village Dist Matrix
#'
#' Find the distance between each village in the cluster set with all the other
#' villages I think (?)
create_village_dist_mat = function(cluster_villages, rct_cluster_selection) {
  villages_sf = st_as_sf(cluster_villages,
    coords = c("target.lon", "target.lat"),
    crs = wgs.84
  ) %>%
    st_transform(kenya.proj4)
  dist_mat = st_distance(
    rct_cluster_selection %>% st_transform(kenya.proj4),
    villages_sf
  ) %>%
    units::drop_units()  %>%
    set_colnames(cluster_villages$target.village.id) %>%
    as_tibble() %>%
    mutate(
      dist.from.cluster.id = rct_cluster_selection$cluster.id,
      cluster.id = cluster_villages$cluster.id[1],
      cluster.group = rct_cluster_selection$cluster.group
    )
  long_dist_mat = dist_mat %>%
    tidyr::gather(target.village.id, dist, -c(dist.from.cluster.id, cluster.id, cluster.group))

  comb_long_dist_mat = left_join(
    long_dist_mat,
    long_dist_mat %>%
      select(dist.from.cluster.id, target.village.id, dist),
    by = c("cluster.id" = "dist.from.cluster.id", "target.village.id"),
    suffix = c(".to.other.pot", ".to.pot")
  ) %>%
    filter(dist.from.cluster.id != cluster.id) %>% 
    select(-dist.from.cluster.id) %>% 
    group_by(target.village.id, cluster.group) %>% 
    filter(row_number(dist.to.other.pot) == 1)  %>%
    ungroup() 

  wide_comb_dist_mat = comb_long_dist_mat %>%
    spread(cluster.group, dist.to.other.pot) %>%
    mutate(
      valid.target.village = bracelet.airtime >= 3860 & control.ink >= 3000,
      target.village.id = as.integer(target.village.id)
    )
  return(wide_comb_dist_mat)
}


cluster_survey_data_dm = cluster.survey.data %>%
  filter(!is.na(target.lat), !is.na(target.lon))  %>%
  arrange(cluster.id, target.village.id) %>%
  group_by(cluster.id)  %>%
  group_nest() %>%
  mutate(
    data = map2(data, cluster.id, ~mutate(.x, cluster.id = .y)),
    dist_mat = map(data, ~ create_village_dist_mat(.x, rct.cluster.selection))
  ) 

clust_survey_matched_dm = cluster_survey_data_dm %>%
  mutate(
    output_data = map2(
      data,
      dist_mat,
      ~ left_join(.x, .y, by = c("cluster.id", "target.village.id"))
    )
  ) %>%
  select(-cluster.id) %>%
  unnest(output_data) %>%
  select(-data, -dist_mat) %>%
  mutate(
    dist.to.pot = dist.to.pot.x,
    valid.target.village = valid.target.village.x
    ) %>%
  select(-dist.to.pot.x, -dist.to.pot.y, -valid.target.village.x, -valid.target.village.y) 


clust_survey_matched_dm_school = clust_survey_matched_dm %>%
  left_join(
    .,
    rct.targetable.schools %>%
      filter(selected.targeted == TRUE) %>%
      select(
        selected.targeted,
        pot.cluster.id,
        cluster.dist.cat
      ),
      by = c("cluster.id" = "pot.cluster.id")
  ) %>%
  group_by(cluster.id)  %>%
  {
    if (first(.$village.dist.cat) == "far") arrange(., desc(dist.to.pot)) else arrange(., dist.to.pot)
  } %>% 
  filter(!duplicated(target.village_name))  %>%
  ungroup() %>%
  mutate(
    target.actual.village.dist.cat = ifelse(dist.to.pot <= (2500/2), "close", "far"),
    target.valid.dist.to.pot = target.actual.village.dist.cat == village.dist.cat,
    target.in.range = dist.to.pot <= 2500,
    dist.switchable = !target.valid.dist.to.pot & valid.target.village,
    valid.target.village = valid.target.village & target.in.range, 
    target.village.pop.size = ifelse(target.pop_households <= median(target.pop_households, na.rm = TRUE), "small", "big")
    )  %>%
  filter(!is.na(target.village.pop.size)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster.pop.size.strata = unique(target.village.pop.size) %>% { ifelse(length(.) == 2, "mixed", .) },
         actual.cluster.dist.cat = unique(target.actual.village.dist.cat) %>% { ifelse(length(.) == 2, "mixed", .) }) %>% 
  ungroup()


rct_clean_clust = clust_survey_matched_dm_school %>%
  select(-dist.switchable) %>%
  left_join(
    clust_survey_matched_dm_school %>% 
      mutate(backup.cluster = FALSE) %>% # cluster.id %in% rct.backup.clusters$cluster.id) %>% 
      group_by(cluster.id, backup.cluster, cluster.pop.size.strata, actual.cluster.dist.cat) %>% 
      summarize(any.buffer.valid.village = any(valid.target.village),
                dist.switchable = any(dist.switchable)) %>% 
      ungroup,
    by = "cluster.id"
  ) %>%
  left_join(
    rct.targetable.schools %>%
      filter(selected.targeted == TRUE) %>%
      select(pot.cluster.id, village.dist.cat, cluster.dist.cat),
    by = c("cluster.id" = "pot.cluster.id")
  ) %>%
  mutate(village.dist.cat = village.dist.cat.x) %>%
  select(-village.dist.cat.x, -village.dist.cat.y)  %>%
  mutate(
    village.dist.cat = ifelse(dist.switchable & !backup.cluster, ifelse(village.dist.cat == "close", "far", "close"), village.dist.cat),
    any.buffer.valid.village = (dist.switchable & !backup.cluster) | any.buffer.valid.village
  ) 



# Display RCT clusters ----------------------------------------------------
 
rct.cluster.selection %>%
  filter(!any.buffer.valid.village) %>%
  select(cluster.id, cluster.group) %>%
  left_join(cluster.survey.data) %>%
  select(cluster.id, matches("valid"), bracelet.airtime, control.ink, dist.to.pot, village.dist.cat, cluster.group, dist.switchable)

# Village selection -------------------------------------------------------


rct_clean_clust %>% 
  filter(valid.target.village, cluster.id %in% rct.cluster.selection$cluster.id)  %>%
  group_by(cluster.id) %>%
  sample_n(1) %>%
  ungroup() %>%
  select(
    cluster.id, target.village.id, dist.to.pot
  )

cluster.survey.data %>% 
  filter(valid.target.village, cluster.id %in% rct.cluster.selection$cluster.id)  %>%
  group_by(cluster.id) %>%
  sample_n(1) %>%
  ungroup() %>%
  select(
    cluster.id, target.village.id, dist.to.pot
  )

randomize_village_targeting = function(seed, cluster_survey_data, rct_cluster_selection) {
  set.seed(seed)
  cluster_survey_data %>%
    filter(valid.target.village, cluster.id %in% rct_cluster_selection$cluster.id) %>%
    group_by(cluster.id) %>%
    sample_n(1) %>%
    ungroup() %>%
    select(
      cluster.id, target.village.id, dist.to.pot, village.dist.cat
    ) %>%
    left_join(
      rct_cluster_selection %>%
        as_tibble() %>%
        select(
          cluster.id,
          cluster.group
        ),
        by = "cluster.id"
    ) %>%
    mutate(
      random_treat = case_when(
        cluster.group == "control.ink" ~ sample(c("control", "ink"), 1),
        cluster.group == "bracelet.airtime" ~ sample(c("bracelet", "calendar"), 1)
      )
    )
}

seeds = 1:100

sim_df = map_dfr(
  seeds,
  ~ randomize_village_targeting(.x, rct_clean_clust, rct.cluster.selection),
  .id = "seed"
)


sim_df %>%
  ggplot(aes(
    x = dist.to.pot,
    fill = village.dist.cat
  )) +
  geom_density(alpha = 0.5) +
  facet_wrap(~random_treat)

old_sim_df = map_dfr(
  seeds,
  ~ randomize_village_targeting(.x, cluster.survey.data, rct.cluster.selection),
  .id = "seed"
)

sim_df

sim_df
old_sim_df

rct_clean_clust %>%
  select(contains('assign'))