
wgs.84 = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 = "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

datetime.format = "%b %d, %Y %I:%M:%S %p"
takeup.datetime.type = col_datetime(datetime.format)
takeup.date.type = col_date(datetime.format)
raw.data.path = . %>% here("data", "raw-data", .)

rct.counties = c("Busia", "Siaya", "Kakamega")
busia.subcounties = c("butula", "nambale", "teso south", "teso north") 
siaya.subcounties = c("gem", "ugenya", "ugunja")

ke.lvl2.adm.data = read_rds(here("data", "adm", "KEN_adm2.rds"))

counties.adm.data = ke.lvl2.adm.data[ke.lvl2.adm.data$NAME_1 %in% rct.counties, ] #, "Vihiga"), ]


yes_no_factor = function(.col, .yes.no = 1:2) {
    .col %>% 
    factor(levels = c(.yes.no, 97:98), labels = c("yes", "no", "prefer not say", "DK"))
} 

multi_factor = function(.col, labels, levels, ...) {
  if (missing(levels)) {
    levels <- c(seq_along(labels), 97:99)
    labels = labels %>%
        c("prefer not say", "DK", "other")
  } 
  map(str_split(.col, " "), factor, levels = levels, labels = labels, ...)
}

prepare_baseline_endline_data = function(.data, .cluster.strat.data) {
  .data %>% 
    mutate(who_worms = multi_factor(who_worms, 
                                    labels = c("child", "adult", "sick", "healthy", "pregnant", "old", "everyone")), 
           effect_worms = multi_factor(effect_worms, 
                                       labels = c("stomachache", "malnourishment", "stop child growing", "tired", "diarrhea", 
                                                  "stop child school")),
           how_spread = multi_factor(how_spread, 
                                     labels = c("drinking dirty water", "open def", "swim/bath dirty water", "no shoes", "no washing hands",
                                                "sex", "uncooked food")),
           stop_worms = multi_factor(stop_worms, 
                                     labels = c("wash hands", "wearing shoes", "using toilets", "drink clean water", "medicine", "clean home")),
           when_treat = multi_factor(when_treat, 
                                     labels = c("every week", "every month", "every 2 months", "every 3 months", "every 6 months", 
                                                "every year", "never", "when symptoms", "hw says")),
           have_phone = factor(have_phone, levels = 0:1, labels = c("No", "Yes")),
           school = factor(school,
                           levels = c(1:15, 99, 97),
                           labels = c("Never gone to school", paste("Primary", 1:8), paste("Secondary", 1:4), "College", "University", "Other", "Prefer Not to Say")),
           floor = factor(floor, levels = c(1:3, 99, 97), labels = c("Cement", "Earth", "Tiles", "Other", "Prefer Not to Say")),
           ethnicity = factor(ethnicity, levels = c(1:7, 99, 97),
                              labels = c("Luo", "Luhya", "Kisii", "Kalengin", "Kikukyu", "Kamba", "Iteso", "Other", "Prefer Not to Say"))) %>% 
    mutate_at(vars(worms_affect, neighbours_worms_affect), list(~ yes_no_factor(., .yes.no = 1:0))) %>% 
    mutate_at(vars(spread_worms), yes_no_factor) %>% 
    mutate_at(vars(school, floor, ethnicity), list(~ recode_factor(., "Other" = NA_character_, "Prefer Not to Say" = NA_character_) %>% 
                                                     fct_drop)) %>% 
    mutate(ethnicity2 = fct_lump(ethnicity, other_level = "Other Ethnicities", prop = 0.05)) %>% 
    left_join(select(.cluster.strat.data, wave, county, cluster.id, assigned.treatment, dist.pot.group), c("wave", "county", "cluster.id")) 
}

prepare_endline_data = function(.data, .census.data, .cluster.strat.data) {
  .data %>% 
    filter(across(c(present, interview, consent), ~ !is.na(.x) & .x == 1)) %>% 
    arrange(KEY.individ, SubmissionDate) %>% 
    # new method - order early if not NA, otherwise order by date
    # arrange(KEY.individ, is.na(dworm_rate), SubmissionDate) %>%
    group_by(KEY.individ) %>% 
    filter(row_number() == 1) %>% # If more than one entry for an individual, take first one (there are 22 such individuals)
    ungroup() %>% 
    prepare_baseline_endline_data(.cluster.strat.data) %>% 
    mutate(across(c(know_deworm, chv_visit, flyer, any.sms.reported, hh_bracelet, hh_cal, cal_value),
                  yes_no_factor, .yes.no = 1:0)) %>% 
    mutate(across(c(treat_begin, days_available, treat_end), factor, levels = c(1, 98), c("knows", "DK"))) %>% 
    mutate(
      treat_begin_date = ymd(sprintf("2016-%d-%d", month_treat_begin, day_treat_begin)),
      treat_end_date = ymd(sprintf("2016-%d-%d", month_treat_end, day_treat_end)),
      #where_offered = labelled(where_offered, c("somewhere else" = 0, "home" = 3, "DK" = 98)) %>% as_factor)
      find_out = multi_factor(find_out, 
                              levels = c(1:9, 99), 
                              labels = c("friend", "family", "chv", "elder", "church", "flyer", "poster", "enumerator", "baraza",
                                         "other")),
      gift_choice = factor(gift_choice, levels = 1:3, labels = c("bracelet", "calendar", "neither"))
    ) %>%
    left_join(select(.census.data, KEY.individ, dist.to.pot, sms.ctrl.subpop, true.monitored, monitored), by = "KEY.individ")  
       # text_content = factor(text_content, levels = c(1:3, 99), labels = c("reminders", "when/where", "social info", "other")))
}
prepare_baseline_data = function(.data, .cluster.strat.data) {
  .data %>% 
    prepare_baseline_endline_data(.cluster.strat.data) %>% 
    mutate(more_less = factor(more_less, levels = 1:4, labels = c("more", "less", "no diff", "DK")),
           treated_when = factor(treated_when, 
                                 levels = c(1:9, 97), 
                                 labels = c("1-2 mon", "3-5 mon", "6-7 mon", "8-9 mon", "10-11 mon",
                                            "1 year", "2 year", "3 year", "4 year or more", 
                                            "prefer no say")),
           who_treated = factor(who_treated, 
                                levels = 1:3, 
                                labels = c("child", "adult", "both")),
           dworm_proportion = factor(dworm_proportion, 
                                      levels = c(1:6, 97:98),
                                      labels = c("few", "nearly half", "half", "more than half", "many", "all", 
                                                 "prefer not say", "DK")),
           ink_more_less = factor(ink_more_less, levels = c(1:3, 97:98), labels = c("more", "less", "same", 
                                                                                    "prefer not say", "DK"))) %>% 
    mutate_at(vars(treated, family_treated), yes_no_factor) %>% 
    mutate_at(vars(few_deworm, many_deworm, matches("(praise|stigma)_(immunize|clothes|deworm|defecate)[A-D]$")), 
              list(~ factor(., levels = 0:2, labels = c("no", "yes", "maybe")))) %>% 
    mutate_at(vars(treated_where, where_family_treated), 
              list(~ factor(., 
                          levels = c(1:4, 97:99), 
                          labels = c("bought", "school MDA", "non-school MDA", "hosp/clinic", 
                                     "prefer not say", "DK", "other"))))
}


tu_data_reader = function(file.name, submit.datetime.type = NULL, .other.types = NULL) {
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
identify_closest_cluster = function(.data, data.coords.formula = ~ lon + lat,  key.variable = "KEY") {
  coord.vars <- all.vars(data.coords.formula)

  .data %>%
    filter(rowSums(is.na(model.frame(data.coords.formula, data = ., na.action = NULL))) == 0) %>%
    (function (.located.data) {
      data.sf <- st_as_sf(.located.data, coords = coord.vars, crs = wgs.84) %>%
        st_transform(kenya.proj4)

      villages.sf <- st_as_sf(known.village.locations, coords = c("target.lon", "target.lat"), crs = wgs.84) %>%
        st_transform(kenya.proj4)

      dist.matrix <- st_distance(data.sf, villages.sf) %>%
        units::drop_units()

      map_dfr(seq_len(nrow(dist.matrix)), function(i) {
        dist.col <- dist.matrix[i, ]
        in.range.dist <- dist.col[dist.col <= 2000]

        if (length(in.range.dist) == 0) {
          tibble(min.dist = NA, closest.cluster = NA)
        } else {
          min.dist <- min(in.range.dist)
          closest.idx <- which(dist.col <= min.dist)
          tibble(min.dist = min.dist,
                 closest.village = known.village.locations$target.village.id[closest.idx],
                 closest.cluster = known.village.locations$cluster.id[closest.idx]) %>%
            distinct(closest.cluster, .keep_all = TRUE)
        }
      }) %>%
        bind_cols(.located.data, .) %>%
        select(all_of(c(key.variable, "min.dist", "closest.cluster")))
    })
}
