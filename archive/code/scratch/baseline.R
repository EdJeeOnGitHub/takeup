library(tidyverse)


baseline_data = read_rds(file.path("data", "takeup_baseline_data.rds"))

all_endline_data = read_rds(file.path("data", "all_endline.rds"))




all_endline_data %>%  
  select(religion, ethnicity)

all_endline_data %>%
  count(ethnicity) %>%
  arrange(-n)

endline_data %>%
  select(contains("clan"))

baseline_data %>%
  select(religion, ethnicity) %>%
  unique()

all_endline_data

# Baseline + Endline
# merge endline:
  # - phone
  # - primary
  # - floor

  # top two from each
  # - ethnicity: proportion majority ethnicity
  # - religion: proportion majority religion



# Clean up know everyone can be infected
# knowledge variables + add more sample

all_endline_data %>%
  select(contains('floor'))

all_endline_data %>%
  select(contains('phone'))

all_endline_data %>%
  select(contains('school')) %>%
  count(school) %>%
  arrange(-n)




# 9,935 adults of which 2,250 adults are surveyed at baseline, 
# 3,750 adults surveyed at endline and 
# 3,935 adults selected for text messaging intervention.


print(str_glue("Endline data: {nrow(all_endline_data)}"))
print(str_glue("Baseline data: {nrow(baseline_data)}"))

nrow(baseline.data)

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


n_hh_df = census_data %>%
    group_by(cluster.id) %>%
    summarise(
        n = sum(num.individuals)
    )

census_data %>%
  colnames()

    census_data %>%
      summarise(mean(have_phone == "Yes", na.rm = TRUE))
