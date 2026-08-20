library(tidyverse)

source("rct-design-fieldwork/takeup_rct_assign_clusters.R")

# Study Randomization and Sampling

```{r wave-1-data}
wave.1.hh.data <- hh.census.data %>% 
  filter(!is.na(cluster.id)) %>% 
  filter(wave == 1) %>% 
  select(cluster.id, KEY, tu.hh.id) 

wave.1.individuals <- wave.1.hh.data %>% 
  distinct(KEY) %>% 
  left_join(census.data, "KEY")  %>% 
  select(cluster.id, KEY, KEY.individ, tu.person.id, name1st, check1st, name_mid, check_mid, name2nd, check2nd, clan, clan_check, two_digits, two_digits_check, age, have_phone) 
```

```{r, fig.width=8}
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
```

## Stratified Cluster Randomization

http://blogs.worldbank.org/impactevaluations/tools-of-the-trade-doing-stratified-randomization-with-uneven-numbers-in-some-strata

```{r}
cluster.density %>% 
  left_join(num.hh.data, "cluster.id") %>% 
  ggplot(aes(n.census, med.dist)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Number of Households", y = "Median Distance Between Households")

cluster.density %>% 
  left_join(mean.phone.own.data, "cluster.id") %>% 
  ggplot(aes(mean.phone.own, med.dist)) +
  geom_point() +
  geom_smooth() +
  labs(x = "Proportion of Phone Ownership", y = "Median Distance Between Households")
```

```{r stratied-cluster-randomization, eval=FALSE}
arms <- c("control", "ink", "bracelet", "calendar")


cluster.strat.data <- hh.census.data %>% 


  filter(!is.na(cluster.id)) %>% 
  group_by(cluster.id, county) %>% 
  summarize(num.hh = n()) %>% 
  ungroup %>%
  left_join(vill.pot.dist, "cluster.id") %>% 
  left_join(cluster.density, "cluster.id") %>% 
  left_join(mean.phone.own.data, "cluster.id") %>% 
  mutate(county.group = ifelse(county == "Kakamega", "Kakamega", "Busia-Siaya") %>% factor,
         density.group = ifelse(med.dist <= median(med.dist), "high", "low") %>% factor) %>% 
  group_by(county, dist.pot.group) %>% 
  mutate(assigned.treatment = sample(c(rep(arms, times = n() %/% length(arms)), sample(arms, n() %% length(arms))))) %>% 
  ungroup
```

```{r}
cluster.strat.data <- read_rds("../data/takeup_cluster_randomization_1.0.rds") %>% 
  select(-county) %>%
  mutate(assigned.treatment = factor(assigned.treatment, levels = c("control", "ink", "calendar", "bracelet"))) %>% 
  # left_join(distinct(hh.census.data, cluster.id, wave), "cluster.id")  
  left_join(cluster.wave.county.data, "cluster.id")  

write_rds(cluster.strat.data, "../data/takeup_processed_cluster_strat.rds")  

hh.census.data %<>% 
  left_join(select(cluster.strat.data, wave, county, cluster.id, assigned.treatment, dist.pot.group), 
            c("wave", "county", "cluster.id"))

census.data %<>% 
  # left_join(select(cluster.strat.data, cluster.id, assigned.treatment), "cluster.id")
  left_join(select(cluster.strat.data, wave, county, cluster.id, assigned.treatment), 
            c("wave", "county", "cluster.id"))

cluster.strat.data %>% count(assigned.treatment)
cluster.strat.data %>% count(county, assigned.treatment)
cluster.strat.data %>% count(county, dist.pot.group, assigned.treatment)
```

```{r, eval=FALSE}
cluster.strat.data %>% 
  select(cluster.id, county, assigned.treatment) %>% 
  write_csv("takeup_cluster_treatment_assignment_1.0.csv")
```

```{r all-villages-csv, eval=FALSE}
cluster.strat.data %>% 
  left_join(mutate(cluster.survey.data, cluster.id = as.integer(cluster.id)), "cluster.id") %>%
  left_join(rct.villages, c("cluster.id", "target.village.id")) %>% 
  transmute(cluster.id = as.integer(cluster.id), 
            target.village_name = target.village_name.x, 
            targeted = !is.na(target.village_name.y)) %>% 
  write_csv("takeup_all_village.csv")
```

