library(tidyverse)
source("R/structural/legacy-utils.R")
library(microbenchmark)
library(truncnorm)
library(testthat)
library(cmdstanr)
library(rstan)
library(posterior)
library(tidybayes)
library(nleqslv)



r_data = read_rds(
  "temp-data/r-predicted-data.rds"
)

source("optim/optim-functions.R")

r_pred_df = r_data$pred_df %>%
  mutate(
    u_sd = sqrt(total_error_sd^2 - 1)
  )

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)
# stick to monitored sms.treatment group
# remove sms.treatment.2
monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data %>% 
mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group,
    sms_treatment = factor(sms.treatment.2))


# fit_file = "data/stan_analysis_data/dist_fit71_STRUCTURAL_LINEAR_U_SHOCKS-1.csv"

# fit = as_cmdstan_fit(fit_file)


# To find fixed point need:
# benefit_cost
# mu_rep
# total_error_sd
# u_sd



# bc_draws = fit$draws(
#   variables = c(
#     "structural_cluster_benefit_cost", 
#     "cluster_dist_cost", 
#     "obs_cluster_mu_rep",
#     "total_error_sd[1]", 
#     "u_sd[1]")
# )


# bc_draws %>%
#   saveRDS("temp-data/temp-bc-draws.rds")

bc_draws = read_rds("temp-data/temp-bc-draws.rds")   
  

# rm(fit)
# gc()
# stop()



bc_draw_df = bc_draws %>%
  as_draws_df()

bc_rvar_df = bc_draw_df %>%
  spread_rvars(
    structural_cluster_benefit_cost[cluster_id], 
    cluster_dist_cost[cluster_id],
    obs_cluster_mu_rep[cluster_id],
    total_error_sd, 
    u_sd
    ) %>%
    left_join(
      analysis_data %>%
          select(
            cluster_id,
            assigned_dist_group, 
            assigned_treatment = assigned.treatment, 
            dist = cluster.dist.to.pot) %>%
            unique(), 
          by = "cluster_id"
    )  


draw_bc_rvar_df = bc_rvar_df %>%
  unnest_rvars()


draw_bc_rvar_df %>%
  group_by(assigned_treatment, assigned_dist_group) %>%
  summarise(
    across(
      c(structural_cluster_benefit_cost, obs_cluster_mu_rep, total_error_sd, u_sd),
      mean
    )
  )  %>%
filter(assigned_dist_group == "close") %>%
filter(assigned_treatment == "control") %>%
pivot_longer(-c(assigned_treatment, assigned_dist_group))  %>%
ungroup() %>%
select(-assigned_treatment, -assigned_dist_group)

draw_bc_rvar_df %>%
  group_by(assigned_treatment, assigned_dist_group) %>%
  summarise(
    across(
      c(structural_cluster_benefit_cost, obs_cluster_mu_rep, total_error_sd, u_sd),
      mean
    )
  ) %>%
  write_csv(
    "temp-data/mean-bc-rvar.csv"
  )

## Distribution of vstars in control close

control_vstar = draw_bc_rvar_df %>%
    filter(assigned_dist_group == "close" & assigned_treatment == "control")  %>%
  mutate(
      R_v_star = pmap(
        list(
          dist, 
          structural_cluster_benefit_cost,
          obs_cluster_mu_rep, 
          total_error_sd, 
          u_sd
        ),
        ~find_v_star(
          distance = ..1, 
          b = ..2, 
          mu_rep = ..3, 
          total_error_sd = ..4, 
          u_sd = ..5, 
          bounds = c(-Inf, Inf)
      )
  )
)

control_vstar %>%
    unnest_wider(R_v_star) %>%
    ggplot(aes(
        x = v_star
    )) +
    geom_histogram(bins = 30) +
    labs(
        title = "Distribution of W* in Control Close",
        x = "W*",
        y = "Count"
    ) +
    theme_minimal() +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "red"
    )
ggsave("temp-data/dist-w*-close.png", width = 8, height = 6)




draw_bc_rvar_df %>%
  pull(dist) %>%
  unique() 

# library(furrr)
# plan(multicore, workers = 8)
draw_bc_rvar_df = draw_bc_rvar_df %>%
  arrange(dist) %>%
  filter(.draw == 1) %>%
  head(100) %>%
  filter(assigned_treatment == "control") %>%
  mutate(
      R_v_star = pmap(
        list(
          dist, 
          structural_cluster_benefit_cost,
          obs_cluster_mu_rep, 
          total_error_sd, 
          u_sd
        ),
        ~find_v_star(
          distance = ..1, 
          b = ..2, 
          mu_rep = ..3, 
          total_error_sd = ..4, 
          u_sd = ..5, 
          bounds = c(-Inf, Inf)
      )
  )
)

draw_bc_rvar_df %>%
  unnest_wider(R_v_star) %>%
  ungroup() %>%
  select(dist, structural_cluster_benefit_cost, obs_cluster_mu_rep, total_error_sd, u_sd, v_star)



library(tidyverse)
source("R/structural/legacy-utils.R")
library(microbenchmark)
library(truncnorm)
library(testthat)
library(cmdstanr)
library(rstan)
library(posterior)
library(tidybayes)
library(nleqslv)

source("optim/optim-functions.R")

distance = seq(from = 0, to = 5000, length.out = 20)



bc_mean_df = read_csv(
    "temp-data/mean-bc-rvar.csv"
) %>%
    filter(assigned_dist_group == "close") %>%
    filter(assigned_treatment == "control")

bc_mean_df %>%
    pivot_longer(c(-assigned_treatment, -assigned_dist_group))

tibble(
    dist = distance
) %>%
    mutate(
        v_star = map(
            dist,
            ~find_v_star(
                distance = .x,
                b = -0.412,
                mu_rep = 0.0876,
                total_error_sd = 1.23,
                u_sd = 0.623,
                bounds = c(-Inf, Inf)
            )
        )
    ) %>%
    unnest_wider(v_star)


draw_bc_rvar_df = draw_bc_rvar_df %>%
  mutate(
      R_v_star = future_pmap(
        list(
          dist, 
          structural_cluster_benefit_cost,
          obs_cluster_mu_rep, 
          total_error_sd, 
          u_sd
        ),
        ~find_v_star(
          distance = ..1, 
          b = ..2, 
          mu_rep = ..3, 
          total_error_sd = ..4, 
          u_sd = ..5, 
          bounds = c(-Inf, Inf)
      )
  )
)
