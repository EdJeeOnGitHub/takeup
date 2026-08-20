library(brms)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(tidybayes)


hier_model = cmdstan_model("scratch/ols-ss-uni.stan", include_paths = "scratch/")


cov_analysis_data = read_csv("temp-data/analysis-cluster-covariate-data.csv") %>%
  mutate(assigned_dist_group = dist.pot.group) %>%
  mutate(
    cluster.id = cluster_id
  ) %>%
  mutate(assigned.treatment = factor(assigned.treatment, levels = c("control", "ink", "calendar", "bracelet"))) %>%
  mutate(assigned_treatment = assigned.treatment)  %>%
  mutate(wt = 1) %>%
  mutate(
      county = factor(county),
      cluster.id = factor(cluster.id),
      assigned_treatment = assigned.treatment,
      assigned_dist_group = dist.pot.group,
      signal = if_else(assigned_treatment %in% c("ink", "bracelet"), "signal", "no signal"),
      signal = factor(signal, levels = c("no signal", "signal")),
      dewormed_num = as.numeric(dewormed)
  ) %>%
  group_by(cluster.id) %>%
  mutate(
    clean_cluster_id = cur_group_id()
  ) %>%
  ungroup()

X_mat =  model.matrix(~ 0 + assigned_treatment*assigned_dist_group, data = cov_analysis_data)
N_signal_idx = cov_analysis_data$signal == "signal"


signal_Xs = rep(0, ncol(X_mat))
no_signal_Xs = rep(0, ncol(X_mat))
signal_idx = c(2, 4, 6, 8)
no_signal_idx =  c(1, 3, 5, 7)


stan_data = list(
    J = cov_analysis_data$clean_cluster_id %>% unique() %>% length(),
    N = nrow(cov_analysis_data),
    Y = cov_analysis_data$dewormed_num,
    X = X_mat,
    j_id = cov_analysis_data$clean_cluster_id,
    signal_idx = N_signal_idx,
    K = ncol(X_mat),
    signal_x_idx = signal_idx,
    no_signal_x_idx = no_signal_idx
)
options(mc.cores = 4)

model_fit = hier_model$sample(
    data = stan_data,
    chains = 4,
    iter_warmup = 400,
    iter_sampling = 400
)

colnames(X_mat) 


cov_name_df = tibble(
    coln = colnames(X_mat)
) %>%
    mutate(
        treatment = str_extract(coln, "(?<=assigned_treatment)\\w+"),
        treatment = if_else(coln == "assigned_dist_groupfar", "control", treatment),
        assigned_dist_group = if_else(str_detect(coln, "assigned_dist_group"), "far", "close"),
        k = 1:n()
    )

cov_name_df

gamma_rvar_df = model_fit %>%
    gather_rvars(
        mu_signal[i],
        mu_no_signal[k]
    ) 
    
    
    %>%
    bind_cols(
        cov_name_df
    ) %>%
    select(-coln)


fit_df = model_fit %>%
    gather_rvars(
        beta[j, k]
    ) %>%
    left_join(
        cov_name_df,
        by = "k"
    ) %>%
    select(-coln)

level_df = fit_df %>%
    group_by(
        treatment,
        assigned_dist_group
    ) %>%
    summarise(
        prob = rvar_mean(.value)
    ) %>%
    ungroup() %>%
    group_by(treatment) %>%
    mutate(
        prob = if_else(
            assigned_dist_group == "close",
            prob,
            prob + prob[assigned_dist_group == "close"]
        )   
    )

wide_level_df = level_df %>%
    pivot_wider(
        names_from = assigned_dist_group,
        values_from = prob
    ) %>%
    mutate(
        diff = far - close
    ) 


level_df %>%
    group_by(assigned_dist_group) %>%
    mutate(
        bra_m_cal = prob[treatment == "bracelet"] - prob[treatment == "calendar"]
    ) %>%
    select(bra_m_cal, assigned_dist_group) %>%
    slice(1) %>%
    pivot_wider(
        names_from = assigned_dist_group,
        values_from = bra_m_cal
    )   %>%
    mutate(
        diff = far - close
    ) %>%
    median_qi(diff)

level_df

wide_level_df %>%
    ungroup() %>%
    mutate(
        bra_m_cal = diff[treatment == "bracelet"] - diff[treatment == "calendar"]
    )

    median_qi(diff)

level_df %>%
    group_by(assigned_dist_group) %>%
    mutate(
        te = if_else(
            treatment == "control",
            prob,
            prob - prob[treatment == "control"]
        )
    ) %>%
    select(treatment, assigned_dist_group, te)  %>%
    pivot_wider(
        names_from = assigned_dist_group,
        values_from = te
    ) 

    summarise(
        treatment = "bracelet - calendar"
    )

level_df %>%
    pivot_wider(
        names_from = assigned_dist_group,
        values_from = prob
    ) %>%
    mutate(
        diff = far - close
    )
    
    bind_rows(
        .,
        . %>%
            mutate(

            )
    )
    mutate(
        bra_m_cal = 
    )


brm_fit = brm(
    bf(
    dewormed_num ~ 0 + (assigned_treatment*assigned_dist_group | signal) + county,
    sigma ~ cluster.id
    ),
    data = cov_analysis_data
)
