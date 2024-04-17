
  if (lnorm_wtp_value_utility_prior) {
    wtp_value_utility ~ lognormal(wtp_value_utility_mean, wtp_value_utility_sd);
  } else {
    wtp_value_utility ~ normal(wtp_value_utility_mean, wtp_value_utility_sd);
  }

  beta_intercept ~ normal(0, beta_intercept_sd);
  
  beta_ink_effect ~ normal(0, beta_ink_effect_sd);
  beta_calendar_effect ~ normal(0, beta_calendar_effect_sd);
  beta_bracelet_effect ~ normal(0, beta_bracelet_effect_sd);
  
  
  base_mu_rep ~ normal(0, mu_rep_sd);
  
  if (num_dist_param > 0) {
    dist_beta_v ~ normal(0, dist_beta_v_sd);
    
    if (use_param_dist_cluster_effects) {
      to_vector(dist_beta_cluster_raw) ~ normal(0, 1);
      dist_beta_cluster_sd ~ normal(0, 0.25);
    }
    
    if (use_param_dist_county_effects) {
      to_vector(dist_beta_county_raw) ~ normal(0, 1);
      dist_beta_county_sd ~ normal(0, 0.25);
    }
  }
  
  if (num_dist_param_quadratic > 0) {
    dist_quadratic_beta_v ~ normal(0, 1);
  }
  
  if (use_cluster_fixed_effect) {
    to_vector(structural_beta_cluster_raw) ~ std_normal();
    structural_beta_cluster_raw ~ std_normal();
    structural_beta_cluster_sd ~ normal(0, structural_beta_cluster_sd_sd);
  }
  
  if (use_county_effects) {
    structural_beta_county_raw ~ normal(0, 0.125);
    structural_beta_county_sd ~ normal(0, structural_beta_cluster_sd_sd);
  }

  raw_u_sd ~ inv_gamma(raw_u_sd_alpha, raw_u_sd_beta);
  raw_cluster_sd_tilde ~ std_normal();