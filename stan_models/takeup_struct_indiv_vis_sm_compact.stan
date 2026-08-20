functions {
  #include struct_section_functions.stan
}

data {
  #include struct_section_data.stan
  vector[num_obs] indiv_standard_dist;
}

transformed data {
  #include struct_section_transformed_data.stan
}

parameters {
  #include struct_section_parameters.stan
}

transformed parameters {
  #include struct_section_transformed_parameters.stan

  vector[num_obs] takeup_linear_predictor;
  vector[num_obs] indiv_dist_cost;
  vector[num_obs] takeup_pr;

  indiv_dist_cost = dist_beta_v[1] * indiv_standard_dist;
  takeup_linear_predictor =
    -structural_cluster_obs_v[obs_cluster_id[included_monitored_obs]] -
    indiv_dist_cost +
    cluster_dist_cost[obs_cluster_id[included_monitored_obs]];
  takeup_pr = Phi_approx(takeup_linear_predictor ./ total_error_sd[1]);
}

model {
  #include wtp_model_section.stan
  #include beliefs_model_sec.stan
  #include dist_model_sec.stan
  #include struct_model_priors.stan

  if (fit_model_to_data) {
    takeup[included_monitored_obs] ~ bernoulli(takeup_pr);
  }
}

generated quantities {
  #include takeup_struct_sm_compact_quantities.stan
}
