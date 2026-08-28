// Sampling-only entry point for the individual-travel-cost/community-fixed-point
// robustness model. This is algebraically identical to
// takeup_struct_private_info.stan, but deterministic observation-level arrays
// are model-local and therefore are not written for every retained draw.

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

model {
  {
    #include struct_section_model_locals.stan

    vector[num_obs] indiv_dist_cost =
      dist_beta_v[1] * indiv_standard_dist;
    vector[num_obs] takeup_linear_predictor =
      -structural_cluster_obs_v[obs_cluster_id[included_monitored_obs]] -
      indiv_dist_cost +
      cluster_dist_cost[obs_cluster_id[included_monitored_obs]];
    vector[num_obs] takeup_pr =
      Phi_approx(takeup_linear_predictor ./ total_error_sd[1]);

    #include wtp_model_section.stan
    #include beliefs_model_sec.stan
    #include dist_model_sec.stan
    #include struct_model_priors.stan

    if (fit_model_to_data) {
      takeup[included_monitored_obs] ~ bernoulli(takeup_pr);
    }
  }
}

generated quantities {
}
