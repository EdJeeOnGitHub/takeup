functions {
  #include struct_section_functions.stan
}

data {
  #include struct_section_data.stan
  // Private vs Community Knowledge
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
    // Additional private/community info stuff
    vector[num_obs] indiv_mu_rep;
    vector[num_obs] indiv_dist_cost;
    vector[num_obs] takeup_linear_predictor;
    vector[num_obs] takeup_pr;
    vector[num_obs] indiv_delta_vstar;

    // HH's V* is the community centroid cost (i.e. V* calculated at cluster level above)
    indiv_delta_vstar = expected_delta(structural_cluster_obs_v[obs_cluster_id[included_monitored_obs]], total_error_sd[1], u_sd[1]);

    indiv_mu_rep = calculate_mu_rep(
        // treatment ids; array[] int
        cluster_incentive_treatment_id[obs_cluster_id[included_monitored_obs]], 
        // distance; vector
        indiv_standard_dist,
        // base mu rep; real
        base_mu_rep, 
        // mu_beliefs_effect; real
        1, 
        // design matrix; matrix
        beliefs_treatment_map_design_matrix, 
        // beta; matrix
        centered_cluster_beta_beliefs[obs_cluster_id[included_monitored_obs], ], 
        // dist_beta; matrix
        centered_cluster_dist_beta_beliefs[obs_cluster_id[included_monitored_obs],],
        // mu rep type; int
        mu_rep_type
    );

    indiv_dist_cost = dist_beta_v[1] * indiv_standard_dist;
    
    takeup_linear_predictor = structural_treatment_effect[cluster_assigned_dist_group_treatment[obs_cluster_id]] - 
                                indiv_dist_cost + 
                                indiv_mu_rep .* indiv_delta_vstar;
    // total error sd doesn't vary with treatment so just use first one
    takeup_pr = Phi_approx(takeup_linear_predictor ./ total_error_sd[1]);
}

model {
  // Model section adjusted for private/community info
#include wtp_model_section.stan
#include beliefs_model_sec.stan
#include dist_model_sec.stan
#include struct_model_priors.stan


  if (fit_model_to_data) {
    takeup[included_monitored_obs] ~ bernoulli(takeup_pr);
  }
}

generated quantities {
  // no generated quantities again
}
