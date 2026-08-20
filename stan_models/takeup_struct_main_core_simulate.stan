functions {
  #include struct_section_functions.stan
}

data {
  #include struct_section_data.stan
  vector<lower=0>[num_clusters] core_cluster_weight;
  int<lower=0, upper=1> use_core_cluster_shock;
  real<lower=0> core_cluster_shock_sd_prior;
  int<lower=0, upper=2> core_lambda_structure;
  real<lower=0> core_lambda_log_ratio_sd_prior;
  int<lower=0, upper=1> core_profile_group_lambda;
  real core_profile_group_log_ratio;
  int<lower=0, upper=1> core_gq_override_lambda;
  vector<lower=0>[num_treatments] core_gq_lambda_override;
  array[num_wtp_obs] int<lower=1, upper=num_clusters> wtp_cluster_id;

  int<lower=0, upper=1> core_sim_zero_lambda;
  vector[num_treatments] core_sim_lambda_log_deviation;
}

transformed data {
  #include struct_section_transformed_data.stan
}

parameters {
  #include struct_section_parameters.stan
  vector[use_core_cluster_shock ? num_clusters : 0]
    core_cluster_shock_raw;
  vector<lower=0>[use_core_cluster_shock ? 1 : 0]
    core_cluster_shock_sd;
  vector[core_lambda_structure == 1 && !core_profile_group_lambda ? 1 : 0]
    core_lambda_group_log_ratio_raw;
  vector[core_lambda_structure == 2 ? num_treatments - 1 : 0]
    core_lambda_arm_log_ratio_raw;
}

model {
}

generated quantities {
  array[num_obs] int<lower=0, upper=1> core_sim_takeup;
  array[num_beliefs_obs] int<lower=0> core_sim_num_knows_1ord;
  array[num_beliefs_obs] int<lower=0> core_sim_num_knows_2ord;
  array[num_wtp_obs] int<lower=-1, upper=1> core_sim_wtp_response;
  array[num_wtp_obs] int<lower=-1, upper=1> core_sim_gift_choice;
  vector[num_treatments] core_sim_signal_lambda;
  vector[num_clusters] core_sim_takeup_probability;

  {
    vector[num_dist_group_treatments] beta =
      rep_vector(0, num_dist_group_treatments);
    vector[num_dist_group_treatments] structural_treatment_effect;
    vector[num_treatments] rep_intercept =
      beliefs_treatment_map_design_matrix * hyper_beta_1ord';
    vector[num_treatments] rep_dist_slope =
      beliefs_treatment_map_design_matrix * hyper_dist_beta_1ord';
    vector[num_clusters] mu_rep;
    vector[num_clusters] benefit;
    real u_sd = raw_u_sd[1];
    real total_error_sd = sqrt(1 + square(u_sd));

    if (core_sim_zero_lambda) {
      core_sim_signal_lambda = rep_vector(0, num_treatments);
    } else {
      core_sim_signal_lambda =
        base_mu_rep * exp(core_sim_lambda_log_deviation);
    }
    beta[1:2] = [beta_intercept, beta_ink_effect]';
    beta[CALENDAR_TREATMENT_INDEX] =
      beta_bracelet_effect + wtp_value_utility * hyper_wtp_mu;
    beta[BRACELET_TREATMENT_INDEX] = beta_bracelet_effect;
    structural_treatment_effect = treatment_map_design_matrix * beta;
    mu_rep = core_sim_signal_lambda[cluster_incentive_treatment_id] .* inv_logit(
      rep_intercept[cluster_incentive_treatment_id] +
      rep_dist_slope[cluster_incentive_treatment_id] .* cluster_standard_dist
    );
    benefit =
      structural_treatment_effect[cluster_assigned_dist_group_treatment] -
      dist_beta_v[1] * cluster_standard_dist;
    for (cluster_index in 1:num_clusters) {
      real cutoff = find_fixedpoint_solution(
        benefit[cluster_index], mu_rep[cluster_index], total_error_sd, u_sd,
        use_u_in_delta, alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps
      );
      core_sim_takeup_probability[cluster_index] =
        Phi_approx(-cutoff / total_error_sd);
    }
    for (obs_index in 1:num_obs) {
      core_sim_takeup[obs_index] = bernoulli_rng(
        core_sim_takeup_probability[obs_cluster_id[obs_index]]
      );
    }

    for (belief_index in 1:num_beliefs_obs) {
      real distance = cluster_standard_dist[
        beliefs_cluster_index[belief_index]
      ];
      real predictor_1 =
        beliefs_treatment_design_matrix[belief_index] * hyper_beta_1ord' +
        beliefs_treatment_design_matrix[belief_index] *
        hyper_dist_beta_1ord' * distance;
      real predictor_2 =
        beliefs_treatment_design_matrix[belief_index] * hyper_beta_2ord' +
        beliefs_treatment_design_matrix[belief_index] *
        hyper_dist_beta_2ord' * distance;
      core_sim_num_knows_1ord[belief_index] = binomial_rng(
        num_recognized[belief_index], inv_logit(predictor_1)
      );
      core_sim_num_knows_2ord[belief_index] = binomial_rng(
        num_recognized[belief_index], inv_logit(predictor_2)
      );
    }

    for (wtp_index in 1:num_wtp_obs) {
      real latent_wtp = normal_rng(hyper_wtp_mu, wtp_sigma);
      real offer = scaled_wtp_offer[wtp_index];
      if (latent_wtp < -offer) {
        core_sim_wtp_response[wtp_index] = -1;
        core_sim_gift_choice[wtp_index] = -1;
      } else if (latent_wtp < 0) {
        core_sim_wtp_response[wtp_index] = 1;
        core_sim_gift_choice[wtp_index] = -1;
      } else if (latent_wtp <= offer) {
        core_sim_wtp_response[wtp_index] = 1;
        core_sim_gift_choice[wtp_index] = 1;
      } else {
        core_sim_wtp_response[wtp_index] = -1;
        core_sim_gift_choice[wtp_index] = 1;
      }
    }
  }
}
