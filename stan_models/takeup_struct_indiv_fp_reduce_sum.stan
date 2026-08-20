// Fit-105 individual fixed-point model specialized for its active branches.
//
// This preserves the original parameter block and likelihood but avoids
// constructing generic individual-by-treatment transformed-parameter arrays.
// The take-up likelihood solves fixed points in batched reduce_sum shards.

functions {
  #include struct_section_functions.stan

  real indiv_fp_partial_sum_lpmf(
      array[] int obs_slice,
      int start,
      int end,
      array[] int takeup,
      array[] int obs_cluster_id,
      vector cluster_standard_dist,
      array[] int cluster_incentive_treatment_id,
      array[] int cluster_assigned_dist_group_treatment,
      vector structural_treatment_effect,
      vector rep_intercept,
      vector rep_dist_slope,
      real base_mu_rep,
      real dist_beta,
      real total_error_sd,
      real u_sd,
      data int use_u_in_delta,
      data real alg_sol_rel_tol,
      data real alg_sol_f_tol,
      data real alg_sol_max_steps) {
    real lp = 0;

    for (slice_index in 1:size(obs_slice)) {
      int obs_index = obs_slice[slice_index];
      int cluster_index = obs_cluster_id[obs_index];
      int treatment_index =
        cluster_incentive_treatment_id[cluster_index];
      real distance = cluster_standard_dist[cluster_index];
      real mu_rep = base_mu_rep * inv_logit(
        rep_intercept[treatment_index] +
        rep_dist_slope[treatment_index] * distance
      );
      real benefit_cost =
        structural_treatment_effect[
          cluster_assigned_dist_group_treatment[cluster_index]
        ] -
        dist_beta * distance;
      real cutoff = find_fixedpoint_solution(
        benefit_cost,
        mu_rep,
        total_error_sd,
        u_sd,
        use_u_in_delta,
        alg_sol_rel_tol,
        alg_sol_f_tol,
        alg_sol_max_steps
      );
      real takeup_probability = Phi_approx(-cutoff / total_error_sd);

      lp += bernoulli_lpmf(takeup[obs_index] | takeup_probability);
    }

    return lp;
  }
}

data {
  #include struct_section_data.stan
}

transformed data {
  #include struct_section_transformed_data.stan

  if (
      use_binomial ||
      use_cost_model != COST_MODEL_TYPE_PARAM_LINEAR ||
      use_cluster_effects ||
      use_cluster_fixed_effect ||
      use_county_effects ||
      use_dist_cluster_effects ||
      use_dist_county_effects ||
      use_param_dist_cluster_effects ||
      use_param_dist_county_effects ||
      use_strata_levels ||
      beliefs_use_stratum_level ||
      beliefs_use_cluster_level ||
      beliefs_use_obs_level ||
      beliefs_use_indiv_intercept ||
      !beliefs_use_dist ||
      !use_homoskedastic_shocks ||
      !use_wtp_model ||
      suppress_reputation ||
      mu_rep_type != 4 ||
      BELIEFS_ORDER != 1 ||
      !lognormal_dist_model ||
      num_dist_group_mix != 1) {
    reject("Data are incompatible with the fit-105 individual-FP specialization.");
  }
}

parameters {
  #include struct_section_parameters.stan
}

model {
  array[num_discrete_dist] vector[num_dist_group_mix] group_dist_mean;
  vector[num_strata] strata_wtp_mu = rep_vector(hyper_wtp_mu, num_strata);
  vector[num_dist_group_treatments] beta = rep_vector(
    0,
    num_dist_group_treatments
  );
  vector[num_dist_group_treatments] structural_treatment_effect;
  vector[num_treatments] rep_intercept =
    beliefs_treatment_map_design_matrix * hyper_beta_1ord';
  vector[num_treatments] rep_dist_slope =
    beliefs_treatment_map_design_matrix * hyper_dist_beta_1ord';
  real u_sd = raw_u_sd[1];
  real total_error_sd = sqrt(1 + square(u_sd));

  for (dist_group_index in 1:num_discrete_dist) {
    group_dist_mean[dist_group_index] =
      hyper_dist_mean_mean +
      hyper_dist_mean_sd * group_dist_mean_raw[dist_group_index];
  }

  if (use_wtp_model) {
    beta[1:2] = [beta_intercept, beta_ink_effect]';
    beta[CALENDAR_TREATMENT_INDEX] =
      beta_bracelet_effect + wtp_value_utility * hyper_wtp_mu;
    beta[BRACELET_TREATMENT_INDEX] = beta_bracelet_effect;
  }
  structural_treatment_effect = treatment_map_design_matrix * beta;

  // Willingness-to-pay model: exact fit-105 branch with no stratum effects.
  #include wtp_model_section.stan

  // Beliefs model: exact global-coefficient branch.
  hyper_beta_1ord[1] ~ normal(0, 1);
  hyper_beta_1ord[2:] ~ normal(0, 0.5);
  hyper_dist_beta_1ord[1] ~ normal(0, 1);
  hyper_dist_beta_1ord[2:] ~ normal(0, 0.5);
  hyper_beta_2ord[1] ~ normal(0, 1);
  hyper_beta_2ord[2:] ~ normal(0, 0.5);
  hyper_dist_beta_2ord[1] ~ normal(0, 1);
  hyper_dist_beta_2ord[2:] ~ normal(0, 0.5);

  if (fit_beliefs_model_to_data) {
    vector[num_beliefs_obs] belief_distance =
      cluster_standard_dist[beliefs_cluster_index];
    vector[num_beliefs_obs] beliefs_latent_predictor_1ord =
      beliefs_treatment_design_matrix * hyper_beta_1ord' +
      (beliefs_treatment_design_matrix * hyper_dist_beta_1ord') .*
      belief_distance;
    vector[num_beliefs_obs] beliefs_latent_predictor_2ord =
      beliefs_treatment_design_matrix * hyper_beta_2ord' +
      (beliefs_treatment_design_matrix * hyper_dist_beta_2ord') .*
      belief_distance;

    for (belief_index in 1:num_beliefs_obs) {
      num_knows_1ord[belief_index] ~ binomial_logit(
        num_recognized[belief_index],
        beliefs_latent_predictor_1ord[belief_index]
      );
      num_knows_2ord[belief_index] ~ binomial_logit(
        num_recognized[belief_index],
        beliefs_latent_predictor_2ord[belief_index]
      );
    }
  }

  // Randomized-distance model: exact one-component/no-random-effects branch.
  for (dist_group_index in 1:num_discrete_dist) {
    if (lognormal_dist_model) {
      group_dist_mean_raw[dist_group_index, 1] ~ std_normal();
    } else {
      group_dist_mean_raw[dist_group_index, 1] ~ std_normal() T[0, ];
    }
    group_dist_sd[dist_group_index] ~ normal(0, hyper_dist_sd_sd);
  }
  if (fit_dist_model_to_data) {
    for (cluster_index in 1:num_clusters) {
      int dist_group = cluster_treatment_map[
        cluster_assigned_dist_group_treatment[cluster_index],
        2
      ];

      if (lognormal_dist_model) {
        cluster_standard_dist[cluster_index] ~ lognormal(
          group_dist_mean[dist_group, 1],
          group_dist_sd[dist_group, 1]
        );
      } else {
        target += normal_lpdf(
          cluster_standard_dist[cluster_index] |
          group_dist_mean[dist_group, 1],
          group_dist_sd[dist_group, 1]
        );
        target += normal_lccdf(
          0 |
          group_dist_mean[dist_group, 1],
          group_dist_sd[dist_group, 1]
        );
      }
    }
  }

  #include struct_model_priors.stan

  if (fit_model_to_data) {
    target += reduce_sum(
      indiv_fp_partial_sum_lpmf,
      included_monitored_obs,
      128,
      takeup,
      obs_cluster_id,
      cluster_standard_dist,
      cluster_incentive_treatment_id,
      cluster_assigned_dist_group_treatment,
      structural_treatment_effect,
      rep_intercept,
      rep_dist_slope,
      base_mu_rep,
      dist_beta_v[1],
      total_error_sd,
      u_sd,
      use_u_in_delta,
      alg_sol_rel_tol,
      alg_sol_f_tol,
      alg_sol_max_steps
    );
  }
}

generated quantities {
}
