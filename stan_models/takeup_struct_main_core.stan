// Fit-105 main community fixed-point model specialized for the paper's active
// branches. Only sampled parameters are written; all deterministic quantities
// are model-local and can be regenerated after sampling with a GQ program.

functions {
  #include struct_section_functions.stan
}

data {
  #include struct_section_data.stan

  // Minimal-model extensions. Unit weights and a disabled shock recover the
  // original fit-105 target exactly.
  vector<lower=0>[num_clusters] core_cluster_weight;
  int<lower=0, upper=1> use_core_cluster_shock;
  real<lower=0> core_cluster_shock_sd_prior;
  array[num_wtp_obs] int<lower=1, upper=num_clusters> wtp_cluster_id;
}

transformed data {
  #include struct_section_transformed_data.stan

  array[num_clusters] int monitored_cluster_size = rep_array(0, num_clusters);
  array[num_clusters] int monitored_cluster_takeup = rep_array(0, num_clusters);
  int core_has_nonunit_weight = 0;

  for (cluster_index in 1:num_clusters) {
    if (abs(core_cluster_weight[cluster_index] - 1) > 1e-12) {
      core_has_nonunit_weight = 1;
    }
  }

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
    reject("Data are incompatible with the fit-105 main-model specialization.");
  }

  // The original likelihood is a product of Bernoulli terms with a common
  // probability within cluster. In the model block, binomial_lupmf explicitly
  // omits the data-only combinatorial term, giving the identical Stan target.
  for (included_index in 1:num_included_monitored_obs) {
    int obs_index = included_monitored_obs[included_index];
    int cluster_index = obs_cluster_id[obs_index];
    monitored_cluster_size[cluster_index] += 1;
    monitored_cluster_takeup[cluster_index] += takeup[obs_index];
  }
}

parameters {
  // Keep the canonical parameter schema so draws remain usable by the
  // existing post-processing and generated-quantities programs.
  #include struct_section_parameters.stan
  vector[use_core_cluster_shock ? num_clusters : 0]
    core_cluster_shock_raw;
  vector<lower=0>[use_core_cluster_shock ? 1 : 0]
    core_cluster_shock_sd;
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
  vector[num_clusters] cluster_mu_rep;
  vector[num_clusters] cluster_benefit_cost;
  vector[num_clusters] core_cluster_shock = rep_vector(0, num_clusters);
  vector[num_clusters] cluster_cutoff;
  vector[num_clusters] cluster_takeup_probability;
  real u_sd = raw_u_sd[1];
  real total_error_sd = sqrt(1 + square(u_sd));

  for (dist_group_index in 1:num_discrete_dist) {
    group_dist_mean[dist_group_index] =
      hyper_dist_mean_mean +
      hyper_dist_mean_sd * group_dist_mean_raw[dist_group_index];
  }

  beta[1:2] = [beta_intercept, beta_ink_effect]';
  beta[CALENDAR_TREATMENT_INDEX] =
    beta_bracelet_effect + wtp_value_utility * hyper_wtp_mu;
  beta[BRACELET_TREATMENT_INDEX] = beta_bracelet_effect;
  structural_treatment_effect = treatment_map_design_matrix * beta;

  cluster_mu_rep = base_mu_rep * inv_logit(
    rep_intercept[cluster_incentive_treatment_id] +
    rep_dist_slope[cluster_incentive_treatment_id] .* cluster_standard_dist
  );
  cluster_benefit_cost =
    structural_treatment_effect[cluster_assigned_dist_group_treatment] -
    dist_beta_v[1] * cluster_standard_dist;
  if (use_core_cluster_shock) {
    core_cluster_shock =
      core_cluster_shock_sd[1] * core_cluster_shock_raw;
    cluster_benefit_cost += core_cluster_shock;
    core_cluster_shock_raw ~ std_normal();
    core_cluster_shock_sd ~ normal(0, core_cluster_shock_sd_prior);
  }

  // Willingness-to-pay model: exact fit-105 branch with no stratum effects.
  #include wtp_model_core_weighted.stan

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
      int cluster_index = beliefs_cluster_index[belief_index];
      target += core_cluster_weight[cluster_index] * binomial_logit_lupmf(
        num_knows_1ord[belief_index] |
        num_recognized[belief_index],
        beliefs_latent_predictor_1ord[belief_index]
      );
      target += core_cluster_weight[cluster_index] * binomial_logit_lupmf(
        num_knows_2ord[belief_index] |
        num_recognized[belief_index],
        beliefs_latent_predictor_2ord[belief_index]
      );
    }
  }

  // Randomized-distance model: exact one-component/no-random-effects branch.
  for (dist_group_index in 1:num_discrete_dist) {
    group_dist_mean_raw[dist_group_index, 1] ~ std_normal();
    group_dist_sd[dist_group_index] ~ normal(0, hyper_dist_sd_sd);
  }
  if (fit_dist_model_to_data) {
    for (cluster_index in 1:num_clusters) {
      int dist_group = cluster_treatment_map[
        cluster_assigned_dist_group_treatment[cluster_index],
        2
      ];
      target += core_cluster_weight[cluster_index] * lognormal_lupdf(
        cluster_standard_dist[cluster_index] |
        group_dist_mean[dist_group, 1],
        group_dist_sd[dist_group, 1]
      );
    }
  }

  #include struct_model_priors.stan

  if (fit_model_to_data) {
    if (multithreaded) {
      cluster_cutoff = map_find_fixedpoint_solution(
        cluster_benefit_cost,
        cluster_mu_rep,
        rep_vector(total_error_sd, num_clusters),
        rep_vector(u_sd, num_clusters),
        use_u_in_delta,
        alg_sol_rel_tol,
        alg_sol_f_tol,
        alg_sol_max_steps
      );
    } else {
      for (cluster_index in 1:num_clusters) {
        cluster_cutoff[cluster_index] = find_fixedpoint_solution(
          cluster_benefit_cost[cluster_index],
          cluster_mu_rep[cluster_index],
          total_error_sd,
          u_sd,
          use_u_in_delta,
          alg_sol_rel_tol,
          alg_sol_f_tol,
          alg_sol_max_steps
        );
      }
    }
    cluster_takeup_probability = Phi_approx(
      -cluster_cutoff / total_error_sd
    );

    for (cluster_index in included_clusters) {
      target += core_cluster_weight[cluster_index] * binomial_lupmf(
        monitored_cluster_takeup[cluster_index] |
        monitored_cluster_size[cluster_index],
        cluster_takeup_probability[cluster_index]
      );
    }
  }
}

generated quantities {
}
