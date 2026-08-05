// Fit-105 main community fixed-point model specialized for the paper's active
// branches. Only sampled parameters are written; all deterministic quantities
// are model-local and can be regenerated after sampling with a GQ program.

functions {
  #include struct_section_functions.stan
  // Heavy-tailed type robustness with bounded Newton solves; inactive in the
  // exact Gaussian branch.
  #include core_student_t_functions.stan
  #include core_asymmetric_observability_functions.stan
}

data {
  #include struct_section_data.stan

  // Minimal-model extensions. Unit weights and a disabled shock recover the
  // original fit-105 target exactly.
  vector<lower=0>[num_clusters] core_cluster_weight;
  int<lower=0, upper=1> use_core_cluster_shock;
  real<lower=0> core_cluster_shock_sd_prior;
  // 0 = common; 1 = Any Signal (Ink/Bracelet) vs No Signal
  // (Control/Calendar); 2 = arm-specific.
  int<lower=0, upper=2> core_lambda_structure;
  // In both flexible structures this is the prior SD of a pairwise
  // log-lambda ratio.
  real<lower=0> core_lambda_log_ratio_sd_prior;
  int<lower=0, upper=1> core_profile_group_lambda;
  real core_profile_group_log_ratio;
  int<lower=0, upper=1> core_gq_override_lambda;
  vector<lower=0>[num_treatments] core_gq_lambda_override;
  int<lower=0, upper=1> core_type_distribution;
  real<lower=2> core_student_t_df;
  real<lower=0> core_type_scale_sq;
  int<lower=2> core_type_mixture_components;
  vector<lower=0>[core_type_mixture_components]
    core_type_mixture_precision;
  simplex[core_type_mixture_components] core_type_mixture_weight;
  array[num_wtp_obs] int<lower=1, upper=num_clusters> wtp_cluster_id;
  // 0 = baseline perfect-observation schedule; 1 = asymmetric reports
  // conditional on recognition; 2 = unrecognized is a fourth null signal.
  int<lower=0, upper=2> core_observation_model;
  // Recognition: 0 = full; 1 = truth-specific intercepts; 2 = condition out.
  int<lower=0, upper=2> core_recognition_structure;
  // Reports: 0 = full multinomial; 1 = no arm-distance slopes; 2 = two-stage.
  int<lower=0, upper=2> core_report_structure;
  // Full-multinomial arm-distance slopes: 0 = independent; 1 = hierarchical.
  int<lower=0, upper=1> core_report_arm_dist_hierarchical;
  real<lower=0> core_report_arm_dist_prior_scale;
  int<lower=0> core_num_peer_response_rows;
  array[core_num_peer_response_rows] int<lower=1, upper=num_clusters>
    core_peer_response_cluster_id;
  array[core_num_peer_response_rows] int<lower=0, upper=1>
    core_peer_true_takeup;
  array[core_num_peer_response_rows] int<lower=0> core_peer_total;
  array[core_num_peer_response_rows] int<lower=0> core_peer_recognized;
  array[core_num_peer_response_rows, 3] int<lower=0>
    core_peer_report_count;
}

transformed data {
  #include struct_section_transformed_data.stan

  array[num_clusters] int monitored_cluster_size = rep_array(0, num_clusters);
  array[num_clusters] int monitored_cluster_takeup = rep_array(0, num_clusters);
  int core_has_nonunit_weight = 0;
  array[num_treatments] int core_is_public_signal = rep_array(0, num_treatments);
  matrix[num_treatments, num_treatments - 1]
    core_signal_lambda_contrast_basis = rep_matrix(
      0,
      num_treatments,
      num_treatments - 1
    );

  // Orthonormal Helmert basis for treatment-specific log-lambda deviations.
  // Every column sums to zero, so base_mu_rep remains their geometric mean.
  for (contrast_index in 1:(num_treatments - 1)) {
    real denominator = sqrt(contrast_index * (contrast_index + 1.0));
    for (treatment_index in 1:contrast_index) {
      core_signal_lambda_contrast_basis[treatment_index, contrast_index] =
        1 / denominator;
    }
    core_signal_lambda_contrast_basis[contrast_index + 1, contrast_index] =
      -contrast_index / denominator;
  }
  core_is_public_signal[2] = 1;
  core_is_public_signal[BRACELET_TREATMENT_INDEX] = 1;

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
  if (num_treatments != 4 || CALENDAR_TREATMENT_INDEX != 3 ||
      BRACELET_TREATMENT_INDEX != 4) {
    reject("Core lambda grouping requires Control/Ink/Calendar/Bracelet order.");
  }
  if (core_profile_group_lambda && core_lambda_structure != 1) {
    reject("Lambda profiling is available only for the grouped structure.");
  }
  if (core_observation_model > 0 &&
      (core_lambda_structure != 0 || core_type_distribution != 0)) {
    reject("Asymmetric observability currently requires common lambda and Gaussian types.");
  }
  if (core_observation_model > 0 && core_num_peer_response_rows == 0) {
    reject("Asymmetric observability requires linked peer-response rows.");
  }
  if (core_observation_model == 0 &&
      (core_recognition_structure != 0 || core_report_structure != 0)) {
    reject("Observation restrictions require an asymmetric observation model.");
  }
  if (core_observation_model == 2 && core_recognition_structure == 2) {
    reject("The unconditional channel cannot condition recognition out.");
  }
  if (core_report_arm_dist_hierarchical &&
      (core_observation_model == 0 || core_report_structure != 0)) {
    reject("Hierarchical report-distance slopes require the full multinomial channel.");
  }
  if (core_num_peer_response_rows > 0) {
    for (row in 1:core_num_peer_response_rows) {
      if (core_peer_recognized[row] > core_peer_total[row] ||
          sum(core_peer_report_count[row]) != core_peer_recognized[row]) {
        reject("Invalid linked peer-response sufficient statistics.");
      }
    }
  }
  if (core_type_distribution == 1 &&
      abs(core_type_scale_sq - (core_student_t_df - 2) /
          core_student_t_df) > 1e-8) {
    reject("Student-t type scale must normalize latent-type variance to one.");
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
  vector[core_lambda_structure == 1 && !core_profile_group_lambda ? 1 : 0]
    core_lambda_group_log_ratio_raw;
  vector[core_lambda_structure == 2 ? num_treatments - 1 : 0]
    core_lambda_arm_log_ratio_raw;
  vector[core_observation_model > 0 && core_recognition_structure < 2 ? 2 : 0]
    core_recognition_intercept;
  vector[core_observation_model > 0 && core_recognition_structure == 0 ? 2 : 0]
    core_recognition_dist_slope;
  matrix[core_observation_model > 0 && core_recognition_structure == 0 ? 2 : 0, num_treatments - 1]
    core_recognition_arm_intercept_raw;
  matrix[core_observation_model > 0 && core_recognition_structure == 0 ? 2 : 0, num_treatments - 1]
    core_recognition_arm_dist_raw;
  matrix[core_observation_model > 0 && core_report_structure < 2 ? 2 : 0, 2]
    core_report_intercept;
  matrix[core_observation_model > 0 && core_report_structure < 2 ? 2 : 0, 2]
    core_report_dist_slope;
  matrix[core_observation_model > 0 && core_report_structure < 2 ? 2 : 0, 2 * (num_treatments - 1)]
    core_report_arm_intercept_raw;
  matrix[core_observation_model > 0 && core_report_structure == 0 ? 2 : 0, 2 * (num_treatments - 1)]
    core_report_arm_dist_raw;
  matrix<lower=0>[core_observation_model > 0 && core_report_structure == 0 && core_report_arm_dist_hierarchical ? 2 : 0, 2]
    core_report_arm_dist_sd;
  vector[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0]
    core_definite_intercept;
  vector[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0]
    core_definite_dist_slope;
  matrix[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0, num_treatments - 1]
    core_definite_arm_intercept_raw;
  vector[core_observation_model > 0 && core_report_structure == 2 ? 1 : 0]
    core_definite_public_signal_dist_slope;
  vector[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0]
    core_accuracy_intercept;
  matrix[core_observation_model > 0 && core_report_structure == 2 ? 2 : 0, num_treatments - 1]
    core_accuracy_arm_intercept_raw;
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
  vector[num_treatments] signal_lambda =
    rep_vector(base_mu_rep, num_treatments);
  vector[num_clusters] cluster_benefit_cost;
  vector[num_clusters] core_cluster_shock = rep_vector(0, num_clusters);
  vector[num_clusters] cluster_cutoff;
  vector[num_clusters] cluster_takeup_probability;
  matrix[num_clusters, 2] core_cluster_recognition = rep_matrix(0, num_clusters, 2);
  matrix[num_clusters, 3] core_cluster_report_taker = rep_matrix(0, num_clusters, 3);
  matrix[num_clusters, 3] core_cluster_report_nontaker = rep_matrix(0, num_clusters, 3);
  matrix[num_clusters, 4] core_cluster_signal_taker = rep_matrix(0, num_clusters, 4);
  matrix[num_clusters, 4] core_cluster_signal_nontaker = rep_matrix(0, num_clusters, 4);
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

  if (core_lambda_structure == 1) {
    real group_log_ratio = core_profile_group_lambda ?
      core_profile_group_log_ratio :
      core_lambda_log_ratio_sd_prior * core_lambda_group_log_ratio_raw[1];
    for (treatment_index in 1:num_treatments) {
      signal_lambda[treatment_index] = base_mu_rep * exp(
        (core_is_public_signal[treatment_index] - 0.5) * group_log_ratio
      );
    }
    if (!core_profile_group_lambda) {
      core_lambda_group_log_ratio_raw ~ std_normal();
    }
  } else if (core_lambda_structure == 2) {
    signal_lambda = base_mu_rep * exp(
      core_lambda_log_ratio_sd_prior / sqrt(2.0) *
      core_signal_lambda_contrast_basis *
      core_lambda_arm_log_ratio_raw
    );
    core_lambda_arm_log_ratio_raw ~ std_normal();
  }

  cluster_mu_rep = signal_lambda[cluster_incentive_treatment_id] .* inv_logit(
    rep_intercept[cluster_incentive_treatment_id] +
    rep_dist_slope[cluster_incentive_treatment_id] .* cluster_standard_dist
  );
  if (core_observation_model > 0) {
    core_recognition_intercept ~ normal(0, 1.5);
    core_recognition_dist_slope ~ normal(0, 0.5);
    to_vector(core_recognition_arm_intercept_raw) ~ normal(0, 0.5);
    to_vector(core_recognition_arm_dist_raw) ~ normal(0, 0.25);
    to_vector(core_report_intercept) ~ normal(0, 1.5);
    to_vector(core_report_dist_slope) ~ normal(0, 0.5);
    to_vector(core_report_arm_intercept_raw) ~ normal(0, 0.5);
    if (core_report_arm_dist_hierarchical) {
      to_vector(core_report_arm_dist_raw) ~ std_normal();
      to_vector(core_report_arm_dist_sd) ~ normal(
        0, core_report_arm_dist_prior_scale
      );
    } else {
      to_vector(core_report_arm_dist_raw) ~ normal(
        0, core_report_arm_dist_prior_scale
      );
    }
    core_definite_intercept ~ normal(0, 1.5);
    core_definite_dist_slope ~ normal(0, 0.5);
    to_vector(core_definite_arm_intercept_raw) ~ normal(0, 0.5);
    core_definite_public_signal_dist_slope ~ normal(0, 0.25);
    core_accuracy_intercept ~ normal(0, 1.5);
    to_vector(core_accuracy_arm_intercept_raw) ~ normal(0, 0.5);

    for (cluster in 1:num_clusters) {
      int treatment = cluster_incentive_treatment_id[cluster];
      for (truth in 1:2) {
        vector[3] conditional_report;
        real recognition_prob = 1;
        if (core_recognition_structure < 2) {
          real recognition_eta = core_recognition_intercept[truth];
          if (core_recognition_structure == 0) {
            recognition_eta += dot_product(
              core_recognition_arm_intercept_raw[truth],
              core_signal_lambda_contrast_basis[treatment]
            ) + (
              core_recognition_dist_slope[truth] + dot_product(
                core_recognition_arm_dist_raw[truth],
                core_signal_lambda_contrast_basis[treatment]
              )
            ) * cluster_standard_dist[cluster];
          }
          recognition_prob = inv_logit(recognition_eta);
        }
        if (core_report_structure == 2) {
          real definite_eta = core_definite_intercept[truth] + dot_product(
            core_definite_arm_intercept_raw[truth],
            core_signal_lambda_contrast_basis[treatment]
          ) + (
            core_definite_dist_slope[truth] +
            core_definite_public_signal_dist_slope[1] *
              core_is_public_signal[treatment]
          ) * cluster_standard_dist[cluster];
          real accuracy_eta = core_accuracy_intercept[truth] + dot_product(
            core_accuracy_arm_intercept_raw[truth],
            core_signal_lambda_contrast_basis[treatment]
          );
          conditional_report = core_two_stage_report_row(
            inv_logit(definite_eta), inv_logit(accuracy_eta), truth
          );
        } else {
          vector[2] report_logit;
          for (category in 1:2) {
            int first = (category - 1) * (num_treatments - 1) + 1;
            int last = category * (num_treatments - 1);
            real report_arm_slope = 0;
            if (core_report_structure == 0) {
              report_arm_slope = dot_product(
                core_report_arm_dist_raw[truth, first:last],
                core_signal_lambda_contrast_basis[treatment]
              );
              if (core_report_arm_dist_hierarchical) {
                report_arm_slope *= core_report_arm_dist_sd[truth, category];
              }
            }
            report_logit[category] = core_report_intercept[truth, category] +
              dot_product(
                core_report_arm_intercept_raw[truth, first:last],
                core_signal_lambda_contrast_basis[treatment]
              ) + (core_report_dist_slope[truth, category] +
                report_arm_slope) * cluster_standard_dist[cluster];
          }
          conditional_report = core_softmax_with_reference(
            report_logit[1], report_logit[2]
          );
        }
        core_cluster_recognition[cluster, truth] = recognition_prob;
        if (truth == 1) {
          core_cluster_report_nontaker[cluster] = conditional_report';
          core_cluster_signal_nontaker[cluster] = core_noisy_signal_row(
            recognition_prob, conditional_report, core_observation_model
          )';
        } else {
          core_cluster_report_taker[cluster] = conditional_report';
          core_cluster_signal_taker[cluster] = core_noisy_signal_row(
            recognition_prob, conditional_report, core_observation_model
          )';
        }
      }
    }
  }
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

  if (fit_beliefs_model_to_data && core_observation_model == 0) {
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
  } else if (fit_beliefs_model_to_data) {
    for (row in 1:core_num_peer_response_rows) {
      int cluster = core_peer_response_cluster_id[row];
      int truth = core_peer_true_takeup[row] + 1;
      vector[3] report_probability = truth == 1 ?
        core_cluster_report_nontaker[cluster]' :
        core_cluster_report_taker[cluster]';
      if (core_recognition_structure < 2) {
        target += core_cluster_weight[cluster] * binomial_lupmf(
          core_peer_recognized[row] |
          core_peer_total[row],
          core_cluster_recognition[cluster, truth]
        );
      }
      target += core_cluster_weight[cluster] * multinomial_lupmf(
        core_peer_report_count[row] | report_probability
      );
    }
    // Retain the second-order belief likelihood; it is not used as the
    // first-order social-information technology in these modes.
    for (belief_index in 1:num_beliefs_obs) {
      int cluster = beliefs_cluster_index[belief_index];
      real predictor = dot_product(
        beliefs_treatment_design_matrix[belief_index], hyper_beta_2ord
      ) + dot_product(
        beliefs_treatment_design_matrix[belief_index], hyper_dist_beta_2ord
      ) * cluster_standard_dist[cluster];
      target += core_cluster_weight[cluster] * binomial_logit_lupmf(
        num_knows_2ord[belief_index] | num_recognized[belief_index], predictor
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
    if (core_observation_model > 0) {
      int num_signals = core_observation_model == 1 ? 3 : 4;
      if (multithreaded) {
        cluster_cutoff = core_noisy_map_find_fixedpoint(
          cluster_benefit_cost,
          base_mu_rep,
          rep_vector(total_error_sd, num_clusters),
          rep_vector(u_sd, num_clusters),
          core_cluster_signal_taker,
          core_cluster_signal_nontaker,
          num_signals,
          use_u_in_delta,
          alg_sol_rel_tol,
          alg_sol_f_tol,
          alg_sol_max_steps
        );
      } else {
        for (cluster in 1:num_clusters) {
          cluster_cutoff[cluster] = core_noisy_find_fixedpoint(
            cluster_benefit_cost[cluster], base_mu_rep, total_error_sd, u_sd,
            core_cluster_signal_taker[cluster]',
            core_cluster_signal_nontaker[cluster]', num_signals,
            use_u_in_delta, alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps
          );
        }
      }
      cluster_takeup_probability = Phi_approx(
        -cluster_cutoff / total_error_sd
      );
    } else if (core_type_distribution == 1) {
      if (multithreaded) {
        cluster_cutoff = core_student_t_map_find_fixedpoint(
          cluster_benefit_cost,
          cluster_mu_rep,
          u_sd,
          use_u_in_delta,
          core_type_scale_sq,
          core_type_mixture_precision,
          core_type_mixture_weight,
          alg_sol_rel_tol,
          alg_sol_f_tol,
          alg_sol_max_steps
        );
      } else {
        for (cluster_index in 1:num_clusters) {
          cluster_cutoff[cluster_index] = core_student_t_find_fixedpoint(
            cluster_benefit_cost[cluster_index],
            cluster_mu_rep[cluster_index],
            u_sd,
            use_u_in_delta,
            core_type_scale_sq,
            core_type_mixture_precision,
            core_type_mixture_weight,
            alg_sol_rel_tol,
            alg_sol_f_tol,
            alg_sol_max_steps
          );
        }
      }
      for (cluster_index in 1:num_clusters) {
        cluster_takeup_probability[cluster_index] = 1 -
          core_student_t_moments(
            cluster_cutoff[cluster_index], u_sd, use_u_in_delta,
            core_type_scale_sq, core_type_mixture_precision,
            core_type_mixture_weight
          )[1];
      }
    } else {
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
    }

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
