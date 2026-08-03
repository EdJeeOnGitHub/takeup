// Compact generated quantities for takeup_struct_main_core.stan. The data and
// parameter blocks must remain schema-identical to the sampling model.

functions {
  #include struct_section_functions.stan

  real core_gq_fixedpoint_residual(
      real cutoff,
      real benefit_cost,
      real mu_rep,
      real total_error_sd,
      real u_sd,
      data int use_u_in_delta) {
    array[1] int x_i = {use_u_in_delta};
    real delta;
    if (use_u_in_delta && u_sd > 0) {
      delta = expected_delta(
        cutoff, total_error_sd, u_sd, {0.0}, x_i
      );
    } else {
      delta = reputational_returns_normal(cutoff);
    }
    return cutoff + benefit_cost + mu_rep * delta;
  }

  // A bounded scalar bisection avoids sporadic algebra-solver failures in
  // extreme counterfactual draws. The paper's stable-equilibrium restriction
  // gives a unique crossing in this range for the retained posterior draws.
  real core_gq_find_fixedpoint_safe(
      real benefit_cost,
      real mu_rep,
      real total_error_sd,
      real u_sd,
      data int use_u_in_delta,
      data real alg_sol_rel_tol,
      data real alg_sol_f_tol,
      data real alg_sol_max_steps) {
    real lower_bound = -5.5;
    real upper_bound = 5.5;
    real lower_residual = core_gq_fixedpoint_residual(
      lower_bound, benefit_cost, mu_rep, total_error_sd, u_sd, use_u_in_delta
    );
    real upper_residual = core_gq_fixedpoint_residual(
      upper_bound, benefit_cost, mu_rep, total_error_sd, u_sd, use_u_in_delta
    );
    real midpoint = -benefit_cost;

    if (lower_residual * upper_residual <= 0) {
      for (step in 1:50) {
        real midpoint_residual;
        midpoint = 0.5 * (lower_bound + upper_bound);
        midpoint_residual = core_gq_fixedpoint_residual(
          midpoint, benefit_cost, mu_rep, total_error_sd, u_sd,
          use_u_in_delta
        );
        if (lower_residual * midpoint_residual <= 0) {
          upper_bound = midpoint;
          upper_residual = midpoint_residual;
        } else {
          lower_bound = midpoint;
          lower_residual = midpoint_residual;
        }
      }
    } else {
      // Defensive fallback for an unstable/no-crossing draw. Damped fixed-point
      // iteration stays finite, allowing downstream finite-draw checks to flag
      // it rather than losing the entire generated-quantities chain.
      for (step in 1:100) {
        real residual = core_gq_fixedpoint_residual(
          midpoint, benefit_cost, mu_rep, total_error_sd, u_sd,
          use_u_in_delta
        );
        midpoint = fmin(5.5, fmax(-5.5, midpoint - 0.5 * residual));
      }
    }
    return midpoint;
  }

  vector core_gq_calculate_roc(
      real benefit_cost,
      real benefit_cost_control,
      real total_error_sd,
      real u_sd,
      real dist_beta,
      real mu_rep,
      real mu_rep_control,
      real mu_rep_deriv,
      data int use_u_in_delta,
      data real alg_sol_rel_tol,
      data real alg_sol_f_tol,
      data real alg_sol_max_steps) {
    array[1] int x_i = {use_u_in_delta};
    real w = core_gq_find_fixedpoint_safe(
      benefit_cost, mu_rep, total_error_sd, u_sd, use_u_in_delta,
      alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps
    );
    real w_control = core_gq_find_fixedpoint_safe(
      benefit_cost_control, mu_rep_control, total_error_sd, u_sd,
      use_u_in_delta, alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps
    );
    vector[2] delta_not_control = expected_delta_deriv(
      w, total_error_sd, u_sd, {0.0}, x_i
    );
    vector[2] delta_control = expected_delta_deriv(
      w_control, total_error_sd, u_sd, {0.0}, x_i
    );
    real roc = calculate_roc(
      w_control, total_error_sd, dist_beta, mu_rep, mu_rep_deriv,
      delta_control[1], delta_control[2]
    );
    real roc_no_visibility =
      -exp(normal_lpdf(w_control | 0, total_error_sd)) * dist_beta;
    real sm_value = calculate_sm(
      w, total_error_sd, dist_beta, mu_rep, mu_rep_deriv,
      delta_not_control[1], delta_not_control[2]
    );
    return [
      w, w_control, delta_control[1], delta_control[2], roc,
      roc_no_visibility, sm_value,
      delta_not_control[1], delta_not_control[2]
    ]';
  }
}

data {
  #include struct_section_data.stan
  vector<lower=0>[num_clusters] core_cluster_weight;
  int<lower=0, upper=1> use_core_cluster_shock;
  real<lower=0> core_cluster_shock_sd_prior;
  array[num_wtp_obs] int<lower=1, upper=num_clusters> wtp_cluster_id;
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
}

model {
}

generated quantities {
  // Rows: Combined, Close, Far. Columns follow the canonical Stan treatment
  // order and are reordered only by paper postprocessing.
  matrix[3, num_treatments] core_compact_takeup_level;
  // Rows: 500m, 1500m, 2500m.
  matrix[3, num_treatments] core_compact_sm;
  matrix[3, num_treatments] core_compact_sm_rescaled;
  real core_compact_cluster_shock_sd = 0;

  {
    array[3] int compact_sm_roc_index = { 6, 16, 26 };
    array[num_discrete_dist] vector[num_dist_group_mix] group_dist_mean;
    vector[num_dist_group_treatments] beta =
      rep_vector(0, num_dist_group_treatments);
    vector[num_dist_group_treatments] structural_treatment_effect;
    vector[num_treatments] rep_intercept =
      beliefs_treatment_map_design_matrix * hyper_beta_1ord';
    vector[num_treatments] rep_dist_slope =
      beliefs_treatment_map_design_matrix * hyper_dist_beta_1ord';
    vector[num_clusters] cluster_shock = rep_vector(0, num_clusters);
    vector[num_treatments] combined_numerator =
      rep_vector(0, num_treatments);
    real total_weight = 0;
    real u_sd = raw_u_sd[1];
    real total_error_sd = sqrt(1 + square(u_sd));
    matrix[num_clusters, num_discrete_dist] all_cluster_dist;

    if (use_core_cluster_shock) {
      core_compact_cluster_shock_sd = core_cluster_shock_sd[1];
      cluster_shock = core_cluster_shock_sd[1] * core_cluster_shock_raw;
    }

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

    // Match the paper's historical counterfactual-distance construction.
    for (cluster_index in 1:num_clusters) {
      int observed_dist_group = cluster_dist_treatment_id[cluster_index];
      real size_weight = sum(cluster_size[cluster_index]);
      total_weight += size_weight;
      for (dist_group_index in 1:num_discrete_dist) {
        if (dist_group_index == observed_dist_group) {
          all_cluster_dist[cluster_index, dist_group_index] =
            cluster_standard_dist[cluster_index];
        } else {
          all_cluster_dist[cluster_index, dist_group_index] = lognormal_rng(
            group_dist_mean[dist_group_index, 1],
            group_dist_sd[dist_group_index, 1]
          );
        }
      }
    }

    for (dist_group_index in 1:num_discrete_dist) {
      vector[num_clusters] counterfactual_dist =
        all_cluster_dist[, dist_group_index];
      for (treatment_index in 1:num_treatments) {
        int cell_index = 0;
        vector[num_clusters] mu_rep;
        vector[num_clusters] benefit_cost;
        vector[num_clusters] cutoff;
        vector[num_clusters] probability;
        real group_numerator = 0;

        for (candidate_index in 1:num_dist_group_treatments) {
          if (cluster_treatment_map[candidate_index, 1] == treatment_index &&
              cluster_treatment_map[candidate_index, 2] == dist_group_index) {
            cell_index = candidate_index;
          }
        }
        if (cell_index == 0) {
          reject("No treatment-by-distance cell found in compact core GQ.");
        }

        mu_rep = base_mu_rep * inv_logit(
          rep_intercept[treatment_index] +
          rep_dist_slope[treatment_index] * counterfactual_dist
        );
        benefit_cost =
          rep_vector(structural_treatment_effect[cell_index], num_clusters) -
          dist_beta_v[1] * counterfactual_dist + cluster_shock;
        for (cluster_index in 1:num_clusters) {
          cutoff[cluster_index] = core_gq_find_fixedpoint_safe(
            benefit_cost[cluster_index],
            mu_rep[cluster_index],
            total_error_sd,
            u_sd,
            use_u_in_delta,
            alg_sol_rel_tol,
            alg_sol_f_tol,
            alg_sol_max_steps
          );
        }
        probability = Phi_approx(-cutoff / total_error_sd);

        for (cluster_index in 1:num_clusters) {
          real size_weight = sum(cluster_size[cluster_index]);
          group_numerator += size_weight * probability[cluster_index];
          if (cluster_dist_treatment_id[cluster_index] == dist_group_index) {
            combined_numerator[treatment_index] +=
              size_weight * probability[cluster_index];
          }
        }
        core_compact_takeup_level[dist_group_index + 1, treatment_index] =
          group_numerator / total_weight;
      }
    }
    core_compact_takeup_level[1] = combined_numerator' / total_weight;

    if (num_roc_distances < compact_sm_roc_index[3]) {
      reject("Compact core GQ requires ROC distances through 2500m.");
    }
    for (compact_dist_index in 1:3) {
      real roc_distance =
        roc_distances[compact_sm_roc_index[compact_dist_index]];
      vector[num_clusters] roc_cluster_dist =
        rep_vector(roc_distance, num_clusters);
      vector[num_clusters] control_mu;
      vector[num_clusters] control_benefit;
      int control_cell = 0;

      for (candidate_index in 1:num_dist_group_treatments) {
        if (cluster_treatment_map[candidate_index, 1] ==
            roc_compare_treatment_id_right) {
          control_cell = candidate_index;
        }
      }
      control_mu = base_mu_rep * inv_logit(
        rep_intercept[roc_compare_treatment_id_right] +
        rep_dist_slope[roc_compare_treatment_id_right] * roc_cluster_dist
      );
      control_benefit = rep_vector(
        structural_treatment_effect[control_cell], num_clusters
      ) - dist_beta_v[1] * roc_cluster_dist + cluster_shock;

      for (treatment_index in 1:num_treatments) {
        int treatment_cell = 0;
        vector[num_clusters] curr_mu;
        vector[num_clusters] curr_mu_deriv;
        vector[num_clusters] curr_mu_share;
        vector[num_clusters] curr_benefit;
        matrix[num_clusters, 9] roc_results;

        for (candidate_index in 1:num_dist_group_treatments) {
          if (cluster_treatment_map[candidate_index, 1] == treatment_index) {
            treatment_cell = candidate_index;
          }
        }
        curr_mu_share = inv_logit(
          rep_intercept[treatment_index] +
          rep_dist_slope[treatment_index] * roc_cluster_dist
        );
        curr_mu = base_mu_rep * curr_mu_share;
        curr_mu_deriv = base_mu_rep * rep_dist_slope[treatment_index] .*
          curr_mu_share .* (1 - curr_mu_share);
        curr_benefit = rep_vector(
          structural_treatment_effect[treatment_cell], num_clusters
        ) - dist_beta_v[1] * roc_cluster_dist + cluster_shock;

        for (cluster_index in 1:num_clusters) {
          roc_results[cluster_index] = core_gq_calculate_roc(
            curr_benefit[cluster_index],
            control_benefit[cluster_index],
            total_error_sd,
            u_sd,
            dist_beta_v[1],
            curr_mu[cluster_index],
            control_mu[cluster_index],
            curr_mu_deriv[cluster_index],
            use_u_in_delta,
            alg_sol_rel_tol,
            alg_sol_f_tol,
            alg_sol_max_steps
          )';
        }
        core_compact_sm[compact_dist_index, treatment_index] =
          mean(roc_results[, 7]);
        core_compact_sm_rescaled[compact_dist_index, treatment_index] =
          core_compact_sm[compact_dist_index, treatment_index] /
          dist_beta_v[1];
      }
    }
  }
}
