// Memory-efficient social-multiplier generated quantities for the individual
// distance robustness models. Keep only treatment means at the three distances
// used in the paper table; never retain distance-by-individual arrays.
//
// roc_distances is 0m, 100m, ..., 5000m, standardized in the R data builder.
array[3] int compact_sm_roc_index = { 6, 16, 26 };

matrix[3, num_treatments] compact_sm;
matrix[3, num_treatments] compact_dist_beta_v;
matrix[3, num_treatments] compact_sm_delta_part;
matrix[3, num_treatments] compact_sm_mu_part;
matrix[3, num_treatments] compact_sm_rescaled;
matrix[3, num_treatments] compact_sm_delta_part_rescaled;
matrix[3, num_treatments] compact_sm_mu_part_rescaled;
// Rows are Combined, Close, and Far. These are sufficient to reconstruct the
// paper's arm ATEs and Far-minus-Close effects without saving cluster arrays.
matrix[3, num_treatments] compact_takeup_level;

{
  matrix[num_clusters, num_dist_group_treatments] structural_cluster_benefit =
    rep_matrix(structural_treatment_effect', num_clusters) +
    rep_matrix(
      structural_beta_cluster + structural_beta_county[cluster_county_id],
      num_dist_group_treatments
    );
  matrix[num_clusters, num_discrete_dist] compact_all_cluster_dist;
  vector[num_treatments] compact_combined_takeup_numerator =
    rep_vector(0, num_treatments);
  vector[num_discrete_dist] compact_dist_group_weight =
    rep_vector(0, num_discrete_dist);
  real compact_total_weight = 0;

  // Reproduce the historical counterfactual-distance construction locally.
  // Observed distances are retained for the realized distance group; the
  // other potential distance is drawn from its fitted lognormal distribution.
  for (cluster_index in 1:num_clusters) {
    real cluster_weight = sum(cluster_size[cluster_index]);
    int observed_dist_group =
      cluster_dist_treatment_id[cluster_index];

    compact_total_weight += cluster_weight;
    compact_dist_group_weight[observed_dist_group] += cluster_weight;
    for (dist_group_index in 1:num_discrete_dist) {
      if (observed_dist_group == dist_group_index) {
        compact_all_cluster_dist[cluster_index, dist_group_index] =
          cluster_standard_dist[cluster_index];
      } else {
        compact_all_cluster_dist[cluster_index, dist_group_index] =
          lognormal_rng(
            group_dist_mean[dist_group_index, 1],
            group_dist_sd[dist_group_index, 1]
          );
      }
    }
  }

  for (dist_group_index in 1:num_discrete_dist) {
    vector[num_clusters] counterfactual_dist =
      compact_all_cluster_dist[, dist_group_index];

    for (treatment_index in 1:num_treatments) {
      int dist_treatment_index = 0;
      vector[num_clusters] curr_mu_rep;
      vector[num_clusters] curr_benefit_cost;
      vector[num_clusters] curr_cutoff;
      vector[num_clusters] curr_probability;
      real dist_group_numerator = 0;

      for (candidate_index in 1:num_dist_group_treatments) {
        if (
            cluster_treatment_map[candidate_index, 1] == treatment_index &&
            cluster_treatment_map[candidate_index, 2] == dist_group_index) {
          dist_treatment_index = candidate_index;
        }
      }
      if (dist_treatment_index == 0) {
        reject("No treatment-by-distance cell found in compact ATE GQ.");
      }

      curr_mu_rep = calculate_mu_rep(
        { treatment_index },
        counterfactual_dist,
        base_mu_rep,
        1,
        beliefs_treatment_map_design_matrix,
        centered_cluster_beta_beliefs,
        centered_cluster_dist_beta_beliefs,
        mu_rep_type
      );
      curr_benefit_cost =
        structural_cluster_benefit[, dist_treatment_index] -
        param_dist_cost(
          counterfactual_dist,
          cluster_linear_dist_cost[, dist_treatment_index],
          cluster_quadratic_dist_cost[, dist_treatment_index]
        );
      curr_cutoff = map_find_fixedpoint_solution(
        curr_benefit_cost,
        curr_mu_rep,
        rep_vector(total_error_sd[treatment_index], num_clusters),
        rep_vector(u_sd[treatment_index], num_clusters),
        use_u_in_delta,
        alg_sol_rel_tol,
        alg_sol_f_tol,
        alg_sol_max_steps
      );
      curr_probability = Phi_approx(
        -curr_cutoff / total_error_sd[treatment_index]
      );

      for (cluster_index in 1:num_clusters) {
        real cluster_weight = sum(cluster_size[cluster_index]);
        dist_group_numerator +=
          cluster_weight * curr_probability[cluster_index];
        if (cluster_dist_treatment_id[cluster_index] == dist_group_index) {
          compact_combined_takeup_numerator[treatment_index] +=
            cluster_weight * curr_probability[cluster_index];
        }
      }
      compact_takeup_level[dist_group_index + 1, treatment_index] =
        dist_group_numerator / compact_total_weight;
    }
  }
  compact_takeup_level[1] =
    compact_combined_takeup_numerator' / compact_total_weight;

  if (num_roc_distances < compact_sm_roc_index[3]) {
    reject("Compact social-multiplier GQ requires roc distances through 2500m.");
  }

  for (compact_dist_index in 1:3) {
    int roc_dist_index = compact_sm_roc_index[compact_dist_index];
    real roc_distance = roc_distances[roc_dist_index];
    vector[num_clusters] roc_cluster_dist =
      rep_vector(roc_distance, num_clusters);
    matrix[num_clusters, num_treatments] curr_net_benefit =
      structural_cluster_benefit[, :num_treatments] -
      param_dist_cost(
        roc_distance,
        cluster_linear_dist_cost[, :num_treatments],
        cluster_quadratic_dist_cost[, :num_treatments]
      );
    matrix[num_clusters, 2] curr_cluster_mu_rep_control =
      calculate_mu_rep_deriv(
        roc_compare_treatment_id_right,
        roc_cluster_dist,
        base_mu_rep,
        1,
        beliefs_treatment_map_design_matrix,
        centered_cluster_beta_1ord,
        centered_cluster_dist_beta_1ord,
        mu_rep_type
      );

    for (treatment_index in 1:num_treatments) {
      matrix[num_clusters, 2] curr_cluster_mu_rep;
      matrix[num_clusters, 9] roc_results;
      vector[num_clusters] denominator;
      vector[num_clusters] delta_part;
      vector[num_clusters] mu_part;

      if (treatment_index > 1) {
        curr_cluster_mu_rep = calculate_mu_rep_deriv(
          treatment_index,
          roc_cluster_dist,
          base_mu_rep,
          1,
          beliefs_treatment_map_design_matrix,
          centered_cluster_beta_1ord,
          centered_cluster_dist_beta_1ord,
          mu_rep_type
        );
      } else {
        curr_cluster_mu_rep = curr_cluster_mu_rep_control;
      }

      roc_results = map_calculate_roc(
        curr_net_benefit[, treatment_index],
        curr_net_benefit[, roc_compare_treatment_id_right],
        rep_vector(total_error_sd[treatment_index], num_clusters),
        rep_vector(
          total_error_sd[roc_compare_treatment_id_right],
          num_clusters
        ),
        rep_vector(u_sd[treatment_index], num_clusters),
        rep_vector(u_sd[roc_compare_treatment_id_right], num_clusters),
        cluster_linear_dist_cost[, treatment_index],
        curr_cluster_mu_rep[, 1],
        curr_cluster_mu_rep_control[, 1],
        curr_cluster_mu_rep[, 2],
        use_u_in_delta,
        alg_sol_rel_tol,
        alg_sol_f_tol,
        alg_sol_max_steps
      );

      denominator =
        1 + curr_cluster_mu_rep[, 1] .* roc_results[, 9];
      delta_part = -dist_beta_v[1] ./ denominator;
      mu_part =
        curr_cluster_mu_rep[, 2] .* roc_results[, 8] ./ denominator;

      compact_sm[compact_dist_index, treatment_index] =
        mean(roc_results[, 7]);
      compact_dist_beta_v[compact_dist_index, treatment_index] =
        dist_beta_v[1];
      compact_sm_delta_part[compact_dist_index, treatment_index] =
        mean(delta_part);
      compact_sm_mu_part[compact_dist_index, treatment_index] =
        mean(mu_part);
      compact_sm_rescaled[compact_dist_index, treatment_index] =
        compact_sm[compact_dist_index, treatment_index] / dist_beta_v[1];
      compact_sm_delta_part_rescaled[compact_dist_index, treatment_index] =
        compact_sm_delta_part[compact_dist_index, treatment_index] /
        dist_beta_v[1];
      compact_sm_mu_part_rescaled[compact_dist_index, treatment_index] =
        compact_sm_mu_part[compact_dist_index, treatment_index] /
        dist_beta_v[1];
    }
  }
}
