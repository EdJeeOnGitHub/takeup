// Log-likelihood totals used for Bayesian data-to-parameter sensitivity.
// These reproduce the four data likelihood blocks in struct_section_model.stan.
// Priors are deliberately excluded.

real sensitivity_log_lik_takeup = 0;
real sensitivity_log_lik_beliefs = 0;
real sensitivity_log_lik_wtp = 0;
real sensitivity_log_lik_distance = 0;

if (fit_model_to_data) {
  if (use_binomial) {
    for (included_cluster_index in 1:num_included_clusters) {
      int cluster_index = included_clusters[included_cluster_index];
      sensitivity_log_lik_takeup += binomial_lpmf(
        cluster_takeup_count[cluster_index, 1] |
        cluster_size[cluster_index, 1],
        structural_cluster_takeup_prob[cluster_index]
      );
    }
  } else {
    for (included_obs_index in 1:num_included_monitored_obs) {
      int obs_index = included_monitored_obs[included_obs_index];
      sensitivity_log_lik_takeup += bernoulli_lpmf(
        takeup[obs_index] |
        structural_cluster_takeup_prob[obs_cluster_id[obs_index]]
      );
    }
  }
}

if (fit_beliefs_model_to_data) {
  vector[num_beliefs_obs] sensitivity_beliefs_predictor_1ord =
    calculate_beliefs_latent_predictor(
      beliefs_treatment_design_matrix,
      centered_obs_beta_1ord,
      centered_obs_dist_beta_1ord,
      cluster_standard_dist[beliefs_cluster_index]
    );
  vector[num_beliefs_obs] sensitivity_beliefs_predictor_2ord =
    calculate_beliefs_latent_predictor(
      beliefs_treatment_design_matrix,
      centered_obs_beta_2ord,
      centered_obs_dist_beta_2ord,
      cluster_standard_dist[beliefs_cluster_index]
    );

  for (belief_index in 1:num_beliefs_obs) {
    if (belief_observed[belief_index]) {
      sensitivity_log_lik_beliefs += binomial_logit_lpmf(
        num_knows_1ord[belief_index] |
        num_recognized[belief_index],
        sensitivity_beliefs_predictor_1ord[belief_index]
      );
      sensitivity_log_lik_beliefs += binomial_logit_lpmf(
        num_knows_2ord[belief_index] |
        num_recognized[belief_index],
        sensitivity_beliefs_predictor_2ord[belief_index]
      );
    }
  }
}

if (fit_wtp_model_to_data) {
  int sensitivity_wtp_stratum_pos = 1;

  for (stratum_index in 1:num_strata) {
    int curr_wtp_stratum_size = wtp_strata_sizes[stratum_index];
    int sensitivity_wtp_stratum_end =
      sensitivity_wtp_stratum_pos + curr_wtp_stratum_size - 1;

    for (wtp_obs_index in sensitivity_wtp_stratum_pos:sensitivity_wtp_stratum_end) {
      if (wtp_response[wtp_obs_index] == -1) {
        if (gift_choice[wtp_obs_index] == -1) {
          sensitivity_log_lik_wtp += log(Phi_approx(
            (-scaled_wtp_offer[wtp_obs_index] - strata_wtp_mu[stratum_index]) /
            wtp_sigma
          ));
        } else {
          sensitivity_log_lik_wtp += log1m(Phi_approx(
            (scaled_wtp_offer[wtp_obs_index] - strata_wtp_mu[stratum_index]) /
            wtp_sigma
          ));
        }
      } else {
        sensitivity_log_lik_wtp += log(
          gift_choice[wtp_obs_index] *
          (
            Phi_approx(
              (
                gift_choice[wtp_obs_index] * scaled_wtp_offer[wtp_obs_index] -
                strata_wtp_mu[stratum_index]
              ) / wtp_sigma
            ) -
            Phi_approx(-strata_wtp_mu[stratum_index] / wtp_sigma)
          )
        );
      }
    }

    sensitivity_wtp_stratum_pos = sensitivity_wtp_stratum_end + 1;
  }
}

if (fit_dist_model_to_data) {
  for (cluster_index in 1:num_clusters) {
    int curr_assigned_dist_group =
      cluster_treatment_map[
        cluster_assigned_dist_group_treatment[cluster_index], 2
      ];
    vector[num_dist_group_mix] sensitivity_group_dist_mix_lp;

    for (group_dist_mix_index in 1:num_dist_group_mix) {
      real sensitivity_dist_location =
        group_dist_mean[curr_assigned_dist_group, group_dist_mix_index] +
        county_dist_effect[
          curr_assigned_dist_group, cluster_county_id[cluster_index]
        ] +
        cluster_dist_effect[curr_assigned_dist_group, cluster_index];

      if (lognormal_dist_model) {
        sensitivity_group_dist_mix_lp[group_dist_mix_index] =
          lognormal_lpdf(
            cluster_standard_dist[cluster_index] |
            sensitivity_dist_location,
            group_dist_sd[curr_assigned_dist_group, group_dist_mix_index]
          ) +
          log(group_dist_mix[
            curr_assigned_dist_group, group_dist_mix_index
          ]);
      } else {
        sensitivity_group_dist_mix_lp[group_dist_mix_index] =
          normal_lpdf(
            cluster_standard_dist[cluster_index] |
            sensitivity_dist_location,
            group_dist_sd[curr_assigned_dist_group, group_dist_mix_index]
          ) +
          normal_lccdf(
            0 |
            sensitivity_dist_location,
            group_dist_sd[curr_assigned_dist_group, group_dist_mix_index]
          ) +
          log(group_dist_mix[
            curr_assigned_dist_group, group_dist_mix_index
          ]);
      }
    }

    sensitivity_log_lik_distance +=
      log_sum_exp(sensitivity_group_dist_mix_lp);
  }
}
