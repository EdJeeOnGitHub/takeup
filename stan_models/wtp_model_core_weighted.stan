// Fit-105 WTP likelihood with cluster likelihood weights. Priors are not
// weighted. With unit weights this is identical to wtp_model_section.stan.

hyper_wtp_mu ~ normal(0, 2);
strata_wtp_mu_tau ~ normal(0, 1);
wtp_sigma ~ normal(0, 1);

if (fit_wtp_model_to_data) {
  int wtp_stratum_pos = 1;

  for (stratum_index in 1:num_strata) {
    int curr_wtp_stratum_size = wtp_strata_sizes[stratum_index];
    int wtp_stratum_end = wtp_stratum_pos + curr_wtp_stratum_size - 1;

    for (wtp_obs_index in wtp_stratum_pos:wtp_stratum_end) {
      real wtp_lp;
      if (core_has_nonunit_weight) {
        // Algebraically identical but stable in the extreme tails reached by
        // some weighted optimizations.
        if (wtp_response[wtp_obs_index] == -1) {
          if (gift_choice[wtp_obs_index] == -1) {
            wtp_lp = normal_lcdf(
              -scaled_wtp_offer[wtp_obs_index] |
              strata_wtp_mu[stratum_index], wtp_sigma
            );
          } else {
            wtp_lp = normal_lccdf(
              scaled_wtp_offer[wtp_obs_index] |
              strata_wtp_mu[stratum_index], wtp_sigma
            );
          }
        } else {
          real endpoint =
            gift_choice[wtp_obs_index] * scaled_wtp_offer[wtp_obs_index];
          real lower_endpoint = fmin(0, endpoint);
          real upper_endpoint = fmax(0, endpoint);
          wtp_lp = log_diff_exp(
            normal_lcdf(
              upper_endpoint | strata_wtp_mu[stratum_index], wtp_sigma
            ),
            normal_lcdf(
              lower_endpoint | strata_wtp_mu[stratum_index], wtp_sigma
            )
          );
        }
      } else if (wtp_response[wtp_obs_index] == -1) {
        // Preserve the historical arithmetic exactly for unit-weight fits.
        if (gift_choice[wtp_obs_index] == -1) {
          wtp_lp = log(Phi_approx(
            (-scaled_wtp_offer[wtp_obs_index] -
             strata_wtp_mu[stratum_index]) / wtp_sigma
          ));
        } else {
          wtp_lp = log(1 - Phi_approx(
            (scaled_wtp_offer[wtp_obs_index] -
             strata_wtp_mu[stratum_index]) / wtp_sigma
          ));
        }
      } else {
        wtp_lp = log(
          gift_choice[wtp_obs_index] *
          (Phi_approx(
             (gift_choice[wtp_obs_index] *
              scaled_wtp_offer[wtp_obs_index] -
              strata_wtp_mu[stratum_index]) / wtp_sigma
           ) - Phi_approx(-strata_wtp_mu[stratum_index] / wtp_sigma))
        );
      }
      target += core_cluster_weight[wtp_cluster_id[wtp_obs_index]] * wtp_lp;
    }

    wtp_stratum_pos = wtp_stratum_end + 1;
  }
}
