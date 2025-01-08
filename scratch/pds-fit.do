log using lpm-regressions, replace
clear


* Functions ----------------------------------------------------------------------
capture program drop hyp_tests
program define hyp_tests
	display "Coef Value: " _b[`1']
	display "Two-Sided Test Against 0" 
	test _b[`1'] = 0
	local sign_val = sign(_b[`1'])
	display "H0: coef <= 0 p-value = " ttail(r(df_r), `sign_val'*sqrt(r(F)))
	display "H0: coef >= 0 p-value = " 1 - ttail(r(df_r), `sign_val'*sqrt(r(F)))
end


capture program drop hyp_tests_chi
program define hyp_tests_chi
	display "Coef Value: " _b[`1']
	display "Two-Sided Test Against 0" 
	test _b[`1'] = 0
	local sign_val = sign(_b[`1'])
	display "H0: coef <= 0 p-value = " 1 - normal(`sign_val'*sqrt(r(chi2)))
	display "H0: coef >= 0 p-value = "  normal(`sign_val'*sqrt(r(chi2)))
end



* Data ----------------------------------------------------------------------
import delimited "temp-data/analysis-cluster-covariate-data.csv"


list dewormed in 1/5

generate dewormed_num = (dewormed == "TRUE")
generate female_num = (female == "FALSE")
generate have_phone_num = (have_phone_lgl == "TRUE")
generate dist_pot_far = (distpotgroup == "far")


encode assignedtreatment, generate(assigned_treatment_fac)
recode assigned_treatment_fac (3 = 1 control) (4 = 2 ink) (2 = 3 calendar) (1 = 4 bracelet), gen(assigned_treatment)
encode county, generate(county_fac)





global cov_vars n_per_cluster      completed_primary         ///      
	 floor_tile_cement         ethnicity_luhya           ///
	 religion_christianity     treated_lgl               ///
	 adults_can_get_worms      correct_when_treat        ///
	 externality_omnibus       know_medicine_stops_worms ///
	 know_children_get_worms   sick_worms_only           ///
	 treated_past_year         praise_immuniz            ///
	 praise_dewor              stigma_immuniz            ///
	 stigma_dewor              female_num                ///
	 have_phone_num            agecensus                 



rename assigned_treatment a_treat
rename dist_pot_far dpf

* treat x dist group, county fe & clustered
reg dewormed_num i.a_treat##i.dpf  i.county_fac, cluster(clusteridx)

* P1, Negative Effect of Distance
hyp_tests 1.dpf
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#1.dpf



* treat x distance cts, county fe & clustered
reg dewormed_num i.a_treat##c.standard_clusterdisttopot i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests standard_clusterdisttopot
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#c.standard_clusterdisttopot

* treat x dist group, controls + county fe & clustered
reg dewormed_num i.a_treat##i.dpf $cov_vars i.county_fac, cluster(clusteridx)

* P1, Negative Effect of Distance
hyp_tests 1.dpf
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#1.dpf

* treat x distance cts, controls + county fe & clustered
reg dewormed_num i.a_treat##c.standard_clusterdisttopot $cov_vars i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests standard_clusterdisttopot
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#c.standard_clusterdisttopot



* Double LASSO, treat x dist group, clustered
pdslasso dewormed_num  i.a_treat##dpf   ($cov_vars i.county_fac), cluster(clusteridx) partial(i.county_fac)

* P1, Negative Effect of Distance
hyp_tests_chi 1.dpf
* P2, Positive Effect Bracelets
hyp_tests_chi 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests_chi 4.a_treat#1.dpf



* Double LASSO, treat x continuous dist group, clustered
pdslasso dewormed_num  i.a_treat##c.standard_clusterdisttopot   ($cov_vars  i.county_fac), cluster(clusteridx) partial(i.county_fac)


* P1, Negative Effect of Distance
hyp_tests_chi standard_clusterdisttopot
* P2, Positive Effect Bracelets
hyp_tests_chi 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests_chi 4.a_treat#c.standard_clusterdisttopot

* Double LASSO, treat x dist cts, clustered
pdslasso dewormed_num  i.a_treat##c.standard_clusterdisttopot ($cov_vars i.county_fac), cluster(clusteridx) pnotpen(i.county_fac)



* Double LASSO on Distance
pdslasso dewormed_num  dpf ($cov_vars i.a_treat i.county_fac), cluster(clusteridx) pnotpen(i.county_fac i.a_treat)


global lasso_cov_vars floor_tile_cement ///
	female_num ///
	have_phone_num ///
	agecensus 

* Manual PDS LASSO, Dist cts
reg dewormed_num i.a_treat##c.standard_clusterdisttopot $lasso_cov_vars i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests standard_clusterdisttopot
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#c.standard_clusterdisttopot


* Manual PDS LASSO, Dist Group
reg dewormed_num i.a_treat##dpf $lasso_cov_vars i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests 1.dpf
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#1.dpf


log close
translate lpm-regressions.smcl lpm-regressions.pdf

log using main-spec
* Main Spec - Discrete Distance
* treat x dist group, county fe & clustered
reg dewormed_num i.a_treat##i.dpf  i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests 1.dpf
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#1.dpf

* treat x dist group, controls + county fe & clustered
reg dewormed_num i.a_treat##i.dpf $cov_vars i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests 1.dpf
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#1.dpf


reg dewormed_num i.a_treat##dpf $lasso_cov_vars i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests 1.dpf
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#1.dpf




* Main Spec - Dist CTS
reg dewormed_num i.a_treat##c.standard_clusterdisttopot  i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests standard_clusterdisttopot
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#c.standard_clusterdisttopot


reg dewormed_num i.a_treat##c.standard_clusterdisttopot $cov_vars  i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests standard_clusterdisttopot
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#c.standard_clusterdisttopot

reg dewormed_num i.a_treat##c.standard_clusterdisttopot $lasso_cov_vars  i.county_fac, cluster(clusteridx)
* P1, Negative Effect of Distance
hyp_tests standard_clusterdisttopot
* P2, Positive Effect Bracelets
hyp_tests 4.a_treat
* P3, Positive Coef on Far x Bracelet
hyp_tests 4.a_treat#c.standard_clusterdisttopot


log close
translate main-spec.smcl main-spec.pdf
