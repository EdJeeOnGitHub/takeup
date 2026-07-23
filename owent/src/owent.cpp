// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/owens_t.hpp>

// [[Rcpp::export]]
Rcpp::NumericVector owens_t_vec(Rcpp::NumericVector h, Rcpp::NumericVector a) {
    int n = h.size();
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        out[i] = boost::math::owens_t(h[i], a[i]);
    }
    return out;
}
