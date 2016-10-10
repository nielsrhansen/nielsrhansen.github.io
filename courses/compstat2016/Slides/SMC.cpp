// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
List bootARfast(NumericVector y, int n, double alpha, double sigma) {
  int t = y.size();
  NumericMatrix x(n, t), u(n, t);
  NumericVector uu(n), d;
  IntegerVector resamp(n), indices = seq_len(n) - 1;
  NumericVector e1 = rnorm(n * t, 0, 1);
  NumericMatrix e(n, t, e1.begin());  
  x(_, 0) = e(_, 0) / sqrt(1 - alpha * alpha);
  for(int s = 0; s < t - 1; s++) {
    d = x(_, s);
    d = (d - y(s)) / sigma;
    uu = exp( - (d * d) / 2);
    u(_, s) = uu;
    resamp = RcppArmadillo::sample(indices, n, TRUE, uu);
    for(int k = 0; k < n; k++) {
      x(k, s + 1) = alpha * x(resamp(k), s) + e(k, s + 1);
    }
  }
  d = x(_, t - 1);
  d = (d - y(t - 1)) / sigma;
  u(_, t - 1) = exp( - (d * d) / 2);
  List ret; ret["x"] = x; ret["u"] = u;
  return ret;
}


