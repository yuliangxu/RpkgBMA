
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include "BMA_header.h"


void BMA::set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

double BMA::square(double y){
  return y*y;
}

double BMA::adjust_acceptance(double accept,double sgm,double target = 0.1){
  double y = 1. + 1000.*(accept-target)*(accept-target)*(accept-target);
  if (y < .9)
    y = .9;
  if (y > 1.1)
    y = 1.1;
  sgm *= y;
  return sgm;
}

arma::uvec BMA::arma_setdiff(arma::uvec x, arma::uvec y){
  
  x = arma::unique(x);
  y = arma::unique(y);
  
  for (arma::uword j = 0; j < y.n_elem; j++) {
    arma::uvec q1 = arma::find(x == y[j]);
    if (!q1.empty()) {
      x.shed_row(q1(0));
    }
  }
  
  return x;
}
