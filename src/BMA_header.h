// include guard
#ifndef __BMA_H_INCLUDED__
#define __BMA_H_INCLUDED__

//=================================
// forward declared dependencies
class BMA;

//=================================
// included dependencies
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <Rcpp/Benchmark/Timer.h>

//=================================
// the actual class
class BMA 
{
public:
  void set_seed(double seed);
  double square(double y);
  double adjust_acceptance(double accept,double sgm,double target );
  arma::uvec arma_setdiff(arma::uvec x, arma::uvec y);
};

#endif // __BMA_H_INCLUDED__ 



