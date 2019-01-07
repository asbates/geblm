#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


// draw multivariate normal sample
// [[Rcpp::export]]
Eigen::VectorXd mv_normal(Eigen::VectorXd mu,
                               Eigen::MatrixXd sigma){
  int r = mu.size();
  Eigen::VectorXd standard_normal(r);
  Eigen::VectorXd result(r);

  standard_normal = as<Eigen::VectorXd>(Rcpp::rnorm(r, 0, 1));
  result = mu + sigma.llt().matrixL() * standard_normal;
  return result;
}


// draw univariate gamma sample
// [[Rcpp::export]]
double uv_gamma(double shape, double rate){
  double result = R::rgamma(shape, 1.0 / rate);
  return result;
}

