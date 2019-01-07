#include <RcppEigen.h>
#include <cmath>
#include "draw_samples.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


// [[Rcpp::export]]
Eigen::MatrixXd lm_proper_cpp(Eigen::MatrixXd X,
                            Eigen::VectorXd y,
                            Eigen::VectorXd beta_prior_mean,
                            Eigen::MatrixXd beta_prior_cov,
                            double tau_prior_shape,
                            double tau_prior_rate,
                            int iterations,
                            int burnin,
                            int thin,
                            Eigen::VectorXd start_beta){
  int n = X.rows();
  int r = X.cols();

  Eigen::MatrixXd xtx = X.transpose() * X;
  Eigen::MatrixXd xty = X.transpose() * y;
  Eigen::VectorXd beta_hat = xtx.inverse() * xty;
  double sse = (y - X * beta_hat).squaredNorm();
  Eigen::MatrixXd cov_inv = beta_prior_cov.inverse();

  double tau_post_shape = tau_prior_shape + (n / 2.0);
  Eigen::VectorXd diff_beta(r);

  Eigen::MatrixXd beta_out(iterations, r);
  Eigen::VectorXd sigma_out(iterations);

  Eigen::VectorXd temp_beta = start_beta;
  double temp_tau = 1.0;
  double temp_rate = 0.0;

  Eigen::MatrixXd V_inv = temp_tau * xtx * cov_inv;
  Eigen::MatrixXd V(V_inv.rows(), V_inv.cols());

  for(int i=0; i < burnin; i++){

    if (i % 1000 == 0) checkUserInterrupt();

    diff_beta = temp_beta - beta_hat;
    temp_rate = tau_prior_rate + 0.5 * (sse + diff_beta.transpose() * xtx * diff_beta);
    temp_tau = uv_gamma(tau_post_shape, temp_rate);
    V_inv = temp_tau * xtx + cov_inv;
    V = V_inv.inverse();
    temp_beta = mv_normal(V * (temp_tau * xty + cov_inv * start_beta),
                               V);
  }

  for(int i=0; i < iterations; i++){

    if (i % 1000 == 0) checkUserInterrupt();

    for(int j=0; j < thin; j++){
      diff_beta = temp_beta - beta_hat;
      temp_rate = tau_prior_rate + 0.5 * (sse + diff_beta.transpose() * xtx * diff_beta);
      temp_tau = uv_gamma(tau_post_shape, temp_rate);
      V_inv = temp_tau * xtx + cov_inv;
      V = V_inv.inverse();
      temp_beta = mv_normal(V * (temp_tau * xty + cov_inv * start_beta),
                                 V);
    }
    beta_out.row(i) = temp_beta.transpose();
    sigma_out(i) = 1.0 / std::pow(temp_tau, 0.5);
  }

  Eigen::MatrixXd result(beta_out.rows(), beta_out.cols() + 1);
  result << beta_out, sigma_out;
  return result;

  //return List::create(
  //  Named("beta") = beta_out,
 //   Named("sigma") = sigma_out
 // );

}
