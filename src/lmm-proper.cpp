#include <RcppEigen.h>
#include "draw_samples.h"
#include "conv-checks.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
List lmm_proper_cpp(Eigen::MatrixXd x,
                    Eigen::MatrixXd z,
                    Eigen::VectorXd y,
                    Eigen::VectorXd beta_prior_mean,
                    Eigen::MatrixXd beta_prior_cov,
                    Eigen::VectorXd lambda_prior_shape,
                    Eigen::VectorXd lambda_prior_rate,
                    int iterations,
                    int burnin,
                    int thin,
                    Eigen::VectorXd start_theta){

  int n = x.rows();
  int p = x.cols();
  int q = z.cols();

  Eigen::MatrixXd xt = x.transpose();
  Eigen::MatrixXd zt = z.transpose();
  Eigen::MatrixXd xtx = xt * x;
  Eigen::MatrixXd xtz = xt * z;
  Eigen::MatrixXd xty = xt * y;
  Eigen::MatrixXd ztx = zt * x;
  Eigen::MatrixXd ztz = zt * z;
  Eigen::MatrixXd zty = zt * y;

  Eigen::VectorXd lambda_post_shape(2);
  lambda_post_shape(0) = lambda_prior_shape(0) + n * 0.5;
  lambda_post_shape(1) = lambda_prior_shape(1) + q * 0.5;

  Eigen::MatrixXd beta_cov_inv = beta_prior_cov.inverse();

  Eigen::MatrixXd w(n, p + q);
  w.leftCols(p) = x;
  w.rightCols(q) = z;

  Eigen::MatrixXd theta_out(iterations, p + q);
  Eigen::MatrixXd sigma_out(iterations, 2);
  Eigen::VectorXd temp_lambda(2);
  temp_lambda << 0.0, 0.0;
  Eigen::VectorXd temp_theta = start_theta;

  Eigen::MatrixXd iq = Eigen::MatrixXd::Identity(q, q);

  Eigen::MatrixXd V_inv(p + q, p + q);
  Eigen::MatrixXd V(p + q, p + q);
  Eigen::VectorXd eta(p + q);
  Eigen::VectorXd diff(n);

  conv_check_lmm_proper(z, lambda_prior_shape);

  for(int i = 0; i < burnin; i++){

    if (i % 1000 == 0) checkUserInterrupt();

    diff = y - w * temp_theta;
    temp_lambda(0) = uv_gamma(lambda_post_shape(0),
                lambda_prior_rate(0) + 0.5 * diff.squaredNorm());
    temp_lambda(1) = uv_gamma(lambda_post_shape(1),
                lambda_prior_rate(1) + 0.5 * temp_theta.tail(q).squaredNorm());

    V_inv.topLeftCorner(p, p) = temp_lambda(0) * xtx + beta_cov_inv;
    V_inv.topRightCorner(p, q) = temp_lambda(0) * xtz;
    V_inv.bottomLeftCorner(q, p) = temp_lambda(0) * ztx;
    V_inv.bottomRightCorner(q, q) = temp_lambda(0) * ztz + temp_lambda(1) * iq;
    V = V_inv.inverse();
    eta.head(p) = temp_lambda(0) * xty + beta_cov_inv * beta_prior_mean;
    eta.tail(q) = temp_lambda(0) * zty;
    temp_theta = mv_normal(V * eta, V);
  }

  for(int i = 0; i < iterations; i++){

    if (i % 1000 == 0) checkUserInterrupt();

    for(int j = 0; j < thin; j++){
      diff = y - w * temp_theta;
      temp_lambda(0) = uv_gamma(lambda_post_shape(0),
                  lambda_prior_rate(0) + 0.5 * diff.squaredNorm());
      temp_lambda(1) = uv_gamma(lambda_post_shape(1),
                  lambda_prior_rate(1) + 0.5 * temp_theta.tail(q).squaredNorm());

      V_inv.topLeftCorner(p, p) = temp_lambda(0) * xtx + beta_cov_inv;
      V_inv.topRightCorner(p, q) = temp_lambda(0) * xtz;
      V_inv.bottomLeftCorner(q, p) = temp_lambda(0) * ztx;
      V_inv.bottomRightCorner(q, q) = temp_lambda(0) * ztz + temp_lambda(1) * iq;
      V = V_inv.inverse();
      eta.head(p) = temp_lambda(0) * xty + beta_cov_inv * beta_prior_mean;
      eta.tail(q) = temp_lambda(0) * zty;
      temp_theta = mv_normal(V * eta, V);
    }

    theta_out.row(i) = temp_theta;
    sigma_out.row(i) = 1.0 / temp_lambda.array().sqrt();

  }


  Eigen::MatrixXd beta_out = theta_out.leftCols(p);
  Eigen::MatrixXd u_out = theta_out.rightCols(q);

  return List::create(
    Named("beta") = beta_out,
    Named("u") = u_out,
    Named("sigma") = sigma_out
  );

}
