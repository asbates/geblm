
#ifndef DRAW_SAMPLES_H
#define DRAW_SAMPLES_H

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;

Eigen::VectorXd mv_normal(Eigen::VectorXd mu, Eigen::MatrixXd sigma);
double uv_gamma(double shape, double rate);

#endif
