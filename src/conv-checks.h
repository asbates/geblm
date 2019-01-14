// This is a c++ header file. It's needed to be able to call the functions
// in conv-checks.cpp in another file. e.g. glmm.cpp

#ifndef CONV_CHECKS_H
#define CONV_CHECKS_H

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;

void conv_check_lmm_improper(Eigen::MatrixXd x,
                              Eigen::VectorXd y,
                              Eigen::MatrixXd z,
                              Eigen::VectorXd a,
                              Eigen::VectorXd b);

void conv_check_lmm_proper(Eigen::MatrixXd z,
                           Eigen::VectorXd a);

#endif
