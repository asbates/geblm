#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
bool conv_check_lmm_improper(Eigen::MatrixXd x,
                              Eigen::VectorXd y,
                              Eigen::MatrixXd z,
                              Eigen::VectorXd a,
                              Eigen::VectorXd b){

  int n = x.rows();
  int p = x.cols();
  int q = z.cols();

  double ae = a(0);
  double au = a(1);
  double be = b(0);
  double bu = b(1);

  Eigen::MatrixXd xt = x.transpose();
  Eigen::MatrixXd ident_n = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd px_perp = ident_n - x * (xt * x).inverse() * xt;

  Eigen::MatrixXd w(n, p + q);
  w.leftCols(p) = x;
  w.rightCols(q) = z;
  Eigen::MatrixXd wtw = w.transpose() * w;
  Eigen::MatrixXd w_pinv = wtw.completeOrthogonalDecomposition().pseudoInverse();
  Eigen::MatrixXd theta_hat = w_pinv * w.transpose() * y;
  double diff_norm = (y - w * theta_hat).squaredNorm();

  int t = (px_perp * z).colPivHouseholderQr().rank();

  // do all the convergence checks pass?
  // if any conditions are not met, pass_check is set to false
  bool pass_check = true;

  if ( not ( ((au < bu) and bu == 0) or (bu > 0) ) ) pass_check = false;

  if ( not(t + 2*au > 2) ) pass_check = false;

  if ( not (n + 2*ae > p + t + 2) ) pass_check = false;

  if ( not (2*be + diff_norm > 0) ) pass_check = false;

  if (pass_check){

    Rcout << "Conditions for geometric convergence are satisfied." << std::endl;


  } else {

    warning("Conditions for geometric convergence are not satisfied.");

  }

  return pass_check;

}


bool conv_check_lmm_proper(Eigen::MatrixXd z,
                            Eigen::VectorXd a){
  double ae = a(0);
  double au = a(1);

  double n = z.rows();

  int rank_z = z.colPivHouseholderQr().rank();

  bool pass_strong_condition1 = true;
  bool pass_strong_condition2 = true;

  bool pass_first_checks = true;
  bool pass_second_checks = true;

  if ( not (n > rank_z + 2) ) pass_strong_condition1 = false;

  if ( not (rank_z > 2) ) pass_strong_condition2 = false;

  if (pass_strong_condition1 and pass_strong_condition2){
    Rcout << "Conditions for geometric convergence are satisfied." << std::endl;
  } else {

    pass_first_checks = false;

    bool pass_weak_condition1 = true;
    bool pass_weak_condition2 = true;
    Eigen::VectorXd ae_compare(2);
    ae_compare << 1.0 - 0.5 * (n - rank_z), 0.0;

    Eigen::VectorXd au_compare(2);
    au_compare << 1.0 - 0.5 * rank_z, 0.0;

    if ( not(ae > ae_compare.maxCoeff()) ) pass_weak_condition1 = false;

    if ( not(au > au_compare.maxCoeff()) ) pass_weak_condition2 = false;

    if (pass_weak_condition1 and pass_weak_condition2){

      Rcout << "Conditions for geometric convergence are satisfied." << std::endl;

    } else{

      pass_second_checks = false;
      warning("Conditions for geometric convergence are not satisfied.");
    }

  }

  return pass_first_checks or pass_second_checks;
}
