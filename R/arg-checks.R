

#' @import checkmate

check_args_lm <- function(data,
                          beta_prior_mean,
                          beta_prior_cov,
                          tau_prior_shape,
                          tau_prior_rate){
  checkmate::assert_data_frame(data, any.missing = FALSE)
  checkmate::assert_numeric(beta_prior_mean, any.missing = FALSE)
  checkmate::assert_matrix(beta_prior_cov, mode = "numeric", any.missing = FALSE)
  checkmate::assert_number(tau_prior_shape, lower = 0)
  checkmate::assert_number(tau_prior_rate, lower = 0)
}


check_args_lmm_improper <- function(data,
                                    tau_prior_shape,
                                    tau_prior_rate){
  checkmate::assert_data_frame(data, any.missing = FALSE)
  checkmate::assert_numeric(tau_prior_shape, any.missing = FALSE, len = 2)
  checkmate::assert_numeric(tau_prior_rate, any.missing = FALSE, len = 2)
}
