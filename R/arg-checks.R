

#' @import checkmate

# nocov start
check_args_lmm_improper <- function(data,
                                    tau_prior_shape,
                                    tau_prior_rate){
  checkmate::assert_data_frame(data, any.missing = FALSE)
  checkmate::assert_numeric(tau_prior_shape, any.missing = FALSE, len = 2)
  checkmate::assert_numeric(tau_prior_rate, any.missing = FALSE, len = 2)
}
# nocov end
