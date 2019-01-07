#' Fit a Bayesian linear model with Normal-Gamma prior via Gibbs Sampling.
#'
#' Fits a Bayesian linear regression model with a Normal prior on beta and a
#' Gamma prior on the precision tau. Model fitting is done via a Gibbs sampler.
#'
#'
#' @param data A data frame containing the variables to be used in the model.
#' @param formula A formula object describing the model fit.
#' @param beta_prior_mean Mean vector for the prior on \eqn{\beta}.
#' @param beta_prior_cov Covariance matrix for the prior on \eqn{\beta}.
#' @param tau_prior_shape Shape parameter for the prior on \eqn{\tau}.
#' @param tau_prior_rate Rate parameter for the prior on \eqn{\tau}.
#' @param iterations Number of MCMC iterations
#' @param burnin Number of burn in iterations
#' @param thin Number of thinning iterations
#' @param start_beta Beta starting vector. The default is \eqn{\beta}-hat.
#'
#' @return A data frame containing MCMC samples from the posterior distribution.
#'
#'
#' @import stats
#' @import checkmate
#'
#' @export
lm_proper <- function(data,
                        formula,
                        beta_prior_mean,
                        beta_prior_cov,
                        tau_prior_shape,
                        tau_prior_rate,
                        iterations = 5000,
                        burnin = 5000,
                        thin = 1,
                        start_beta = NULL){
  check_args_lm(data,
                beta_prior_mean,
                beta_prior_cov,
                tau_prior_shape,
                tau_prior_rate)

  iterations <- checkmate::asInt(iterations, lower = 1)
  burnin <- checkmate::asInt(burnin, lower = 1)
  thin <- checkmate::asInt(thin, lower = 1)

  fit <- stats::lm(formula, data, x = TRUE)
  x <- fit$x
  response <- all.vars(formula)[1]
  predictors <- stats::variable.names(fit)
  y <- data[, response]

  if (is.null(start_beta)){
    start_beta <- fit$coefficients
  } else {
    checkmate::assert_numeric(start_beta,
                              any.missing = FALSE,
                              len = dim(x)[1])
  }

  cpp_result <- lm_proper_cpp(x,
                            y,
                            beta_prior_mean,
                            beta_prior_cov,
                            tau_prior_shape,
                            tau_prior_rate,
                            iterations,
                            burnin,
                            thin,
                            start_beta)

  result <- as.data.frame(cpp_result)
  names(result) <- c(predictors, "sigma")
  return(result)
}
