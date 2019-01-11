#' Fit a Bayesian Generalized Linear Mixed Model With Improper Prior
#'
#' Fists a Generalized Linear Model with a flat prior on \eqn{\beta} and
#' improper Gamma priors on the precision for \eqn{u} and \eqn{e}.
#'
#' @param data A data frame containing the variables to be used in the model.
#' @param formula A formula object describing the model fit.
#' @param burnin Number of burn in iterations
#' @param iterations Number of MCMC iterations
#' @param thin Number of thinning iterations
#' @param lambda_prior_shape A length 2 vector containing the shape for the prior
#'  on \eqn{\lambda_e} and \eqn{\lambda_u}, in that order.
#' @param lambda_prior_rate A length 2 vector containing the rate for the prior
#'  on \eqn{\lambda_e} and \eqn{\lambda_u}, in that order.
#' @param start_theta Starting vector for \eqn{\theta = (\beta' u')'}.
#'  Default is the frequentist guess.
#'
#' @return A list containing MCMC samples from the posterior distribution of
#'  \eqn{\beta}, \eqn{u}, \eqn{\sigma_e}, and \eqn{\sigma_u}.
#'
#'
#' @import lme4
#' @import stats
#' @import checkmate
#'
#' @export

lmm_improper <- function(data,
                         formula,
                         burnin = 5000,
                         iterations = 5000,
                         thin = 1,
                         lambda_prior_shape,
                         lambda_prior_rate,
                         start_theta){

  checkmate::assert_data_frame(data, any.missing = FALSE)

  if (!inherits(formula, "formula"))
    stop("'formula' argument must be a formula.")

  check_args_lmm_improper(data,
                          tau_prior_shape,
                          tau_prior_rate)

  fit <- lme4::lmer(formula, data = data)
  x <- lme4::getME(fit, "X")
  z <- as.matrix( lme4::getME(fit, "Z") )
  rand_effect <- lme4::ranef(fit)[[1]]
  y <- lme4::getME(fit, "y")

  p <- ncol(x)
  q <- 1

  if ( missing(lambda_prior_shape) ){
    lambda_prior_shape <- c(0, -0.5)
  } else {
    checkmate::assert_numeric(lambda_prior_shape,
                              len = 2,
                              any.missing = FALSE)
  }

  if ( missing(lambda_prior_rate) ){
    lambda_prior_rate <- c(0, 0)
  } else {
    checkmate::assert_numeric(lambda_prior_rate,
                              len = 2,
                              any.missing = FALSE)
  }

  if ( missing(start_theta) ){
    beta <- lme4::fixef(fit)
    u <- rand_effect[ ,1]
    start_theta <- c(beta, u)
  } else {
    checkmate::assert_numeric(start_theta,
                              len = p + q,
                              any.missing = FALSE)
  }

  burnin <- checkmate::asInt(burnin, lower = 1)
  iterations <- checkmate::asInt(iterations, lower = 1)
  thin <- checkmate::asInt(thin, lower = 1)

  result <- lmm_improper_cpp(x,
                             z,
                             y,
                             lambda_prior_shape,
                             lambda_prior_rate,
                             iterations,
                             burnin,
                             thin,
                             start_theta)

  result$beta <- as.data.frame(result$beta)
  result$sigma <- as.data.frame(result$sigma)

  names(result$beta) <- stats::variable.names(x)
  names(result$u) <- paste0(names(u), "[", rownames(u), "]")
  names(result$sigma) <- c("sigma[e]", "sigma[u]")

  return(result)

}
