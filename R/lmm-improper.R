#' Fit a linear mixed model with an improper prior.
#'
#' Fits a Bayesian linear mixed model with a flat prior on beta and improper
#'     gamma prior on the precision lambda for the fixed effects and error term.
#'     A standard reference prior for lambda is used by default. The formula
#'     syntax is the same as \code{\link[lme4:lmer]{lme4::lmer}} which is used
#'     to obtain the model matrices. Only a singleintercept is supported for
#'     the random effects.
#'
#' @section Warning:
#' We leave it to the user to ensure only a single intercept is specified
#'     in the model formula.
#'
#' @references Roman, J. C. and Hobert, J. P. (2012).
#'     Convergence analysis of the Gibbs sampler for Bayesian general linear
#'     mixed models with improper priors. Annals of Statistics, 40 2823–2849.
#'
#' @param data A data frame.
#' @param formula A formula describing the model fit. Passed to
#'     \code{\link[lme4:lmer]{lme4::lmer}} to construct model matrices.
#' @param burnin Number of burn in iterations
#' @param iterations Number of sampling iterations
#' @param thin Number of thinning iterations
#' @param lambda_prior_shape (Optional) Shape parameter for the prior on the
#'     precision of the error term and random effects, in that order. Defaults
#'     to \code{c(0, -0.5)}.
#' @param lambda_prior_rate (Optional) Rate parameter for the prior on the
#'     precision of the error term and random effects, in that order. Defaults
#'     to \code{c(0, 0)}.
#' @param start_theta (Optional) Starting vector for
#'     \eqn{\theta = (\beta' u')'}, the concatenation of the fixed and random
#'     effects coefficients. Defaults to the frequentist estimate.
#'
#' @return A list containing MCMC samples from the posterior distributions of
#'  the fixed effects, random effects, and standard deviations for the error
#'  term and fixed effects.
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

  fit <- lme4::lmer(formula, data = data)
  x <- lme4::getME(fit, "X")
  z <- as.matrix( lme4::getME(fit, "Z") )
  rand_effect <- lme4::ranef(fit)[[1]]
  y <- lme4::getME(fit, "y")

  p <- ncol(x)
  q <- ncol(z)

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
  result$u <- as.data.frame(result$u)
  result$sigma <- as.data.frame(result$sigma)

  names(result$beta) <- stats::variable.names(x)
  names(result$u) <- paste0(names(rand_effect), "[", rownames(rand_effect), "]")
  names(result$sigma) <- c("sigma[e]", "sigma[u]")

  return(result)

}
