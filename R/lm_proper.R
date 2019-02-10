#' Fit a linear model with a proper prior.
#'
#' Fits a Bayesian linear regression model with a normal prior on beta and a
#'     gamma prior on the precision tau. 'Wide' priors are used by default.
#'
#'
#' @param data A data frame.
#' @param formula A formula describing the model fit. Passed to
#'     \code{stats::lm()} to construct the model matrix.
#' @param burnin Number of burn in iterations.
#' @param iterations Number of sampling iterations.
#' @param thin Number of thinning iterations.
#' @param beta_prior_mean (Optional) Mean vector for the prior on beta. Defaults
#'     to zero.
#' @param beta_prior_cov (Optional) Covariance matrix for the prior on beta.
#'     Defaults to \code{diag(p) * 100} where p is the number of predictors.
#' @param tau_prior_shape (Optional) Shape parameter for the prior on tau.
#'     Defaults to 0.001.
#' @param tau_prior_rate (Optional) Rate parameter for the prior on tau.
#'     Defaults to 0.001.
#' @param start_beta (Optional) Beta starting vector. Defaults to beta hat,
#'     the frequentist estimate.
#'
#' @return An object of class \code{geblm} containing samples from the
#'     posterior distributions of beta and sigma.
#'
#' @import stats
#' @import checkmate
#'
#' @examples
#' fit <- lm_proper(mtcars, mpg ~ wt)
#'
#' @export
lm_proper <- function(data,
                      formula,
                      burnin = 5000,
                      iterations = 5000,
                      thin = 1,
                      beta_prior_mean,
                      beta_prior_cov,
                      tau_prior_shape,
                      tau_prior_rate,
                      start_beta){

  checkmate::assert_data_frame(data, any.missing = FALSE)

  if (!inherits(formula, "formula"))
    stop("'formula' argument must be a formula.")


  fit <- stats::lm(formula, data, x = TRUE, y = TRUE)
  x <- fit$x
  y <- fit$y
  p <- ncol(x)

  if (missing(beta_prior_mean)){
    beta_prior_mean <- rep(0, p)
  } else {
    checkmate::assert_numeric(beta_prior_mean,
                              len = p,
                              any.missing = FALSE)
  }

  if (missing(beta_prior_cov)){
    beta_prior_cov <- diag(p) * 100
  } else {
    checkmate::assert_matrix(beta_prior_cov,
                             mode = "numeric",
                             nrows = p,
                             ncols = p,
                             any.missing = FALSE)
  }

  if (missing(tau_prior_shape)){
    tau_prior_shape <- 0.001
  } else {
    checkmate::assert_number(tau_prior_shape, lower = 0)
  }

  if (missing(tau_prior_rate)){
    tau_prior_rate <- 0.001
  } else {
    checkmate::assert_number(tau_prior_rate, lower = 0)
  }

  if (missing(start_beta)){
    start_beta <- fit$coefficients
  } else {
    checkmate::assert_numeric(start_beta,
                              len = p,
                              any.missing = FALSE)
  }

  iterations <- checkmate::asInt(iterations, lower = 1)
  burnin <- checkmate::asInt(burnin, lower = 1)
  thin <- checkmate::asInt(thin, lower = 1)

  result <- lm_proper_cpp(x,
                          y,
                          beta_prior_mean,
                          beta_prior_cov,
                          tau_prior_shape,
                          tau_prior_rate,
                          iterations,
                          burnin,
                          thin,
                          start_beta)

  result$beta <- as.data.frame(result$beta)
  result$sigma <- as.data.frame(result$sigma)
  names(result$beta) <- stats::variable.names(fit)
  names(result$sigma) <- "sigma"

  result$conv_checks <- TRUE
  result$model_type <- "lm_proper"
  result$x <- x

  class(result) <- c("geblm", class(result))

  return(result)
}
