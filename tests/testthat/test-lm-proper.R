context("test-lm-proper")

beta <- c(10, 15)
sigma <- 2
df <- fake_data_lm_proper(beta, sigma)

fit <- lm_proper(df,
                 y ~ -1 + x1 + x2,
                 beta_prior_mean = c(0,0),
                 beta_prior_cov = diag(2) * 100,
                 tau_prior_shape = 0.001,
                 tau_prior_rate = 0.001,
                 iterations = 1000,
                 burnin = 1000
)

bayes_guess <- colMeans(fit)
beta_guess <- bayes_guess[1:2]
sigma_guess <- bayes_guess[3]

test_that("estimates are close to true values", {
  expect_equivalent(beta_guess, beta, tolerance = 0.5)
  expect_equivalent(sigma_guess, sigma, tolerance = 0.5)
})


test_that("output is a named data frame of correct dimensions",{
  expect_named(fit, c("x1", "x2", "sigma"))
  expect_is(fit, "data.frame")
  expect_length(fit$x1, 1000)
})

test_that("output doesn't have missing values", {
  has_na <- any(apply(fit, 2, is.na))
  expect_false(has_na)
})


