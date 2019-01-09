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

beta_guess <- colMeans(fit$beta)
sigma_guess <- colMeans(fit$sigma)

test_that("estimates are close to true values", {
  expect_equivalent(beta_guess, beta, tolerance = 0.5)
  expect_equivalent(sigma_guess, sigma, tolerance = 0.5)
})


test_that("output has correct format",{

  expect_named(fit$beta, c("x1", "x2"))
  expect_named(fit$sigma, "sigma")

  expect_is(fit, "list")
  expect_is(fit$beta, "data.frame")
  expect_is(fit$sigma, "data.frame")

  expect_length(fit$beta$x1, 1000)
  expect_length(fit$sigma$sigma, 1000)
})

test_that("output doesn't have missing values", {
  beta_has_na <- any( apply(fit$beta, 2, is.na) )
  sigma_has_na <- any( apply(fit$sigma, 2, is.na) )
  expect_false(beta_has_na)
  expect_false(sigma_has_na)
})


