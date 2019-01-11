context("lmm improper")

beta <- c(10, 15)
sigma <- c(2,2)
fake_data <- fake_data_lmm_improper(beta, sigma[1], sigma[2])
df <- fake_data$df
u <- fake_data$u

fit <- lmm_improper(df,
                    y ~ X1 + X2 - 1 + (1|group),
                    burnin = 1000,
                    iterations = 1000)

beta_guess <- colMeans(fit$beta)
u_guess <- colMeans(fit$u)
sigma_guess <- colMeans(fit$sigma)


test_that("estimates are close to true value", {
  expect_equivalent(beta_guess, beta, tolerance = 0.5)
  expect_equivalent(u_guess, u)
  expect_equivalent(sigma_guess, sigma)

})

