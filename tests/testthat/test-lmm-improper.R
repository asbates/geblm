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
  expect_equivalent(u_guess, u, tolerance = 0.5)
  expect_equivalent(sigma_guess, sigma, tolerance = 0.5)

})

test_that("output has correct format", {
  expect_named(fit$beta, c("X1", "X2"))
  expect_named(fit$u, paste0("(Intercept)",
                             "[",
                             levels(as.factor(df$group)),
                             "]"))
  expect_named(fit$sigma, c("sigma[e]", "sigma[u]"))

  expect_is(fit, "list")
  expect_is(fit$beta, "data.frame")
  expect_is(fit$u, "data.frame")
  expect_is(fit$sigma, "data.frame")

  expect_length(fit$beta$X1, 1000)
  expect_length(fit$u[,1], 1000)
  expect_length(fit$sigma[,1], 1000)
})


test_that("output doesn't have missing values", {
  beta_has_na <- any( apply(fit$beta, 2, is.na) )
  u_has_na <- any( apply(fit$u, 2, is.na) )
  sigma_has_na <- any( apply(fit$sigma, 2, is.na) )
  expect_false(beta_has_na)
  expect_false(u_has_na)
  expect_false(sigma_has_na)
})

test_that("convergence check outputs message on success", {
  expect_output(lmm_improper(df,
                              y ~ X1 + X2 - 1 + (1|group),
                              burnin = 1,
                              iterations = 1),
                 "Conditions for geometric convergence are satisfied.")
})

test_that("convergence check outputs warning on failure", {
  # ((au < bu) and bu == 0) or (bu > 0) )
  a <- c(0,0)
  b <- c(0,0)
  expect_warning(lmm_improper(df,
                              y ~ X1 + X2 - 1 + (1|group),
                              burnin = 1,
                              iterations = 1,
                              lambda_prior_shape = a,
                              lambda_prior_rate = b),
                 "Conditions for geometric convergence are not satisfied.")
})
