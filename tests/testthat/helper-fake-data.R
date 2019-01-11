
# make fake data for linear model
fake_data_lm_proper <- function(beta, sigma, seed = 42){
  set.seed(seed)
  x <- matrix(rnorm(100), byrow = TRUE, ncol = 2)
  y <- x %*% beta + rnorm(2, sd = sigma)
  df <- data.frame(y, x)
  colnames(df) <- c("y", "x1", "x2")
  return(df)
}


# make fake data for linear mixed model with improper prior

fake_data_lmm_improper <- function(beta, sigma_e, sigma_u, seed = 15){

  set.seed(seed)
  n <- 50
  p <- length(beta)
  q <- 4

  # design matrices
  x <- matrix(rnorm(p * n), ncol = p)
  z <- matrix(nrow = n, ncol = q)
  g1 <- 10
  g2 <- 10
  g3 <- 10
  g4 <- 20
  z[ , 1] <- c(rep(1, g1), rep(0, g2 + g3 + g4))
  z[ , 2] <- c(rep(0, g1), rep(1, g2), rep(0, g3 + g4))
  z[ , 3] <- c(rep(0, g1 + g2), rep(1, g3), rep(0, g4))
  z[ , 4] <- c(rep(0, g1 + g2 + g3), rep(1, g4))


  # random effects
  u <- rnorm(q, sd = sigma_u)

  # errors
  e <- rnorm(n, sd = sigma_e)

  # response
  y <- x %*% beta + z %*% u + e

  # put it all together
  df <- data.frame(y, x,
                   group = c(rep("group1", g1),
                             rep("group2", g2),
                             rep("group3", g3),
                             rep("group4", g4)))

  list(df = df, u = u)
}
