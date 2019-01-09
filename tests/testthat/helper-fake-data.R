
# make fake data for linear model
fake_data_lm_proper <- function(beta, sigma, seed = 42){
  set.seed(seed)
  x <- matrix(rnorm(100), byrow = TRUE, ncol = 2)
  y <- x %*% beta + rnorm(2, sd = sigma)
  df <- data.frame(y, x)
  colnames(df) <- c("y", "x1", "x2")
  return(df)
}
