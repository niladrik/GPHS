## This is the implementation of the Algorithm in paper by Nishihara, Adams and Murray (2014)

## Algorithm 1:

ESS <- function(x, mu, sigma, log.L){

  # choose ellipse
  v <- rmnorm(1, mu, sigma)
  u <- runif(1, 0, 1)

  # set log-likelihood threshold
  log.y <- log.L(x) + log(u)

}






















f <- function()

