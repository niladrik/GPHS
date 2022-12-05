## matern covariance function

matern = function(d , phi, nu = 3){ # previous value of nu = 3
  ifelse(d > 0, (sqrt(2 * nu) * d * phi)^nu / (2^(nu - 1) * gamma(nu)) *
           besselK(x = d * phi * sqrt(2 * nu), nu = nu), 1.0)
}


# density of multivariate normal distribution
dmvnorm <- function(y, mu, sigma){
  p = length(y)
  # inverse and det of sigma
  result = inv_and_logdet(sigma)
  return(exp(-tcrossprod(y - mu, result[[1]] * (y - mu)) * 0.5 - result[[2]] - p * log(2*pi) / 2))
}


