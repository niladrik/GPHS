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


#### Inverse gamma prior  ####

# the prior function
p <- function(theta){
  # length scale hyperparameter having discrete uniform prior over 20 points
  l = theta[[1]]

  # variance of the signal function having inverse gamma(1, 1) prior
  tau = theta[[2]]
  # IG(alpha, beta)
  alpha = 1
  beta = 1
  return( beta^alpha * tau^(-alpha - 1) * exp(-beta/ tau) / (gamma(alpha) * 20))
}

