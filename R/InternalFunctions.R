## matern covariance function

matern = function(d , phi, nu = 3){ # previous value of nu = 3
  ifelse(d > 0, (sqrt(2 * nu) * d * phi)^nu / (2^(nu - 1) * gamma(nu)) *
           besselK(x = d * phi * sqrt(2 * nu), nu = nu), 1.0)
}


## square exponential function:
sqExp <- function(d, phi){
  return(exp(-d^2 * phi / 2))
}

# density of multivariate normal distribution
dmvnorm_own <- function(y, mu = rep(0, length(y)), sigma = diag(1, length(y)), tau = 0.1){
  p = length(y)
  # inverse and det of sigma
  result = mySolve(tau * sigma)
  return(exp(-crossprod(y - mu, result[[1]] %*% (y - mu)) * 0.5 - result[[2]] - p * log(2*pi) / 2))
}

# covariance matrix

Sigma <- function(data, l){
  # number of columns
  n.data = ncol(data)

  # extracting x and y
  y = matrix(data[, 1], ncol = 1)
  x = matrix(data[, 2:n.data], ncol = n.data - 1)

  n = nrow(x)

  p_s = ncol(x) ## dimension of the space

  D = list()

  for (i in 1:p_s)
  {
    D[[i]] = as.matrix(dist(x[, i], diag = TRUE, upper = TRUE))
  }

  new_cov = matrix(1, n, n)
  for (i in 1:p_s)  new_cov = new_cov * matern(d = D[[i]], 1/l)

  return(new_cov)
}


# function to generate data

f_data <- function(x){
  func <- 0
  for(j in 1:5000){
    func <- func + j^(-1.7) * sin(j) * cos(pi * (j - 0.5) * x)
  }
  return(func)
}


#### Inverse gamma prior  ####

# the prior function
p_ig <- function(theta){
  # length scale hyperparameter having discrete uniform prior over 20 points
  l = theta[[1]]

  # variance of the signal function having inverse gamma(1, 1) prior
  tau = theta[[2]]
  # IG(alpha, beta)
  alpha = 1
  beta = 1
  return( beta^alpha * tau^(-alpha - 1) * exp(-beta/ tau) / (gamma(alpha) * 20))
}

#### Uniform prior ####
p_uni <- function(theta, l.min, l.max){
  return(1/(l.min - l.max))
}

## Inversion of matrix for symmetric matrix
mySolve <- function(A){
  eig_result = eigen(A, symmetric = TRUE)
  lam = eig_result$values
  V = as.matrix(eig_result$vectors)
  lam.new = lam
  lam.new[lam < 1e-11] = 1e-10
  lam.inv = 1/lam.new
  cholDeco = sqrt(lam.new) * V
  return(list(inv = tcrossprod(sqrt(lam.inv) * V), det = prod(lam.new), cholDeco = cholDeco))
}
