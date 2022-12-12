## matern covariance function

matern = function(d , phi, nu = 3){ # previous value of nu = 3
  ifelse(d > 0, (sqrt(2 * nu) * d * phi)^nu / (2^(nu - 1) * gamma(nu)) *
           besselK(x = d * phi * sqrt(2 * nu), nu = nu), 1.0)
}


## covariance matrix

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

