#### Functions for future extension ####

#' Calculate the inverse and log(determinant) of symmetric matrix
#'
#' @param Sigma is a symmetric matrix whose inverse and determinant are to be calculated
#'
#' @return A list containing the inverse and the log(determinant) of the matrix
#'
#'
#' @examples
#' mat <- diag(c(2, 4))
#' inv_and_logdet(mat)

inv_and_logdet = function(Sigma){
  # compatibility check
  if(!isSymmetric.matrix(Sigma)){
    stop("Sigma should be symmetric")
  }
  A = chol(Sigma) ## Cholesky decomposition
  d = 2 * sum(log(diag(A))) ## determinant
  result = list()
  result[[1]] = chol2inv(A) ## inverse from Cholesky decomposition (X'X)^{-1}
  result[[2]] = d
  return(result)
}

## square exponential function:
sqExp <- function(d, phi){
  return(exp(-d^2 * phi / 2))
}


## function to generate data

f_data <- function(x){
  func <- 0
  for(j in 1:5000){
    func <- func + j^(-1.7) * sin(j) * cos(pi * (j - 0.5) * x)
  }
  return(func)
}


## Inverse gamma prior
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

## Uniform prior
p_uni <- function(theta, l.min, l.max){
  return(1/(l.min - l.max))
}


