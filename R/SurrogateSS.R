# This is implementation of algorithm 4 in Adams and Murray(2010)

#' Surrogate data model - Slice sampling
#'
#' @param theta hyperparameter of interest
#' @param f latent variable
#' @param sigma scale parameter
#' @param l conditional likelihood function
#' @param p prior distribution of the covariance hyperparameter
#' @param data the data in hand
#' @param niter maximum number of iterations in an update
#'
#' @return updated theta and f
#' @export
#'
#' @examples
SurrogateSS <- function(theta, f, sigma, data, l = NULL, p, niter){
  # checking if l is NULL then assigning the multivariate normal density to it
  if(is.null(l)){
    l = dmvnorm_own
  }
  # fixing S by hand to a constant
  alpha = 0.1
  S = alpha * diag(1, nrow = length(f))

  # calculating sigma
  Sig = Sigma(data, theta)

  # draw surrogate data
  g = mvrnorm(1, f, S)

  # calculating R
  # invert sigma matrix
  Sig_inv = mySolve(Sig)$inv

  # invert S matrix
  S_inv = diag(1/alpha, nrow = length(f))

  R_return = mySolve(Sig_inv + S_inv)

  R = R_return$inv

  # calculating m
  m = R %*% S_inv %*% g

  # Cholesky decomposition of R
  L = R_return$cholDeco

  # inverse of L
  Linv = mySolve(L)$inv

  # compute the latent variable
  eta = Linv %*% (f - m)

  # create a bracket
  v = runif(1, 0, sigma)

  theta.min = theta - v
  theta.max = theta + sigma

  u = runif(1, 0, 1)

  # determining threshold
  y = u * l(f) * dmvnorm_own(g, rep(0, length(g)), Sig + S, 1) * p(theta) #, l.min, l.max)

  # count the number of iterations
  count = 0

  repeat{
    # counting the number of iterations
    count = count + 1

    # drawing the proposal
    theta.p = runif(1, theta.min, theta.max)

    # compute the new Sigma matrix
    Sigma.p = Sigma(data, theta.p)

    # new R
    Rp_res = mySolve(mySolve(Sigma.p)$inv + S_inv)

    R.p = Rp_res$inv
    # calculating new m
    m.p = R.p %*% S_inv %*% g

    # compute new L matrix
    L.p = Rp_res$cholDeco

    # compute the function
    f.p = as.vector(L.p %*% eta + m.p)

    if(l(f.p) * dmvnorm_own(g, rep(0, length(g)), Sigma.p + S, 1) * p(theta.p) > y){
      return(list(f = f.p, theta = theta.p))
    }else if(theta.p < theta){
      # shrinking the bracket minimum
      theta.min = theta.p
    }else{
      # shrinking the bracket maximum
      theta.max = theta.p
    }

    # checking the number of iterations till now
    if(count >= niter){
      return(list(f = f.p, theta = theta.p))
    }
  }
  print(count)

}


#' Calculate the inverse and log(determinant) of symmetric matrix
#'
#' @param Sigma is a symmetric matrix whose inverse and determinant are to be calculated
#'
#' @return A list containing the inverse and the log(determinant) of the matrix
#'
#' @export
#'
#' @examples
#' mat <- diag(c(2, 4))
#' inv_and_logdet(mat)

inv_and_logdet = function(Sigma){
  A = chol(Sigma) ## Cholesky decomposition
  d = 2 * sum(log(diag(A))) ## determinant
  result = list()
  result[[1]] = chol2inv(A) ## inverse from Cholesky decomposition (X'X)^{-1}
  result[[2]] = d
  return(result)
}

# density of multivariate normal distribution
#' Title
#'
#' @param y vector of quantiles
#' @param mu mean vector, default is ```rep(0, length(y))```
#' @param sigma covariance matrix, default is ```diag(1, length(y))```
#' @param tau constant multiplier for the covariance matrix
#'
#' @return It returns the density function for the multivariate normal distribution
#' @export
#'
#' @examples
#' dmvnorm_own(x=c(0,0))
#' dmvnorm_own(x=c(0,0), mean=c(1,1))
dmvnorm_own <- function(y, mu = rep(0, length(y)), sigma = diag(1, length(y)), tau = 0.1){
  p = length(y)
  # inverse and det of sigma
  result = mySolve(tau * sigma)
  return(exp(-crossprod(y - mu, result[[1]] %*% (y - mu)) * 0.5 - result[[2]] - p * log(2*pi) / 2))
}

