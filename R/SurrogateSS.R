# This is implementation of algorithm 4 in Adams and Murray(2010)

#' Surrogate data model - Slice sampling
#'
#' @param theta hyperparameter of interest
#' @param f latent variable
#' @param sigma scale parameter
#' @param Sig covariance matrix of distribution of f
#' @param S covariance matrix of the noisy version of the true latent variable
#' @param l conditional likelihood function
#' @param p distribution of the covariance hyperparameter
#'
#' @return updated theta and f
#' @export
#'
#' @examples
SurrogateSS <- function(theta, f, sigma, Sig, S, l, p){
  # draw surrogate data
  g = mvrnorm(1, f, S)

  # calculating R
  R = solve(solve(Sig) + solve(S))

  # calculating m
  m = R %*% solve(S) %*% g

  # Cholesky decomposition of R
  L = chol(R)

  # inverse of L
  Linv = solve(L)

  # compute the latent variable
  eta = Linv %*% (f - m)

  # create a bracket
  v = runif(0, sigma)

  theta.min = theta - v
  theta.max = theta + sigma

  u = runif(0, 1)

  # determining threshold
  y = u * l(f) * dmvnorm(g, 0, Sig + S) * p(theta)
  repeat{
    # drawing the proposal
    theta.p = runif(1, theta.min, theta.max)

    # compute function
    f.p = L * eta + m

    if(l(f.p) * dmvnorm(g, 0, Sig + S) * p(theta.p) > y){
      return(list(f = f.p, theta = theta.p))
    } else if(theta.p < theta){
      # shrinking the bracket minimum
      theta.min = theta.p
    }else{
      # shrinking the bracket maximum
      theta.max = theta.p
    }
  }


}


#' Calculate the inverse and log(determinant) of symmetric matrix
#'
#' @param Sigma is a symmetric matrix whose inverse and determinant are to be calculated
#'
#' @return A list containing the inverse and the log(determinant)
#'
#' @examples
#' inv_and_logdet(matrix(c(2,0,0,4), nrow = 2, ncol = 2))
#'
inv_and_logdet = function(Sigma)
{
  A = chol(Sigma) ## Cholesky decomposition
  d = 2 * sum(log(diag(A))) ## determinant
  result = list()
  result[[1]] = chol2inv(A) ## inverse from Cholesky decomposition (X'X)^{-1}
  result[[2]] = d
  return(result)
}

dmvnorm <- function(x, mu, sigma){
  # inverse and det of sigma
  result = inv_and_logdet(sigma)
  return(exp(-tcrossprod(y - mu, result[[1]] * (y - mu)) * 0.5 - result[[2]] - p * log(2*pi) / 2))
}

