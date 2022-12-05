# This is implementation of algorithm 4 in Adams and Murray(2010)

#' Surrogate data model - Slice sampling
#'
#' @param theta hyperparameter of interest
#' @param f latent variable
#' @param sigma scale parameter
#' @param l conditional likelihood function
#' @param p distribution of the covariance hyperparameter
#' @param data the data in hand
#'
#' @return updated theta and f
#' @export
#'
#' @examples
SurrogateSS <- function(theta, f, sigma, data, l, p){
  # fixing S by hand to a constant
  alpha = 0.1
  S = alpha * diag(1, nrow = length(f))

  # calculating sigma
  sigma = Sigma_mat(data, theta)

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

    # updating Sigma
    Sig.p = Sigma_mat(data, theta.p)

    # compute function
    f.p = L * eta + m

    if(l(f.p) * dmvnorm(g, 0, Sig.p + S) * p(theta.p) > y){
      return(list(f = f.p, theta = theta.p))
    }else if(theta.p < theta){
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
#' @return A list containing the inverse and the log(determinant) of the matrix
#'
#' @export
#'
#' @examples
#' mat <- diag(c(2, 4))
#' inv_and_logdet(mat)

inv_and_logdet = function(Sigma)
{
  A = chol(Sigma) ## Cholesky decomposition
  d = 2 * sum(log(diag(A))) ## determinant
  result = list()
  result[[1]] = chol2inv(A) ## inverse from Cholesky decomposition (X'X)^{-1}
  result[[2]] = d
  return(result)
}

