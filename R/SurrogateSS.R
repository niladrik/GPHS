# This is implementation of algorithm 4 in Adams and Murray(2010)

#' Surrogate data model - Slice sampling
#'
#' @param theta hyperparameter of interest
#' @param f latent variable
#' @param sigma scale parameter
#' @param Sig covariance matrix of distribution of f
#' @param S covariance matrix of the noisy version of the true latent variable
#'
#' @return updated theta and f
#' @export
#'
#' @examples
SurrogateSS <- function(theta, f, sigma, Sig, S){
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
}
