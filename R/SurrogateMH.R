## This is the implementation of Surrogate data Metropolis Hastings algorithm as given in Adams and Murray (2010)
#
#' Surrogate data model - Metropolis Hastings
#'
#' @param theta hyperparameter of interest
#' @param f latent variable
#' @param q proposal distribution
#' @param Sig covariance matrix of distribution of f
#' @param S covariance matrix of the noisy version of the true latent variable
#'
#' @return updated theta and f
#' @export
#'
#' @examples
SurrogateMH <- function(theta, f, q, Sig, S){
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

