## This is the implementation of Surrogate data Metropolis Hastings algorithm as given in Adams and Murray (2010)
#
#' Surrogate data model - Metropolis Hastings
#'
#' @param theta hyperparameter of interest
#' @param data the data in hand
#' @param p the prior distribution
#' @param f latent variable
#'
#' @return updated theta and f
#' @export
#'
#' @examples
#'
#'
#'
SurrogateMH <- function(theta, f, data, p){
  # fixing S by hand to a constant
  alpha = 0.1
  S = alpha * diag(1, nrow = length(f))

  # calculating sigma
  Sig = Sigma(data, theta)

  # draw surrogate data
  g = mvrnorm(1, f, S)

  # calculating R
  # invert sigma matrix
  Sig_inv = solve(Sig)

  # invert S matrix
  S_inv = diag(1/alpha, nrow = length(f))

  R = solve(S_inv + Sig_inv)

  # calculating m
  m = R %*% solve(S) %*% g

  # Cholesky decomposition of R
  L = chol(R)

  # inverse of L
  Linv = solve(L)

  # compute the latent variable
  eta = Linv %*% (f - m)

  # propose new theta
  # taking log-normal(log(theta), 1) distribution to be the proposal distribution
  theta.p = exp(rnorm(1, log(theta), 1))

  # if(theta >= 1){
  #   # taking the proposal distribuion to be Uniform(theta-1, theta+1)
  #   theta.p = runif(1, theta - 1, theta + 1)
  # }else{
  #   theta.p = 1
  # }

  # compute the function
  f.p = as.vector(L %*% eta + m)
  u = runif(1)

  # calculating new Sigma
  Sig.p = Sigma(data, theta.p)

  # calculate the ratio
  lik.new = l(f.p) * dmvnorm(g, rep(0, length(g)), Sig.p + S) * p(theta.p) * dnorm(log(theta.p), log(theta), 1) / theta.p
  lik.old = l(f) * dmvnorm(g, rep(0, length(g)), Sig + S) * p(theta) * dnorm(log(theta), log(theta.p), 1) / theta
  if(lik.new / lik.old > u){
    return(list(f = f.p, theta = theta.p))
  }else{
    return(list(f = f, theta = theta))
  }
}

