## This is the implementation of Surrogate data Metropolis Hastings algorithm as given in Adams and Murray (2010)
#
#' Surrogate data model - Metropolis Hastings
#'
#' @param theta hyperparameter of interest
#' @param f latent variable
#' @param data the data in hand
#' @param p the prior distribution
#' @param l the conditional likelihoof function of data conditioned on the latent variable f
#'
#' @return updated theta and f
#' @export
#'
#' @examples
#'
#'
#'
SurrogateMH <- function(theta, f, data, p, l){
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

  R.result = mySolve(S_inv + Sig_inv)
  R = R.result$inv

  # calculating m
  m = R %*% S_inv %*% g

  # Cholesky decomposition of R
  L = R.result$cholDeco

  # inverse of L
  Linv = mySolve(L)$inv

  # compute the latent variable
  eta = Linv %*% (f - m)

  # propose new theta
  # taking log-normal(log(theta), 1) distribution to be the proposal distribution
  theta.p = exp(rnorm(1, log(theta), 1))

  # compute the new Sigma matrix
  Sigma.p = Sigma(data, theta.p)

  # new R
  Rp_result = mySolve(mySolve(Sigma.p)$inv + S_inv)
  R.p = Rp_result$inv

  # calculating new m
  m.p = R.p %*% S_inv %*% g

  # compute new L matrix
  L.p = Rp_result$cholDeco # chol(R.p)

  # compute the function
  f.p = as.vector(L.p %*% eta + m.p)

  # generating the threshold randomly
  u = runif(1)

  # calculate the ratio
  lik.new = l(f.p) * dmvnorm_own(g, rep(0, length(g)), Sigma.p + S) * p(theta.p) * dnorm(log(theta.p), log(theta), 1) / theta.p
  lik.old = l(f) * dmvnorm_own(g, rep(0, length(g)), Sig + S) * p(theta) * dnorm(log(theta), log(theta.p), 1) / theta
  if(lik.new / lik.old > u){
    return(list(f = f.p, theta = theta.p))
  }else{
    return(list(f = f, theta = theta))
  }
}

