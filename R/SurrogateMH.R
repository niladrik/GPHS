#' Metropolis Hastings update for Surrogate data model
#' @description
#' This is the implementation of the Metropolis-Hastings algorithm to update the hyperparameter
#'  and latent variable from the joint posterior in Surrogate data model as given in the paper by Adams and Murray (2010)
#' @param theta a scalar, hyperparameter of interest
#' @param f a scalar or vector, latent variable
#' @param data a matrix or data frame, the dataset in hand
#' @param p a function, the prior distribution
#' @param l a function, the conditional likelihood function of data conditioned on the latent variable f
#'
#' @return
#' Return a list containing
#' \item{theta}{updated value of theta}
#' \item{f}{updated value of f}
#'
#' @export
#' @references
#' Murray, I., & Adams, R. P. (2010).
#'  Slice sampling covariance hyperparameters of latent Gaussian models.
#'   \emph{Advances in neural information processing systems, 23.}
#'   \doi{10.48550/ARXIV.1006.0868}
#' @examples
#' set.seed(12345)
#' # generate the data
#' x <- seq(1, 1e-3, -0.1)
#' ff <- x
#' y <- ff + rnorm(length(x), 0, sd(ff)/10)
#' training_data = cbind(y, x)
#' # prior function
#' prior <- function(x, mu = 0, sigma = 1){
#'   return(dnorm(log(x), mu, sigma))
#' }
#' # starting points
#' theta.0 = runif(1)
#' f.0 = training_data[, 1] / 2
#'
#' # We choose the conditional distribution of the data|f as N(f, tau*I)
#' l <- function(x, mu = rep(0, length(x)), sigma = diag(1, length(x)), tau = 0.1){
#'   return(dmvnorm_own(x, mu, sigma, tau))
#' }
#'
#' SurrogateMH(theta = theta.0, f = f.0,
#'              data = training_data, p = prior, l = l)
#'
SurrogateMH <- function(theta, f, data, p, l){
  # compatibility checks
  if(is.null(theta)){
    stop("theta cannot be null")
  }
  if(is.null(f)){
    stop("f cannot be null")
  }
  if(is.null(data)){
    stop("data cannot be null")
  }
  if(length(theta) != 1){
    stop("theta needs to be a scalar")
  }
  if(length(f) != nrow(data)){
    stop("dimension of f and data are not compatible")
  }
  if(!is.function(p)){
    stop("p needs to be a function")
  }
  if(!is.function(l)){
    stop("l needs to be a function")
  }
  if(is.data.frame(data)){
    data = as.matrix(data)
  }
  if(!is.matrix(data)){
    stop("data should be a matrix or data frame")
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

  R.result = mySolve(S_inv + Sig_inv)
  R = R.result$inv

  # calculating m
  m = R %*% S_inv %*% g

  # Cholesky decomposition of R
  L = R.result$sqrtDeco

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
  L.p = Rp_result$sqrtDeco # chol(R.p)

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

