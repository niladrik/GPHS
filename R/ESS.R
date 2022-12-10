## Algorithm 1:

#' Elliptical Slice Sampling
#' @description
#' This is the implementation of the Elliptical Slice sampling algorithm as described in the paper by Nishihara, Adams and Murray (2014)
#' @param x current state
#' @param mu mean of Gaussian process
#' @param sigma covariance matrix of Gaussian process
#' @param log.L log-likelihood function
#' @param niter maximum number of iterations in one update, default is 100
#'
#' @return It will return the new state x'
#' @export
#'
#' @references
#' Nishihara, R., Murray, I., & Adams, R. P. (2014). Parallel MCMC with generalized elliptical slice sampling.
#' \emph{The Journal of Machine Learning Research, 15(1)}, 2087-2112.
#' \doi{10.48550/ARXIV.1210.7477}
#'
#' @examples
#' set.seed(12345)
#' x <- seq(1, 1e-3, -0.1)
#' ff <- x
#' y <- ff + rnorm(length(x), 0, sd(ff)/10)
#' training_data = cbind(y, x)
#' # conditional likelihood
#' log.L = function(y, mu = NULL, sigma = NULL){
#' p = length(y)
#' if(is.null(mu)){
#'   mu = rep(0, p)
#' }
#' if(is.null(sigma)){
#'    sigma = diag(1, p)
#' }
#' d = mvtnorm::dmvnorm(y, mu, sigma)
#' log.val = log(d)
#' return(log.val)
#' }
#' # initial value is randomly chosen to be 5
#' # here we consider f = log(l) ~ N(mu, sigma)
#' x.start = 5
#' # calculating the length of x.start
#' n = length(x.start)
#' ESS(x = x.start, mu = rep(0, n),
#'     sigma = diag(0.1, n), log.L = log.L, niter = 100)
#'
ESS <- function(x, mu, sigma, log.L, niter = 100){

  # some compatibility checks
  if(niter <= 0){
    stop("niter needs to be positive!")
  }
  if(any(sigma <= 0)){
    stop("scale parameter needs to be positive!")
  }
  if(is.null(x)){
    stop("initial state cannot be missing!")
  }

  # choose ellipse
  v <- MASS::mvrnorm(1, mu, sigma)
  u <- runif(1, 0, 1)

  # set log-likelihood threshold
  log.y <- log.L(x) + log(u)

  # drawing an initial proposal
  theta <- runif(1, 0, 2*pi)

  # defining a bracket
  theta.min = theta - 2*pi
  theta.max = theta

  # defining a flag variable for the while loop
  flag <- 1

  # counter
  count = 0

  # slice sampling the loop
  while(flag){

    # counting the number of iterations
    count = count + 1

    # calculating a new x in the proposed angle
    x.n = (x - mu) * cos(theta) + (v - mu) * sin(theta) + mu

    # checking if the proposed x is in the slice
    if(log.L(x.n) > log.y){
      flag = 0
      return(x.n)
    }

    # shrinking the bracket and proposing a new angle
    if(theta < 0){
      theta.min = theta
    } else{
      theta.max = theta
    }
    theta = runif(1, theta.min, theta.max)

    # checking the number of iterations till now
    if(count >= niter){
      print("Exceeded the maximum number of iterations!")
      flag = 0
      return(x.n)
    }
  }
}
