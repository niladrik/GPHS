## This is the implementation of the Algorithm in paper by Nishihara, Adams and Murray (2014)

## Algorithm 1:

#' Title Elliptical Slice Sampling
#'
#' @param x current state
#' @param mu mean of Gaussian process
#' @param sigma sd of Gaussian process
#' @param log.L log-likelihood function
#'
#' @return It will return the new state x'
#' @export
#'
#' @examples
#'
#'
#'
#'
ESS <- function(x, mu, sigma, log.L, niter){

  # choose ellipse
  v <- mvrnorm(1, mu, sigma)
  u <- runif(1, 0, 1)

  # set log-likelihood threshold
  log.y <- log.L(x) + log(u)

  # drawing an initial proposal
  theta <- runif(0, 2*pi)

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
