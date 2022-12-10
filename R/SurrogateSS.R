#' Slice sampling for Surrogate data model
#' @description
#' This is the implementation of the Slice sampling algorithm to update the hyperparameter
#'  and latent variable from the joint posterior in Surrogate data model as given in the paper by Adams and Murray (2010)
#' @param theta hyperparameter of interest
#' @param f latent variable
#' @param sigma scale parameter
#' @param l conditional likelihood function
#' @param p prior density function of the covariance hyperparameter
#' @param data the data in hand
#' @param niter maximum number of iterations in an update
#'
#' @return  Return a list containing
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
SurrogateSS <- function(theta, f, sigma, data, l, p, niter){
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
  if(niter < 0){
    stop("niter cannot be negative or 0")
  }


  # # checking if l is NULL then assigning the multivariate normal density to it
  # if(is.null(l)){
  #   l = dmvnorm_own
  # }

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
  L = R_return$sqrtDeco

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
    L.p = Rp_res$sqrtDeco

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
  # compatibility check
  if(!isSymmetric.matrix(Sigma)){
    stop("Sigma should be symmetric")
  }
  A = chol(Sigma) ## Cholesky decomposition
  d = 2 * sum(log(diag(A))) ## determinant
  result = list()
  result[[1]] = chol2inv(A) ## inverse from Cholesky decomposition (X'X)^{-1}
  result[[2]] = d
  return(result)
}

#' Density of multivariate normal distribution
#'
#' @param y vector of quantiles
#' @param mu mean vector, default is ```rep(0, length(y))```
#' @param sigma covariance matrix, default is ```diag(1, length(y))```
#' @param tau constant multiplier (positive) for the covariance matrix, default is $0.1$
#'
#' @return It returns the density function for the multivariate normal distribution
#' @export
#'
#' @examples
#' dmvnorm_own(y = c(0,0))
#' dmvnorm_own(y = c(0,0), mu = c(1,1))
#' dmvnorm_own(y = c(0,0), mu = c(1,1), sigma = diag(2, 2), tau = 1)
dmvnorm_own <- function(y, mu = rep(0, length(y)), sigma = diag(1, length(y)), tau = 0.1){
  p = length(y)
  # compatibility checks
  if(length(mu) != p){
    stop("y and mu are not compatible")
  }
  if(!isSymmetric(sigma, tol = 1e-5)){
    stop("sigma needs to be a symmetric matrix")
  }
  if(nrow(sigma) != p){
    stop("y and sigma are not compatible")
  }
  if(tau <= 0){
    stop("tau needs to be positive")
  }

  # inverse and det of sigma
  result = mySolve(tau * sigma)
  return(exp(-crossprod(y - mu, result[[1]] %*% (y - mu)) * 0.5 - result[[2]] - p * log(2*pi) / 2))
}

#' Inversion of a symmetric matrix
#'
#' @param A a symmetric matrix
#'
#' @return A list containing the inverse, determinant and the square root matrix $Q$ of $A$, that is, it holds $QQ=A$
#'
#' @export
#'
#' @examples
#' A = diag(2, 2)
#' mySolve(A)
#'
mySolve <- function(A){
  # compatibility check
  if(!isSymmetric.matrix(A)){
    stop("A must be symmetric")
  }
  # eigen decomposition
  eig_result = eigen(A, symmetric = TRUE)

  # eigen values
  lam = eig_result$values

  # eigen vector matrix
  V = eig_result$vectors

  # fixing the tolerance for eigen values
  tol = 1e-7

  # extracting the eigen values which are above the tolerance level
  pos = lam > tol

  # calculating the inverse
  inv = tcrossprod(V[, pos] %*% diag(1/lam[pos]), V[, pos])

  # determinant
  det = prod(lam[pos])

  # square root decomposition
  sqrt.mat = tcrossprod(V[, pos] %*% diag(sqrt(lam[pos])), V[, pos])

  return(list(inv = inv, det = det, sqrtDeco = sqrt.mat))
}
