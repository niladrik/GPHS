% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SurrogateMH.R
\name{SurrogateMH}
\alias{SurrogateMH}
\title{Metropolis Hastings update for Surrogate data model}
\usage{
SurrogateMH(theta, f, data, p, l)
}
\arguments{
\item{theta}{a scalar, hyperparameter of interest}

\item{f}{a scalar or vector, latent variable}

\item{data}{a matrix or data frame, the dataset in hand}

\item{p}{a function, the prior distribution}

\item{l}{a function, the conditional likelihood function of data conditioned on the latent variable f}
}
\value{
Return a list containing
\item{theta}{updated value of theta}
\item{f}{updated value of f}
}
\description{
This is the implementation of the Metropolis-Hastings algorithm to update the hyperparameter
and latent variable from the joint posterior in Surrogate data model as given in the paper by Adams and Murray (2010)
}
\examples{
set.seed(12345)
# generate the data
x <- seq(1, 1e-3, -0.1)
ff <- x
y <- ff + rnorm(length(x), 0, sd(ff)/10)
training_data = cbind(y, x)
# prior function
prior <- function(x, mu = 0, sigma = 1){
  return(dnorm(log(x), mu, sigma))
}
# starting points
theta.0 = runif(1)
f.0 = training_data[, 1] / 2

# We choose the conditional distribution of the data|f as N(f, tau*I)
l <- function(x, mu = rep(0, length(x)), sigma = diag(1, length(x)), tau = 0.1){
  return(dmvnorm_own(x, mu, sigma, tau))
}

SurrogateMH(theta = theta.0, f = f.0,
             data = training_data, p = prior, l = l)

}
\references{
Murray, I., & Adams, R. P. (2010).
Slice sampling covariance hyperparameters of latent Gaussian models.
\emph{Advances in neural information processing systems, 23.}
\doi{10.48550/ARXIV.1006.0868}
}
