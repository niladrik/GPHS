
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- # GPHS -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/niladrik/GPHS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/niladrik/GPHS/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## GPHS: Hyperparameter selection in Gaussian Process models

The R package ‘GPHS’ provides the users with various algorithms for MCMC
update of the hyperparameters in covariance function in Gaussian process
models. It includes the Elliptical Slice Sampling algorithm as discussed
in [Nishihara, Murray and Adams
(2014)](https://www.cs.princeton.edu/~rpa/pubs/nishihara2014generalized.pdf)
and Metropolis-Hastings algorithm and Slice sampling algorithm for the
Surrogate data model as discussed in [Murray et. al
(2010)](https://arxiv.org/abs/1001.0175). In addition to these
functions, the package also provides the users with a function to invert
a matrix using (1) Cholesky decomposition and (2) Eigen decomposition,
and a function to calculate the density of a multivariate normal
distribution.

### Installation

You can install the most recent version of ‘GPHS’ package from
[GitHub](https://github.com/niladrik/GPHS) using the following commands:

``` r
# Plain installation
devtools::install_github("niladrik/GPHS")
```

``` r
# For installation with vignette
devtools::install_github("niladrik/GPHS", build_vignettes = TRUE)
```

### Examples

Here we would give an example of the function `ESS`.  

We define the log likelihood function

``` r
library(GPHS)
set.seed(12345)

# defining the log likelihood function
logl = function(y, mu = rep(0, length(y)), sigma = diag(1, length(y))){
d = mvtnorm::dmvnorm(y, mu, sigma)
log.val = log(d)
return(log.val)
}
```

Then we fix some random starting points for the parameter of interest:
`x`, and we run the `ESS` update 100 times.

``` r
n = 100
## vector of generated x
x_vec = vector("numeric", n + 1)
## starting point
x_vec[1] = 5
## dimension of current state
p = length(x_vec[1])
## ESS update
for(i in 1:n){
  x_vec[i+1] = ESS(x = x_vec[i], mu = rep(0, p), sigma = diag(1, p), log.L = logl, niter = 100)
}
```

The trace plot of the samples drawn is:
<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" /> For
the illustration of the other functions, please refer to the vignette.

### References

Murray, I., & Adams, R. P. (2010). Slice sampling covariance
hyperparameters of latent Gaussian models. *Advances in neural
information processing systems, 23.*
<https://doi.org/10.48550/arXiv.1006.0868>

Murray, I., Adams, R., & MacKay, D. (2010, March). Elliptical slice
sampling. In *Proceedings of the thirteenth international conference on
artificial intelligence and statistics* (pp. 541-548). JMLR Workshop and
Conference Proceedings. <https://doi.org/10.48550/arxiv.1001.0175>

Nishihara, R., Murray, I., & Adams, R. P. (2014). Parallel MCMC with
generalized elliptical slice sampling. *The Journal of Machine Learning
Research, 15(1)*, 2087-2112. <https://doi.org/10.48550/arXiv.1210.7477>
