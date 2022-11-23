
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- # GPHS -->
<!-- badges: start -->
<!-- badges: end -->

## GPHS: Hyperparameter selection in Gaussian Process models

The R package ‘GPHS’ provides functions for MCMC update of the
hyperparameters in covariance function in Gaussian process models. It
includes the Elliptical Slice Sampling algorithm as discussed in
[Nishihara, Murray and Adams
(2014)](https://www.cs.princeton.edu/~rpa/pubs/nishihara2014generalized.pdf)
and the Surrogate data model Metropolis-Hastings algorithm as discussed
in [Murray et. al (2010)](https://arxiv.org/abs/1001.0175). It also
includes a function to provide the user with a comparison of the
procedures using graphical tools.

### Installation

You can install the most recent version of ‘GPHS’ package from
[GitHub](https://github.com/niladrik/GPHS) using the following commands:

``` r
# Plain installation
devtools::install_github("niladrik/gphs") # gphs package
# For installation with vignette
devtools::install_github("niladrik/gphs", build_vignettes = TRUE)
```

### Remainder of work
