
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
and Metropolis-Hastings algorithm and Slice sampling algorithm for the
Surrogate data model as discussed in [Murray et. al
(2010)](https://arxiv.org/abs/1001.0175). It also includes a function to
provide the user with a comparison of the procedures using graphical
tools.

### Installation

You can install the most recent version of ‘GPHS’ package from
[GitHub](https://github.com/niladrik/GPHS) using the following commands:

``` r
# Plain installation
devtools::install_github("niladrik/GPHS") # GPHS package
# For installation with vignette
devtools::install_github("niladrik/GPHS", build_vignettes = TRUE)
```

### Remainder of work

In the remainder of the semester, I plan to complete the function for
surrogate data model and that for graphical comparison of the
algorithms. I will include compatibility checks for each functions and
then test for the correctness. Also, I will have to write a vignette and
include examples for the functions.
