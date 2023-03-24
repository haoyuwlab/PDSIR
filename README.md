
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PDSIR

<!-- badges: start -->
<!-- badges: end -->

The PDSIR R package implements an efficient data augmentation MCMC
(DA-MCMC) algorithm for fitting the general stochastic epidemic model to
incidence data. This package contains the code used in the paper
“Uniformly Ergodic Data-Augmented MCMC for Fitting the General
Stochastic Epidemic Model to Incidence Data” by R. Morsomme and J. Xu
available on ArXiv. The novelty of our DA-MCMC algorithm is the *joint*
update of the high-dimensional latent data in a Metropolis-Hastings
step. Compared to existing single-site sampler, our block sampler for
the latent data significantly improves the mixing of the Markov chain.

## Installation

You can install the development version of PDSIR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rmorsomme/PDSIR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PDSIR)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
