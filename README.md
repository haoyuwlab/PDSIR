
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PDSIR

<!-- badges: start -->
<!-- badges: end -->

The `PDSIR` R package implements an efficient data augmentation MCMC
(DA-MCMC) algorithm for exact Bayesian inference under the semi-Markov
stochastic susceptible-infectious-removed (SIR) model, given discretely
observed counts of infections. The novelty of this DA-MCMC algorithm is
the *joint* update of the high-dimensional latent data. In a
Metropolis-Hastings step, the latent data are jointly proposed from a
surrogate process carefully designed to closely resemble the target
process and from which we can efficiently generate epidemics consistent
with the observed data. This yields a MCMC algorithm that explores the
high-dimensional latent space efficiently, mixes significantly better
than single-site samplers, and scales to outbreaks with thousands of
infections.

## Installation

You can install the development version of PDSIR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rmorsomme/PDSIR")
```

## Proof of concept

We employ the DA-MCMC on artifical data; semi-Marko model

This package contains the code used in the paper “Uniformly Ergodic
Data-Augmented MCMC for Fitting the General Stochastic Epidemic Model to
Incidence Data” by R. Morsomme and J. Xu available on ArXiv. We use it
to fit a semi-Markov susceptible-infectious-removed model to the
2013-2015 outbreak of Ebola Haemorrhagic Fever in Gu'eck'edou, Guinea.

This is a basic example which shows you how to solve a common problem:

``` r
library(PDSIR)
## basic example code

# setup
S0 <- 500 # initial number of susceptible individuals
I0 <- 5
t_end <- 6

iota_dist <- "weibull"
theta <- list(R0 = 2, lambda = 1, shape = 2) # reproduction rate R_0; parameters of the Weibull distributed infection periods

theta <- complete_theta(theta, iota_dist, S0)

# Simulate artificial data
SIR <- simulate_SEM(S0, I0, t_end, theta, iota_dist)

  # Observed data
#draw_trajectories(SIR, plot_id, path, t_end)
#  Y     <- observed_data(SIR, K)
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
