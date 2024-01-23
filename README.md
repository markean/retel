
<!-- README.md is generated from README.Rmd. Please edit that file -->

# retel

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/markean/retel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/markean/retel/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

retel implements the regularized exponentially tilted empirical
likelihood method. The proposed method removes the convex hull
constraint using a novel regularization technique, providing a suitable
pseudo-likelihood for Bayesian inference.

The following functions enable users to set estimating functions by
providing data and parameters:

- `etel()` computes exponentially tilted empirical likelihood without
  regularization.
- `retel()` computes regularized exponentially tilted empirical
  likelihood with regularization parameters.

This repository accompanies the research paper titled ‘Regularized
Exponentially Tilted Empirical Likelihood for Bayesian Inference,’
available on [arXiv](https://arxiv.org/abs/2312.17015). The
`retel-paper/` folder contains code and outputs from the paper. This
work was supported by the U.S. National Science Foundation under Grants
No. [SES-1921523](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1921523)
and
[DMS-2015552](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2015552).

## Installation

You can install the latest stable release of retel from CRAN.

``` r
install.packages("retel")
```

### Development version

You can install the development version of retel from GitHub.

``` r
# install.packages("pak")
pak::pak("markean/retel")
```

## Usage

``` r
library(retel)

# Generate data
set.seed(63456)
x <- rnorm(100)

# Define an estimating function (ex. mean)
fn <- function(x, par) {
  x - par
}

# Set parameter value
par <- 0

# Set regularization parameters
mu <- 0
Sigma <- 1
tau <- 1

# Call the retel function to compute the log-likelihood ratio. The return value
# contains the optimization results as the attribute 'optim'.
retel(fn, x, par, mu, Sigma, tau)
#> [1] -0.06709306
#> attr(,"optim")
#> 
#> Call:
#> 
#> nloptr(x0 = rep(0, p), eval_f = eval_obj_fn, eval_grad_f = eval_gr_obj_fn, 
#>     opts = opts, g = g, mu = mu, Sigma = Sigma, n = n, tau = tau)
#> 
#> 
#> Minimization using NLopt version 2.7.1 
#> 
#> NLopt solver status: 1 ( NLOPT_SUCCESS: Generic success return value. )
#> 
#> Number of Iterations....: 4 
#> Termination conditions:  xtol_rel: 1e-04 
#> Number of inequality constraints:  0 
#> Number of equality constraints:    0 
#> Optimal value of objective function:  0.999330716387232 
#> Optimal value of controls: -0.03738174
```
