# BMGGM
[![Build Status](https://travis-ci.org/jlin-vt/BMGGM.svg?branch=master)](https://travis-ci.org/jlin-vt/BMGGM)

Learning multiple Gaussian graphical models by Bayesian

## Installation
This package can be installed using the `devtools` package in R:

```r
library(devtools)
devtools::install_github("jlin-vt/BMGGM")
library(BMGGM)
```

## A simple example

To get started, the user is recommended to generate some synthetic data.
```r 
set.seed(50)
p <- 10
K <- 4
n <- 400
dat <- GenerateData(p, K, n)
```

The second step is to set the options for MCMC.
```r 
options <- list()
options$burnin <- 10000
options$nmc <- 10000
```

You also need to intialize the priors.
```r 
PriorPar <- list()
PriorPar$a <- 1
PriorPar$b <- 5
PriorPar$a0 <- 1
PriorPar$b0 <- 10
PriorPar$eps <- 10000
PriorPar$delta <- 1
PriorPar$c <- 100
PriorPar$Theta <- matrix(0.2, K, K)
```

Intialize the updates for the parameters.
```r 
InitVal <- list()
InitVal$sigma2 <- 1
InitVal$mu <- rep(0, p * K)
InitVal$Beta <- matrix(runif((p * K) * (p * K)), p * K, p * K)
InitVal$adj <- ifelse(InitVal$Beta, 1, 0)
```

Finally, apply MCMC sampler to execute BMGGM:
```r 
# Run
res <- Bmggm(dat, options, PriorPar, InitVal)
adj_save <- res$adj_save
```

The [vignette](https://github.com/jlin-vt/BMGGM/blob/master/vignettes/examples.Rmd) demonstrates example usage of all main functions. 

## Status
The preprint describing the corncob methodology is available [here](https://github.com/jlin-vt/BMGGM/blob/master/vignettes/examples.Rmd). The manuscript has been submitted to _Biometrics_.

## Bug Reports / Change Requests
If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/jlin-vt/BMGGM/issues).

## License
The package is available under the terms of the GNU General Public License v3.0.
