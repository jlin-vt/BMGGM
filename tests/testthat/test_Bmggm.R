library(mvtnorm)
library(glasso)
library(Matrix)
library(JGL)
library(pscl)
library(lattice)
library(MCMCpack)
library(MASS)
library(ROCR)
library(fifer)
library(edgebundleR)
library(igraph)
library(data.table)
library(huge)

# Data
p <- 10
K <- 4
n <- 400
dat <- GenerateData(p, K, n)

# set the options for MCMC
options <- list()
options$burnin <- 100
options$nmc <- 100

# intialize the priors
PriorPar <- list()
PriorPar$a <- 1
PriorPar$b <- 5
PriorPar$a0 <- 1
PriorPar$b0 <- 10
PriorPar$eps <- 10000
PriorPar$delta <- 1
PriorPar$c <- 100
PriorPar$Theta <- matrix(0.2, K, K)

# initial the updates
InitVal <- list()
InitVal$sigma2 <- 1
InitVal$mu <- rep(0, p * K)
InitVal$Beta <- matrix(runif((p * K) * (p * K)), p * K, p * K)
InitVal$adj <- ifelse(InitVal$Beta, 1, 0)

# tic
ptm <- proc.time()

# Run
res <- Bmggm(dat, options, PriorPar, InitVal)
adj_save <- res$adj_save

# tok
Time.used <- proc.time() - ptm
cat("Time used (min) = ", Time.used[3]/60, "\n")