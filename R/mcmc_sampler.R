#' Bayesian multiple Gaussian graphical models by MCMC.
#' 
#' @param dat a list of objets:
#'  n: number of observations.
#'  p: dimension of each pathway.
#'  K: number of pathways.
#'  z_P: indicator vector of genes membership
#'  P: dimension of the data.
#'  
#' @param options a list of objets:
#'  burnin: number of MCMC iterations before burnin.
#'  nmc: number of MCMC iterations after burnin.
#'  
#' @param PriorPar a list of objets:
#'  a: shape1 parameter for Theta for off-digonal block.
#'  b: shape2 parameter for Theta for off-digonal block.
#'  a0: shape1 parameter for Theta for digonal block.
#'  b0: shape2 parameter for Theta for digonal block.
#'  eps: rate parameter for v0^2.
#'  delta: shape parameter for v0^2.
#'  c: the parameter for decision boundary of spike-and-slab.
#'  Theta: a K x K initial graph PPI matrix.
#'  
#' @param InitVal a list of objets:
#'  mu: intercept term.
#'  sigma2: overall noise level, same across groups.
#'  Beta: a P x P initial coefficient matrix. 
#'  adj: a P x P initial adjacency matrix. 
#'  
#' @return a list of objets:     
#'  Beta_save: p x p x K x nmc sample of coefficient matrix.
#'  adj_save: p x p x K x nmc sample of adjacency matrix.
#'  Theta_save: K x K x nmc sample of graph similarity matrix.
#'  
#' @export
Bmggm <- function(dat, options, PriorPar, InitVal) {
  # load data
  z_P <- dat$z_P
  data <- dat$data
  n <- dim(data)[1]
  P <- dim(data)[2]
  p <- P/K
  
  # load options
  burnin <- options$burnin
  nmc <- options$nmc
  
  # load priors
  a <- PriorPar$a
  b <- PriorPar$b
  a0 <- PriorPar$a0
  b0 <- PriorPar$b0
  eps <- PriorPar$eps
  delta <- PriorPar$delta
  c <- PriorPar$c
  Theta <- PriorPar$Theta
  
  # load initials
  sigma2 <- InitVal$sigma2
  mu <- InitVal$mu
  Beta <- InitVal$Beta
  adj <- InitVal$adj
  
  # set up matrices for return values
  Beta_save <- array(NA, c(P, P, nmc))
  adj_save <- Beta_save
  Theta_save <- array(0, c(K, K, nmc))

  # MCMC sampling
  for (iter in 1:(burnin + nmc)) {
    # show iterations
    if (iter%%100 == 0) {
      cat("iter =", iter, "\n")
    }
    
    for (i in 1:(P - 1)) {
      col_start_index <- i + 1
      col_end_index <- P
      
      # update Beta 
      Y <- as.matrix(data[, i])
      X <- as.matrix(data[, col_start_index:col_end_index] - mu[col_start_index:col_end_index])
      X <- cbind(1, X)
      v0 <- rigamma(ncol(X), eps - 1, median((Beta[i, col_start_index:col_end_index]))/2 + delta)
      v0 <- sqrt(v0)
      v1 <- c * v0
      gamma <- c(1, adj[i, col_start_index:col_end_index])
      invSigG <- diag(gamma * (1/v1^2) + (1 - gamma) * (1/v0^2),  ncol = ncol(X))
      KK <- t(X) %*% X/sigma2 + invSigG
      invK <- solve(KK)
      MM <- (t(X) %*% Y)/sigma2
      Beta[i, col_start_index:col_end_index] <- rmvnorm(n = 1, mean = invK %*% MM, sigma = invK)[-1]
      Beta[col_start_index:col_end_index, i] <- Beta[i, col_start_index:col_end_index]
      mu[i] <- rmvnorm(n = 1, mean = invK %*% MM, sigma = invK)[1]
      
      # update G 
      prior <- Theta[z_P[i], z_P[col_start_index:col_end_index]]
      v1 <- v1[-1]
      v0 <- v0[-1]
      log_odds_adj <- -log(v1) - Beta[i, col_start_index:col_end_index]^2/(2 *  v1^2) + log(prior) + log(v0) + Beta[i, col_start_index:col_end_index]^2/(2 * v0^2) - log(1 - prior)
      adj[i, col_start_index:col_end_index] <- ifelse(log_odds_adj > 100, 1, rbinom(ncol(X), size = 1, prob = exp(log_odds_adj)/(exp(log_odds_adj) + 1)))
      adj[col_start_index:col_end_index, i] <- adj[i, col_start_index:col_end_index]
    }  
    
    # compute means
    mu[P] <- mean(data[, P])
    
    # update Theta 
    for (m in 1:K) {
      for (k in m:K) {
        # Get terms that are a sum over all edges on log scale
        sum_over_edges <- 0
        sum_over_non_edges <- 0
        for (i in 1:(P - 1)) {
          for (j in (i + 1):P) {
            if (z_P[i] == m && z_P[j] == k) {
              sum_over_edges <- sum_over_edges + adj[i, j]
              sum_over_non_edges <- sum_over_non_edges + 1 -  adj[i, j]
            }  
          }  
        }  
        Theta[m, k] <- ifelse(m == k, rbeta(1, sum_over_edges + a0, sum_over_non_edges + b0), rbeta(1, sum_over_edges + a, sum_over_non_edges + b))
        Theta[k, m] <- Theta[m, k]
      } 
    }  
    
    # Retain values for posterior sample 
    if (iter > burnin) {
      Beta_save[, , iter - burnin] <- Beta[, ]
      adj_save[, , iter - burnin] <- adj[, ]
      Theta_save[, , iter - burnin] <- Theta
    }
  }  
  return(list(Beta_save = Beta_save, adj_save = adj_save))
}