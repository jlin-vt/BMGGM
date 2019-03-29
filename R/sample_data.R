#' Generate Data by Gaussian Graphical Model
#' 
#' @param p dimension of each pathway.
#' @param K number of pathways.
#' @param n number of observations.
#' @param network_type choice of type of network: 'ar(2)', 'chain', 'random' or 'scale-free'.
#' 
#' @return a list of objets:
#'  data data.
#'  z_P pathway membership of genes.
#'  A Precsion Matrix for the whole network. 
#'  
#' @export
GenerateData <- function(p, K, n, network_type = "ar2") {
  # initialize
  z_P <- c()
  Omegas <- list()
  A <- matrix(0, p * K, p * K)
  
  # set network for each block
  for (k in 1:K) {
    z_P <- c(z_P, rep(k, p))
    if (network_type == "ar2") 
      Omega <- toeplitz(c(1, 0.5, 0.25, rep(0, p - 3)))
    if (network_type == "chain") {
      Omega <- make_ring(p)
      Omega <- as.matrix(get.adjacency(Omega))
      L <- Omega
      L[upper.tri(L)] <- 0
      L <- ifelse(L == 1, sample(c(runif(1, -0.4, -0.1), runif(1, 
                                                               0.1, 0.4)), 1), 0)
      Omega <- L + t(L) + diag(p)
    }
    if (network_type == "random") {
      Omega <- sample_gnm(n = p, m = p)
      Omega <- as.matrix(get.adjacency(Omega))
      L <- Omega
      L[upper.tri(L)] <- 0
      L <- ifelse(L == 1, sample(c(runif(1, -0.4, -0.1), runif(1, 
                                                               0.1, 0.4)), 1), 0)
      Omega <- L + t(L) + diag(p)
    }
    if (network_type == "scale-free") {
      Omega <- sample_pa(n = p, directed = F)
      Omega <- as.matrix(get.adjacency(Omega))
      L <- Omega
      L[upper.tri(L)] <- 0
      L <- ifelse(L == 1, runif(1, 0.5, 0.99), 0)
      Omega <- L + t(L) + diag(p)
    }
    A[which(z_P == k), which(z_P == k)] <- Omega
    Omegas[[k]] <- Omega
  }
  
  # make sure precision matrix must be positive definite
  if (network_type == "scale-free") 
    A <- FixMatrix(A, 1.5)
  if (any(eigen(A)$values < 0)) 
    stop("Precision matrix must be positive definite")
  
  # generdate data based on precision matrix
  Cov_True <- solve(A)
  data <- mvrnorm(n, rep(0, p * K), Cov_True)
  return(list(data = data, z_P = z_P, A = A))
}