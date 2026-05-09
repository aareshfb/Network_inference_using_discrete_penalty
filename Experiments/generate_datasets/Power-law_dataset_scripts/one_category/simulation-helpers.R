#simulating structure for each module
#module structure is sampled using igraph, with the desired power-parameter
#each module has n_m genes (assume equally sized sub-networks), and follow a given degree distribution
#hub-networks and power-law networks

nonzero <- function(x){
  ## function to get a two-column matrix containing the indices of the
  ### non-zero elements in a "dgCMatrix" class matrix
  
  stopifnot(inherits(x, "dgCMatrix"))
  if (all(x@p == 0))
    return(matrix(0, nrow=0, ncol=2,
                  dimnames=list(character(0), c("row","col"))))
  res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
  colnames(res) <- c("row", "col")
  res <- res[x@x != 0, , drop = FALSE]
  return(res)
}

# write a function to plot the degree distribution
plot_degree_distribution = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
       col = 1, main = "Degree Distribution")
}

#function to simulate network edge weights
simulate_edges <- function(n_edge){
  x <- runif(n = n_edge, min = 0.4, max = 1) #simulate edge weights
  y <- sample(0:1, size = n_edge, replace = T) #simulate pos vs neg interaction
  ind <- which(y > 0)
  x[ind] <- -x[ind]
  return(x)
}

simulate_powerlaw_networks <- function(p = 100, alpha = 1.2){
  ## simulate power-law networks with preferential attachment alpha 
  g <- sample_pa(n = p, power = alpha, m = 1)
  # plot(g)
  
  #get adjacency matrix for graph
  mat <- as_adjacency_matrix(g)
  ind <- nonzero(mat)
  
  edge_wts <- simulate_edges(dim(ind)[1])
  
  A <- Matrix(0, nrow = p, ncol = p)
  A[ind] <- edge_wts
  A <- t(A)
  
  #make matrix diagonally dominant 
  scale_factor <- rowSums(abs(A)) #sum absolute values of off-diagonal entries
  diag(A) <- scale_factor + 1
  
  #making matrix symmetric
  A <- (A + t(A)) / 2
  image(A)
  # message("check for Positive definitenss: ", is.positive.semi.definite(as.matrix(A)))
  
  #compute covariance matrix
  cov_matrix <- convert_covariance(A)
  # image(cov_matrix)
  
  return(list("cov" = cov_matrix, "A" = A))

}


perturb_net_weights <- function(module, delta = 1){
  
  ind_nz <- nonzero(module); 
  n <- dim(ind_nz)[1]
  
  #perturb non-zero edge weights in network
  module[ind_nz] <- module[ind_nz] + runif(n, min = -delta, max = delta)
  
  #post-processing to make matrix positive definite
  #make matrix positive definite
  diag(module) <- 0
  scale_factor <- rowSums(abs(module))  #sum absolute values of off-diagonal entries
  diag(module) <- 1 + scale_factor
  
  #making matrix symmetric 
  module <- triu(module) + t(triu(module))
  image(module)
  # message("check for Positive definitenss: ", is.positive.semi.definite(as.matrix(module)))
  return(module)
  
}

#convert the concentration matrices to covariances
convert_covariance <- function(A){
  B <- solve(A)
  d <- diag(B) #diagonal entries are variance terms, should be positive
  C <- sqrt(d%*%t(d)); #scaling matrix for covariance
  cov_matrix <- B / C
  cov_matrix <- drop0(cov_matrix, tol = 1e-4)
  return(cov_matrix)
}


