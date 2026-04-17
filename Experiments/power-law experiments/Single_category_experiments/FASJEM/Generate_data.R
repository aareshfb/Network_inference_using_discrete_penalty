Generate_Data <- function(L, n, p, rdseed) {
  set.seed(rdseed)
  fname <- "/nfs/turbo/umms-ukarvind/shared_data/l0-networks/simulations-hypergraphs/simulation_data/hypergraph_"
  
  Data <- list()
  
  for (i in 1:L) {
    file_path1 <- paste0(fname, L, "/covariance_matrix/clust_covariance_", i, ".csv")
    cov <- as.matrix(read.csv(file_path1))
    cov <- chol(cov)
    
    A <- matrix(rnorm(n * p), n, p)
    Data[[i]] <- A %*% cov
    
    cat("Population", i, "completed\r")
  }
  
  return(Data)
}

Read_Simulation_Data <- function(L, file_ext = ".csv", sep = ",") {
  # Read expression counts data and return list
  Data <- list()
  fpath <- "/nfs/turbo/umms-ukarvind/shared_data/l0-networks/simulations-hypergraphs/simulation_data/hypergraph_"
  fname = '/covariance_matrix/clust_covariance_'
  for (i in 1:L) {
    file_path <- paste0(fpath,L,fname ,i, file_ext)
    Data[[i]] <- as.matrix(read.table(file_path, sep = sep, header = TRUE))
  }
  
  return(Data)
}