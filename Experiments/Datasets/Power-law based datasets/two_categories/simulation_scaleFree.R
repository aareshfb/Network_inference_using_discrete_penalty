# data simulation for testing network inference 

# start with a hypergraph connecting communities,
# simulate data in this order : only shifts in edges (and/or) shifts 
# in connectivities

## simulate networks with hypergraph, with different categories in each population

# generate samples from simulated covariances. compute distance between 
# sample covariance matrices to inform network inference

library(igraph)
library(Matrix)
library(matrixcalc)
library(igraph)
library(glue)
library(ape)
# igraph_options(add.vertex.names=F)
# igraph_options()$add.vertex.names


setwd("/Volumes/umms-ukarvind/shared_data/l0-networks/simulations-category/")

dir.create('simulation_data')
# setwd(outdir)

#------------------------
#set structure of graphs
#------------------------
p = 250 #no of features
M = 5 #no of disconnected subclasses

#------------------------
#------------------------
## perturb networks
perturbNets <- function(mm, cpt, ss){
  ## make perturbations symmetrix
  p <- nrow(mm)
  idx <- sample(ss, size = cpt*p*p)
  mm[idx] <- runif(length(idx), min = -2, max = 2)
  mm <- (mm + t(mm)) / 2
  ## make diagonals larger than off diag elements
  diag(mm) <- diag(mm) + rowSums(abs(mm))
  mm
}

#------------------------
#------------------------

#------------------------
#loop over hypergraph size
#------------------------
kk.list <- c(3, 5, 10, 20)

for(kk in kk.list){
  message(glue('hypergraph size : ', kk))
  
  ## make output directory for files
  outdir <- glue('./simulation_data/hypergraph_{kk}')
  if(!exists(outdir)){
    dir.create(outdir)
    dir.create(glue('{outdir}/true_networks'))
    dir.create(glue('{outdir}/covariance_matrix'))
  }
  
  #simulate base network with a modular structure
  cov_list <- list()
  omega_list <- list()
  
  for(j in 1:M){
    message('simulating module ', j)
    l <- simulate_powerlaw_networks(p = p/M)
    cov_list[[j]] <- l[["cov"]] ; omega_list[[j]] <- l[["A"]]
  }
  
  sigma <- bdiag(cov_list); image(sigma)
  image(omega_list[[j]])
  
  #-----------------------------
  #derive MST through hypergraph
  #-----------------------------
  hg <- matrix(rnorm(n = kk*kk, mean = 10, sd = 3), nrow = kk)
  hg <- (hg + t(hg))/2 #make symmetric distance 
  rownames(hg) <- paste("cluster", 1:kk, sep = ""); colnames(hg) <- paste("cluster", 1:kk, sep = "")
  m <- ape::mst(hg) #learn mst through hyper-graph
  
  #given the mst, use bfs to learn path through graph
  g_mst <- graph_from_adjacency_matrix(m, mode = "undirected")
  plot(g_mst)
  
  ## save adjacency matrix to file
  write.table(as.matrix(m), glue('{outdir}/MST.txt'), quote = F, row.names = F, col.names = F)
  
  pdf(glue('{outdir}/simulation-mst.pdf'), width = 8, height = 8)
  par(mar=c(.5,.5,.5,.5))
  plot(g_mst, vertex.size=3, vertex.label.font=2, 
       vertex.label.cex=.75, layout = layout_with_fr)
  dev.off()
  
  ## find a leaf node to set as root
  degree.mst <- igraph::degree(g_mst)
  plot(degree.mst, pch=20)
  
  message(d.max <- max(degree.mst))
  leaf.nodes <- which(degree.mst == 1)
  # leaf.nodes <- names(degree.mst[degree.mst == 1])
  root <- sample(leaf.nodes, 1)
  message(glue('Root node is {root}'))
  
  mbfs <- bfs(g_mst, root = root, father = T, order = T) 
  bfs_order <- as.vector(mbfs$order)
  bfs_pred <- as.vector(mbfs$father)
  
  #navigate this graph and simulate network for each cluster
  network_list <- vector(mode = "list", length = kk)
  root_node <- bfs_order[1]
  network_list[[root_node]] <- omega_list 
  
  
  #-----------------------------
  #generate sub-population networks
  #-----------------------------
  s = 3
  delta <- 1 #magnitude of perturbation
  
  #use bfs to traverse the hyper-graph properly
  ## generate common graphs, perturb to obtain category-specific networks
  for(node in 2:kk){
    
    j <- bfs_order[node]; i <- bfs_pred[j]; 
    message(i, " ", j)  #, " ", wt)
    net_j <- network_list[[i]] #initialize with predecessor
    
    ind <- sample(1:M, s, replace = F) #randomly sample module to perturb
    
    ### comment out this portion if we want to keep structure same
    ## if parent node is a branching point, simulate new module
    if(degree.mst[i] >= max(2, d.max - 1)){
      message('simulating new module')
      net_j[[ind[1]]] <- simulate_powerlaw_networks(p = p/M, alpha = 1.2)[["A"]]
    }
    
    for(m in 2:s){
      message('module ', ind[m])
      net_j[[ind[m]]] <- perturb_net_weights(net_j[[ind[m]]], delta)  
    }
    network_list[[j]] <- net_j
  }
  
  ## perturb networks for each category
  ## change edge structure compared to global network
  cpt <- .025  ## change 2.5% edges for each category
  ss = setdiff(1:(p*p), (1:p)^2) ## leave out diagonal elements
  c1net <- lapply(network_list, function(net) lapply(net, perturbNets, cpt, ss))
  c2net <- lapply(network_list, function(net) lapply(net, perturbNets, cpt, ss))
  
  c1prec <- lapply(c1net, bdiag)
  c2prec <- lapply(c2net, bdiag)
  
  ## write networks to file
  precision_list <- lapply(network_list, function(omega) bdiag(omega))
  
  for(j in seq_along(c1prec)){
    write.csv(as.matrix(precision_list[[j]]), 
              file.path(outdir, 
                        glue('true_networks/global_clust_precision_', j, '.csv')), 
              quote = F, row.names = F)
    
    write.csv(as.matrix(c1prec[[j]]), 
              file.path(outdir, 
                        glue('true_networks/cat1_clust_precision_', j, '.csv')), 
              quote = F, row.names = F)
    
    write.csv(as.matrix(c2prec[[j]]), 
              file.path(outdir, 
                        glue('true_networks/cat2_clust_precision_', j, '.csv')), 
              quote = F, row.names = F)
  }
  
  #-----------------------------
  #convert modules into covariances
  #-----------------------------
  c1cov <- lapply(c1net, function(net) lapply(net, convert_covariance))
  c1cov <- lapply(c1cov, bdiag)
  
  c2cov <- lapply(c2net, function(net) lapply(net, convert_covariance))
  c2cov <- lapply(c2cov, bdiag)
  
  ## save covariance matrices for generating data
  for(j in seq_along(c1cov)){
    message(j)
    write.csv(as.matrix(c1cov[[j]]), 
              file.path(outdir, 
                        glue('covariance_matrix/cat1_clust_covariance_', j, '.csv')), 
              quote = F, row.names = F)
    
    write.csv(as.matrix(c2cov[[j]]), 
              file.path(outdir, 
                        glue('covariance_matrix/cat2_clust_covariance_', j, '.csv')), 
              quote = F, row.names = F)
  }
}





# figs <- lapply(precision_list, function(x)  image(x))
# for(j in 1:length(precision_list)){
#   
#   pdf(glue('{outdir}/cluster_',j,'_precision.pdf'), width = 10)
#   plot(figs[[j]])
#   dev.off()
#   write.csv(as.matrix(precision_list[[j]]), 
#             file.path(outdir, 
#                       glue('clust_precision_', j, '.csv')), 
#             quote = F, row.names = F)
# }

# write_rds(precision_list, "simulation_data/precisionList.rds")


#-------------------------------
#simulate data for each cluster
#-------------------------------
# generateData <- function(cov_list, category){
#   p <- nrow(cov_list[[1]])
#   n <- 20*p #number of samples
#   k <- length(cov_list)
#   data.clusters <- vector(mode = "list", length = k)
#   
#   for(j in 1:k){
#     message('simulating cluster ', j)
#     #simulate normal variables with this covariance distribution
#     cov_matrix <- cov_list[[j]]
#     U <- chol(cov_matrix)
#     z <- matrix(rnorm(n = n*p), nrow = n)
#     data.clusters[[j]] <- z%*%U
#   }
#   
#   if(!exists(category)){
#     dir.create(category)
#   }
#   
#   for(j in 1:k){
#     write.csv(as.matrix(data.clusters[[j]]), 
#               file.path(outdir, category,
#                         glue('simulation_clust_', j, '.csv')), quote = F, row.names = F)
#   }
# }
# 
# generateData(c1cov, 'cat1')
# generateData(c2cov, 'cat2')

#-----------------------------------------------------
## save graph structures
#-----------------------------------------------------
# outdir <- 'graph-structure'
# dir.create(outdir)
# for(i in seq_along(c1net)){
#   write.csv(as.matrix(precision_list[[i]]), 
#             file.path(outdir,
#                       glue('global_precision_', i, '.csv')),
#             quote = F, row.names = F)
#   
#   write.csv(as.matrix(c1prec[[i]]), 
#             file.path(outdir,
#                       glue('cat1_precision_', i, '.csv')),
#             quote = F, row.names = F)
#   
#   write.csv(as.matrix(c2prec[[i]]), 
#             file.path(outdir,
#                       glue('cat2_precision_', i, '.csv')),
#             quote = F, row.names = F)
#   
# }

#-----------------------------------------------------
#plot graph structures
#-----------------------------------------------------

# c1graphs <- generateGraphs(c1prec)
# c2graphs <- generateGraphs(c2prec)
# 
# g1 <- c1graphs[[1]]; g2 <- c2graphs[[1]]
# g <- g1 - g2
# g <- delete.vertices(g, igraph::degree(g) == 0)
# par(mfrow=c(1,1))
# plot(g, layout = layout_with_fr, vertex.size = .5)
