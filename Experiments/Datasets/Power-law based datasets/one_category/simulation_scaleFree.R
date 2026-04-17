# data simulation for testing network inference 

# start with a hypergraph connecting communities,
# simulate data in this order : only shifts in edges (and/or) shifts 
# in connectivities

# generate samples from simulated covariances. compute distance between 
# sample covariance matrices to inform network inference

library(igraph)
library(Matrix)
library(matrixcalc)
library(igraph)
library(glue)
library(ape)

setwd("/Volumes/umms-ukarvind/shared_data/l0-networks/")

if(!exists('simulation_data')){
  dir.create('simulation_data')
}


#------------------------
#set structure of graphs
#------------------------
p = 250 #no of features
M = 5 #no of disconnected subclasses

#------------------------
#loop over hypergraph size
#------------------------
kk.list <- c(3, 5, 10, 20, 50, 100)

for(kk in kk.list){
  message(glue('hypergraph size : ', kk))
  
  ## make output directory for files
  outdir <- glue('./simulation_data/hypergraph_{kk}')
  dir.create(outdir)
  dir.create(glue('{outdir}/true_networks'))
  dir.create(glue('{outdir}/covariance_matrix'))
  
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
  pdf('hypergraph.pdf', width=6,height=6)
  plot(g_mst, vertex.label.cex= 0.25,
                # vertex.label.color = ifelse(V(g)$name %in% tfs, 'darkmagenta', 'gray'),
                vertex.size = 5, edge.width=2,
                layout = layout_with_fr, main = 'Population Hypergraph')
  dev.off()              
  
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
  ## generate network structure
  for(node in 2:kk){
    j <- bfs_order[node]; i <- bfs_pred[j]; 
    message(i, " ", j)  #, " ", wt)
    net_j <- network_list[[i]] #initialize with predecessor
    
    ind <- sample(1:M, s, replace = F) #randomly sample module to perturb
    
    ## comment out this portion if we want to keep structure same
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
  
  precision_list <- lapply(network_list, function(omega) bdiag(omega))
  
  for(j in 1:length(precision_list)){
    
    write.csv(as.matrix(precision_list[[j]]), 
              file.path(outdir, 
                        glue('true_networks/clust_precision_', j, '.csv')), 
              quote = F, row.names = F)
  }
  
  #-----------------------------
  #convert modules into covariances
  #-----------------------------
  cov_list <- lapply(network_list, function(module_list)
    lapply(module_list, function(module) convert_covariance(module)))
  
  cov_list <- lapply(cov_list, function(x) bdiag(x))
  
  ## save covariance matrices for generating data
  for(j in 1:length(cov_list)){
    message(j)
    write.csv(as.matrix(cov_list[[j]]), 
              file.path(outdir, 
                        glue('covariance_matrix/clust_covariance_', j, '.csv')), 
              quote = F, row.names = F)
  }
  
  
}


