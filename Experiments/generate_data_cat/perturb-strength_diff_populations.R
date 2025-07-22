## vary level of perturbation in categorical networks
## this is the live version! keep SAFE!
library(igraph)
library(Matrix)
library(matrixcalc)
library(igraph)
library(glue)
library(ape)

#------------------------
#set structure of graphs
#------------------------
p = 250 #no of features
M = 5 #no of disconnected sub-classes
psx = 1.0 #modify this value in the for loop below as well #psx: perturb strength  

#setwd("/Volumes/umms-ukarvind/shared_data/l0-networks/simulations-category/")
setwd("Q:/shared_data/l0-networks/simulations-category/simulation-codes")
if (!dir.exists('../simulation_data/vary_T')){
  dir.create('../simulation_data/vary_T')
}

outdir <- '../simulation_data/vary_T/'

#------------------------
#------------------------
## simulate global network with fixed hypergraph size
kk <- 100 ## number of populations
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
dir.create(glue(outdir, 'hypergraph_{kk}_P{psx}'))
write.table(as.matrix(m), glue('{outdir}hypergraph_{kk}_P{psx}/MST.txt'), quote = F, row.names = F, col.names = F)

pdf(glue('{outdir}hypergraph_{kk}_P{psx}/simulation-mst.pdf'), width = 8, height = 8)
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


## loop through perturbation strengths and generate datasets
#------------------------
## perturb a fraction of these edges
## keep level of sparsity same for all networks, so remove same number of edges from global network
#------------------------
## perturb networks
perturbNets <- function(net, pstr){
  ## make perturbations symmetric
  p <- nrow(net)
  nnz <- (length(net@i) -p) / 2 ## number of non-zero edges 
  idx <- sample(ss, size = pstr*nnz)
  net[idx] <- runif(length(idx), min = -2, max = 2)
  net <- (net + t(net)) / 2
  ## make matrix diagnonally dominant
  diag(net) <- diag(net) + rowSums(abs(net))
  net
}

#------------------------
#------------------------
## join modules and introduce perturbations in network as a whole
## perturbed networks no longer block-diagonal
precision_list <- lapply(network_list, function(omega) bdiag(omega))
p <- nrow(precision_list[[1]])
ss = setdiff(1:(p*p), (1:p)^2) ## leave out diagonal elements

## loop through perturbation strength
for(psx in c(1.0)){
  
  message(psx)
  #dir.create(glue(outdir, 'hypergraph_{kk}_P{psx}'))
  dir.create(glue(outdir, 'hypergraph_{kk}_P{psx}/true_networks'))
  dir.create(glue(outdir, 'hypergraph_{kk}_P{psx}/covariance_matrix'))
  
  c1net <- lapply(precision_list, perturbNets, pstr = psx)
  c2net <- lapply(precision_list, perturbNets, pstr = psx)
  
  for(j in seq_along(c1net)){
    write.csv(as.matrix(precision_list[[j]]), 
              file.path(outdir, 
                        glue('hypergraph_{kk}_P{psx}/true_networks/global_clust_precision_', j, '.csv')), 
              quote = F, row.names = F)
    
    write.csv(as.matrix(c1net[[j]]), 
              file.path(outdir, 
                        glue('hypergraph_{kk}_P{psx}/true_networks/cat1_clust_precision_', j, '.csv')), 
              quote = F, row.names = F)
    
    write.csv(as.matrix(c2net[[j]]), 
              file.path(outdir, 
                        glue('hypergraph_{kk}_P{psx}/true_networks/cat2_clust_precision_', j, '.csv')), 
              quote = F, row.names = F)
  }
  
  #-----------------------------
  #convert modules into covariances
  #-----------------------------
  c1cov <- lapply(c1net, convert_covariance)
  c2cov <- lapply(c2net, convert_covariance)

  ## save covariance matrices for generating data
  for(j in seq_along(c1cov)){
    message(j)
    write.csv(as.matrix(c1cov[[j]]), 
              file.path(outdir, 
                        glue('hypergraph_{kk}_P{psx}/covariance_matrix/cat1_clust_covariance_', j, '.csv')), 
              quote = F, row.names = F)
    
    write.csv(as.matrix(c2cov[[j]]), 
              file.path(outdir, 
                        glue('hypergraph_{kk}_P{psx}/covariance_matrix/cat2_clust_covariance_', j, '.csv')), 
              quote = F, row.names = F)
  }
  
}
  

  


