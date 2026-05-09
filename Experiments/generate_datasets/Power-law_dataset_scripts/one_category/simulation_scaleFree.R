# data simulation for testing network inference 

# start with a hypergraph connecting communities,
# simulate data in this order : only shifts in edges (and/or) shifts 
# in connectivities

# generate samples from simulated covariances. compute distance between 
# sample covariance matrices to inform network inference

library(igraph)
library(Matrix)
library(matrixcalc)
library(tidyverse)
library(igraph)
library(glue)

# setwd('/Volumes/umms-ukarvind/vravik/l0-networks/')
setwd("~/Dropbox (University of Michigan)/l0-networks/simulation-compare_l1_vs_l0/r-scripts/")

outdir <- '../simulation_data/'
if(!dir.exists(outdir)){
  dir.create(outdir)
  dir.create(glue(outdir,'/figures'))
  dir.create(glue(outdir,'/true_networks'))
  dir.create(glue(outdir,'/generated_data'))
}


#------------------------
#set structure of graphs
#------------------------
p = 2000 #no of features
M = 10 #no of modules in each network

#simulate base network with a modular structure
cov_list <- list()
omega_list <- list()

for(j in 1:M){
  message('simulating module ', j)
  l <- simulate_powerlaw_networks(p = p/M)
  cov_list[[j]] <- l[["cov"]] ; omega_list[[j]] <- l[["A"]]
}

sigma <- bdiag(cov_list); image(sigma)

#-----------------------------
#derive MST through hypergraph
#-----------------------------
library(ape)
k <- 20 #number of clusters

wts <- matrix(rnorm(n = k*k, mean = 10, sd = 2), nrow = k)
wts <- (wts + t(wts))/2 #make symmetric distance 
rownames(wts) <- paste("cluster", 1:k, sep = "")
colnames(wts) <- paste("cluster", 1:k, sep = "")
m <- ape::mst(wts) #learn mst through hyper-graph
plot(m, graph = 'nsca')

mm <- as.matrix(m)
write.table(mm, file.path(outdir, 'MST-adjacency_matrix.txt'), 
            quote = F, row.names = F, col.names = F)

#given the mst, use bfs to learn path through graph
g_mst <- graph_from_adjacency_matrix(m, mode = "undirected")
plot(g_mst)

pdf(glue('{outdir}/simulation-mst.pdf'), width = 6, height = 6)
par(mar=c(.5,.5,.5,.5))
plot(g_mst, vertex.size=7, vertex.label.font=2, vertex.label = NA,
     vertex.label.cex=.5, layout = layout_with_fr)
dev.off()

## find a leaf node to set as root
degree.mst <- igraph::degree(g_mst)
plot(degree.mst, pch=20)

(dmax <- max(degree.mst))
leaf.nodes <- which(degree.mst == 1)
# leaf.nodes <- names(degree.mst[degree.mst == 1])
root <- sample(leaf.nodes, 1)
message(glue('Root node is {root}'))

mbfs <- bfs(g_mst, root = root, father = T, order = T) 
bfs_order <- as.vector(mbfs$order)
bfs_pred <- as.vector(mbfs$father)

#navigate this graph and simulate network for each cluster
network_list <- vector(mode = "list", length = k)
root_node <- bfs_order[1]
network_list[[root_node]] <- omega_list ## fix this

#-----------------------------
#generate sub-population networks
#-----------------------------
s = 4
delta <- 1 #magnitude of perturbation

## start traversing from a leaf node
## introduce large shifts only when I'm entering new phylogenetic branch
## if degree of parent is 4/5, i introduce new block

#use bfs to traverse the hyper-graph properly
for(node in 2:k){
  
  j <- bfs_order[node]; i <- bfs_pred[j]; 
  message(i, " ", j)  #, " ", wt)
  net_j <- network_list[[i]] #initialize with predecessor
  
  ind <- sample(1:M, s, replace = F) #randomly sample module to perturb
  
  ## comment out this portion if we want to keep structure same
  ## if parent node is a branching point, simulate new module
  if(degree.mst[i] >= 3){
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
figs <- lapply(precision_list, function(x)  image(x))

for(j in 1:length(precision_list)){
    # pdf(glue('{outdir}/figures/cluster_',j,'_precision.pdf'), width = 10)
    # plot(figs[[j]])
    # dev.off()
    write.csv(as.matrix(precision_list[[j]]), 
              file.path(outdir, 
                        glue('true_networks/clust_precision_', j, '.csv')), 
              quote = F, row.names = F)
  
}

# write_rds(precision_list, glue("{outdir}/precisionList.rds"))

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
                      glue('covariance-matrix/clust_covariance_', j, '.csv')), 
            quote = F, row.names = F)
}

#-------------------------------
#simulate data for each cluster
#-------------------------------
#get n = 250 data points for each cluster
n <- 25*p #number of samples

# data.clusters <- vector(mode = "list", length = k)

for(j in 1:k){
  message('simulating cluster ', j)
  #simulate normal variables with this covariance distribution
  cov_matrix <- cov_list[[j]]
  U <- chol(cov_matrix)
  z <- matrix(rnorm(n = n*p), nrow = n)
  dat <- z%*%U
  write.csv(as.matrix(dat), 
            file.path(outdir, 
                      glue('generated_data/simulation_clust_', j, '.csv')),
            quote = F, row.names = F)
}





