library(tidyverse)
library(igraph)
library(Matrix)
library(matrixcalc)
library(igraph)
library(glue)


p_root <- .3 #prob of interaction for active TFs in root node

p_reg_loss <- .1 #chance of population losing one of active regulators
p_reg_gain <- .2  # chance of population gaining a new regulators
p_gain <- .2 # number of edges added w/ new regulator
p_maintain <- .8

p_inhibitor <- .3 # prob of interaction being negative

n_regs <- 30
n_genes <- 200

reg_names <- paste0('reg', 1:n_regs)
gene_names <- c(reg_names, paste0('gene',1:(n_genes - n_regs)))

k_active <- 8
## initialize theta with only select active regulators
(active_regs <- sample(1:n_regs, size=k_active, replace = F))
## create initial conditions file
# ics <- data.frame(Genes = character(), Values = character(), stringsAsFactors = F)
# ics <- rbind(ics, data.frame(Genes = paste(reg_names[active_regs],collapse=',')  ) )

theta <- matrix(0, nrow = n_regs, ncol = n_genes)
## add non-zero prob of edges for active TFs
for(r in active_regs){
  raw_probs <- runif(n_genes)
  theta[r, raw_probs < p_root] <- 1
}

nze <- which(theta > 0)
idx <- sample(nze, p_inhibitor*length(nze), replace = F)
theta[idx] <- -1
rownames(theta) <- reg_names
colnames(theta) <- gene_names

par(mfrow=c(1,1))
image(theta)
hist(theta)

#-----------------------------
#derive MST through hypergraph
#-----------------------------
library(ape)
k <- 8 #number of clusters

set.seed(33)
hg <- matrix(rnorm(n = k*k, mean = 10, sd = 3), nrow = k)
hg <- (hg + t(hg))/2 #make symmetric distance 
rownames(hg) <- paste("cluster", 1:k, sep = ""); colnames(hg) <- paste("cluster", 1:k, sep = "")
m <- ape::mst(hg) #learn mst through hyper-graph
g_mst <- graph_from_adjacency_matrix(m, mode = "undirected")

pdf('spanning-tree-arbitrary.pdf')
plot(g_mst, layout = layout_with_fr)
dev.off()

## find a leaf node to set as root
degree.mst <- igraph::degree(g_mst)
(dmax <- max(degree.mst))
leaf.nodes <- which(degree.mst == 1)
# root <- names(V(g_mst))[sample(leaf.nodes, 1)]
root <- sample(leaf.nodes, 1)
message(glue('Root node is {root}'))

#given the mst, use bfs to learn path through graph
mbfs <- bfs(g_mst, root = root, father = T,
            order = T) 

(bfs_order <- as.vector(mbfs$order))
(bfs_pred <- as.vector(mbfs$father))

#navigate this graph and simulate network for each cluster
network_list <- vector(mode = "list", length = k)
network_list[[root]] <- theta
#-----------------------------
## navigate the MST and update edges in each population
#-----------------------------
for(node in 2:k){
  
  j <- bfs_order[node]; i <- bfs_pred[j]; 
  message(i, " ", j) 
  net_i <- network_list[[i]]
  net_j <- net_i #initialize with predecessor
  
  
  ## add edges for existing regulators
  ar <- which(rowSums(net_i) != 0) # find active regulators in network
  
  ze <- which(net_i[ar, ] == 0)
  ix_gain <- sample(ze, size = p_gain*length(ze), replace = F)
  inh_gain <- sample(ix_gain, size = p_inhibitor*length(ix_gain), replace = F)
  net_j[ar, ][ix_gain] <- 1
  net_j[ar, ][inh_gain] <- -1
  
  ## maintaining / losing links from parent
  ix_lose <- sample(which(net_i != 0), size = (1-p_maintain)*sum(abs(net_i)), replace = F)
  net_j[ix_lose] <- 0
  
  ## update regulatory edges
  if(runif(1) < p_reg_gain){
    gain_reg <- sample(setdiff(1:n_regs, ar), 1)
    idx <- which(runif(n_genes) < p_gain)
    inh <- sample(idx, size=length(idx)*p_inhibitor, replace = F)
    net_j[gain_reg, idx] <- 1
    net_j[gain_reg, inh] <- -1
  }
  
  if(runif(1) < p_reg_loss){
    lose_reg <- sample(ar, 1)
    net_j[lose_reg, ] <- 0
  }
  
  network_list[[j]] <- net_j
}
  

par(mfrow=c(3,4))
sapply(network_list, image)

par(mfrow = c(3,4))
sapply(network_list, hist)
sapply(network_list, function(z) sum(rowSums(z) != 0))

#-----------------------------
## write networks to file
#-----------------------------
## include self-loops for all regulators
network_list <- lapply(network_list, function(z) {diag(z) <- 1; z})
sapply(network_list, diag)

par(mfrow=c(3,4))
sapply(network_list, function(theta) hist(colSums(abs(theta)))) ##number of regulators per gene
sapply(1:k, function(i) write.csv(network_list[[i]], file = glue('network-{i}.csv'), quote = F))

#-----------------------------
#-----------------------------
## translate to boolean rules
# Iterate through each cluster
dir.create('sim_results')
for (k_idx in 1:k) {
  adj_mat <- network_list[[k_idx]]
  rules_df <- data.frame(Gene = character(), Rule = character(), stringsAsFactors = F)
  
  for (j in 1:ncol(adj_mat)) {
    gene_name <- gene_names[j]
    
    # Find row indices for activators (1) and inhibitors (-1)
    act_idx <- which(adj_mat[, j] == 1)
    inh_idx <- which(adj_mat[, j] == -1)
    
    # Get the actual names of the regulators
    activators <- reg_names[act_idx]
    inhibitors <- reg_names[inh_idx]
    
    # --- Construct the Logic String ---
    act_str <- if(length(activators) > 0) paste0("(", paste(activators, collapse = " or "), ")") else ""
    inh_str <- if(length(inhibitors) > 0) paste0("not (", paste(inhibitors, collapse = " or "), ")") else ""
    
    final_rule <- ""
    
    if (act_str != "" && inh_str != "") {
      final_rule <- paste(act_str, "and", inh_str)
    } else if (act_str != "") {
      final_rule <- act_str
    } else if (inh_str != "") {
      final_rule <- inh_str
    } else{
      next
      # final_rule <- gene_name
    }
    
    
    rules_df <- rbind(rules_df, data.frame(Gene = gene_name, Rule = final_rule))  
    
  }
  
  # Save as Tab-Separated File in cluster-specific folder
  # BoolODE expects a header: Gene [TAB] Rule
  file_name <- paste0("rules_cluster_", k_idx, ".txt")
  write.table(rules_df, file = file_name, sep = "\t", quote = F, row.names = F)
  message(paste("Saved:", file_name))
}



