## check if networks are PD 
library(tidyverse)
library(matrixcalc)

setwd('/Volumes/umms-ukarvind/shared_data/l0-networks/simulations-category/simulation_data/hypergraph_20/')
setwd('true_networks/')

plist <- list.files(pattern = 'cat1_*')
plist
cov_list <- lapply(plist, read_csv)
cov_list <- lapply(cov_list, as.matrix)
lapply(cov_list, is.positive.definite)
