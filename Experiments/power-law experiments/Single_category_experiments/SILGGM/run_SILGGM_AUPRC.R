library(SILGGM)
library(PRROC)
#library(JointNets)

source('Generate_data.R') #setwd('simulation-compare_l1_vs_l0/jointGHS/')
source('funcs.R')

args <- commandArgs(trailingOnly = TRUE)

L <- as.integer(args[1])
p <- as.integer(args[2])
nbyp <- as.numeric(args[3])
rd_seed <- as.integer(args[4])
N <- p*nbyp

#p <- 250
#N <- 250*20
#L <- 10
#nbyp <-N/p

#W_fname <- paste0("/nfs/turbo/umms-ukarvind/shared_data/l0-networks/simulations-hypergraphs/simulation_data/hypergraph_", L, "/MST.txt")
#W <- as.matrix(read.table(W_fname))

file_name <- paste0("SILGGM_CPU",nbyp,"_L",L,".csv")

dir_path <- paste0("Output/AUPRC/", rd_seed)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)

# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("rd_seed,L,p,N,Lamb1,Lamb2,Recall,Precision,F1 score, max Error, norm Error, BIC (last pop),runtime,AUPRC", file_name)
}

data <- Generate_Data(L,N,p,rd_seed)
True_Precission <- Read_Simulation_Data(L)

Lamb1 <- c(sqrt(2*log(p/sqrt(N))/N)) # 0.44
Lamb2 <- c("NA") # 0.1

results <- matrix(0, length(Lamb1)*length(Lamb2), 14)
results[, 1] <- rd_seed
results[, 2] <- L
results[, 3] <- p
results[, 4] <- N

total_time <- 0
count <- 1
for (i in 1:length(Lamb1)){
  for (j in 1:length(Lamb2)){
    Errors <- matrix(0,L,5)
    
    # collect all populations here
    all_scores <- list()
    all_truth  <- list()
    
    for (l in 1:L){
      print(count)
      lamb1 <- Lamb1[i]
      lamb2 <- Lamb2[j]
    
      start_time <- proc.time()
      result <- SILGGM(data[[l]])
      end_time <- proc.time()
      total_time <- total_time +end_time[3] - start_time[3]
      
      est_mat <- result$partialCor
      true_mat <- True_Precission[[l]]

      diag(est_mat) <- 0
      diag(true_mat) <- 0
      
      error <- error_cal(est_mat, true_mat)
      Errors[l,] <- c(error$recall, error$precision, error$f1, error$max_error, error$norm_error)
      
      scores <- est_mat
      truth  <- as.integer(true_mat != 0)

      all_scores[[l]] <- scores
      all_truth[[l]]  <- truth
    }
    
    all_scores <- unlist(all_scores)
    all_truth  <- unlist(all_truth)
    
    pr <- pr.curve(
      scores.class0 = all_scores[all_truth == 1],
      scores.class1 = all_scores[all_truth == 0],
      curve = FALSE
    )
    
    AUPRC <- pr$auc.integral
    
    results[count, 5] <- lamb1
    results[count, 6] <- lamb2
    results[count,7:11] <- colMeans(Errors)
    results[count,12] <- "NA" #bic_my(data[[l]], result$theta)
    results[count,13] <- total_time
    results[count,14] <- AUPRC
    
    count <- count + 1
    
  }
}
write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
