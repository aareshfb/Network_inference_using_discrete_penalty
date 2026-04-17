library(SILGGM)
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

dir_path <- paste0("Output/Table/", rd_seed)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)

# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("rd_seed,L,p,N,Lamb1,Lamb2,Recall,Precision,F1 score, max Error, norm Error, BIC (last pop),runtime", file_name)
}

data <- Generate_Data(L,N,p,rd_seed)
True_Precission <- Read_Simulation_Data(L)

Lamb1 <- c(sqrt(2*log(p/sqrt(N))/N),1e-1,1e-2,1e-3,1e-4,1e-5) # 0.44
Lamb2 <- c("NA") # 0.1

results <- matrix(0, length(Lamb1)*length(Lamb2), 13)
results[, 1] <- rd_seed
results[, 2] <- L
results[, 3] <- p
results[, 4] <- N

total_time <- 0
count <- 1
for (i in 1:length(Lamb1)){
  for (j in 1:length(Lamb2)){
    Errors <- matrix(0,L,5)
    for (l in 1:L){
      print(count)
      lamb1 <- Lamb1[i]
      lamb2 <- Lamb2[j]
    
      start_time <- proc.time()
      result <- SILGGM(data[[l]], lambda =lamb1)
      end_time <- proc.time()
      total_time <- total_time +end_time[3] - start_time[3]
      
      error <- error_cal(result$partialCor, True_Precission[[l]])
      Errors[l,] <- c(error$recall, error$precision, error$f1, error$max_error, error$norm_error)
    }
    results[count, 5] <- lamb1
    results[count, 6] <- lamb2
    results[count,7:11] <- colMeans(Errors)
    results[count,12] <- "NA" #bic_my(data[[l]], result$theta)
    results[count,13] <- total_time
    count <- count + 1
    
  }
}
write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

