library(JGL)
#library(JointNets)

source('Generate_data.R') #setwd('simulation-compare_l1_vs_l0/JGL/')
source('funcs.R')

args <- commandArgs(trailingOnly = TRUE)

L <- as.integer(args[1])
p <- 250
nbyp <- 20
rd_seed <- as.integer(args[2])
N <- p*nbyp
Lamb1 <- as.numeric(args[3])

file_name <- paste0("JGL_CPU",nbyp,"_L",L,".csv")

dir_path <- paste0("Output/AUPRC/", rd_seed)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)

# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("rd_seed,L,p,N,Lamb1,Lamb2,Recall,Precision,F1 score, max Error, norm Error, BIC,runtime", file_name)
}

data <- Generate_Data(L,N,p,rd_seed)
True_Precission <- Read_Simulation_Data(L)

Lamb2 <- c(1) # 0.1

results <- matrix(0, length(Lamb1)*length(Lamb2), 13)
results[, 1] <- rd_seed
results[, 2] <- L
results[, 3] <- p
results[, 4] <- N

count <- 1
for (i in 1:length(Lamb1)){
  for (j in 1:length(Lamb2)){
    print(count)
    lamb1 <- Lamb1[i]
    lamb2 <- Lamb2[j]
    
    start_time <- proc.time()
    result <- JGL(Y=data,penalty="fused",lamb1, lamb2, return.whole.theta=TRUE) # 
    end_time <- proc.time()
    
    Errors <- matrix(0,L,5)
    for (l in 1:L){
      error <- error_cal(result$theta[[l]], True_Precission[[l]])
      Errors[l,] <- c(error$recall, error$precision, error$f1, error$max_error, error$norm_error)
    }
    results[count, 5] <- lamb1
    results[count, 6] <- lamb2
    results[count,7:11] <- colMeans(Errors)
    results[count,12] <- bic_my(data, result$theta)
    results[count,13] <- end_time[3] - start_time[3]
    count <- count + 1
  }
}
write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

