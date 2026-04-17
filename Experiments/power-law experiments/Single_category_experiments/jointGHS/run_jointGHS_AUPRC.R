library(jointGHS)
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
file_name <- paste0("jointGHS_CPU",nbyp,"_L",L,".csv")

dir_path <- paste0("Output/AUPRC/", rd_seed)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)

# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("rd_seed,L,p,N,Lamb1,Lamb2,tau_sq,Recall,Precision,F1 score, max Error, norm Error, BIC,runtime", file_name)
}

data <- Generate_Data(L,N,p,rd_seed)
True_Precission <- Read_Simulation_Data(L)

Epsilon <- c(0.001) # 0.44
AIC_EP_list <- c(0.1) # 0.1
Tau_sq_list <- 10^seq(2, -6, length.out = 200)

results <- matrix(0, length(Epsilon)*length(AIC_EP_list)* length(Tau_sq_list), 14)
results[, 1] <- rd_seed
results[, 2] <- L
results[, 3] <- p
results[, 4] <- N

count <- 1
for (i in 1:length(Epsilon)){
  for (j in 1:length(AIC_EP_list)){
    for (k in 1:length(Tau_sq_list)) {
      print(count)

      epsilon <- Epsilon[i]
      aic_epsilon <- AIC_EP_list[j]
      tau_sq <- Tau_sq_list[k]
      
    
      start_time <- proc.time()
      result <- jointGHS(X=data,epsilon = epsilon, AIC_eps = aic_epsilon, tau_sq=tau_sq*rep(1, L),fix_tau = TRUE, AIC_selection = FALSE) # 
      end_time <- proc.time()

      Errors <- matrix(0,L,5)
      for (l in 1:L){
        error <- error_cal(cov2cor(result$theta[[l]]), True_Precission[[l]])
        Errors[l,] <- c(error$recall, error$precision, error$f1, error$max_error, error$norm_error)
      }
      results[count, 5] <- epsilon
      results[count, 6] <- aic_epsilon
      results[count, 7] <- tau_sq
      results[count,8:12] <- colMeans(Errors)
      results[count,13] <- bic_my(data, result$theta)
      results[count,14] <- end_time[3] - start_time[3]
      
      count <- count + 1
    }
  }
}
write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
