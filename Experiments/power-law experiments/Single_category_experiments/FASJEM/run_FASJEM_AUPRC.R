library(fasjem)

source('Generate_data.R') #setwd('simulation-compare_l1_vs_l0/FASJEM/')
source('funcs.R')

args <- commandArgs(trailingOnly = TRUE)

L <- as.integer(args[1])
p <- as.integer(args[2])
nbyp <- as.numeric(args[3])
rd_seed <- as.integer(args[4])
N <- p*nbyp

#p <- 250
#N <- 250*20
#L <- 20 #3,5,20,20,50,100
#nbyp <-N/p
file_name <- paste0("FASJEM_CPU",nbyp,"_L",L,".csv")

dir_path <- paste0("Output/AUPRC/", rd_seed)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)

# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("rd_seed,L,p,N,Lamb,Epsilon0,Recall,Precision,F1 score, max Error, norm Error, BIC,runtime(s)", file_name)
}
cat("Reading data...\n")
data <- Generate_Data(L,N,p,rd_seed)
True_Precission <- Read_Simulation_Data(L)
cat("Data read successfully.\n")

Lamb <- 10^seq(2, -4, length.out = 200) # 5
Epsilon0 <- c(0.1) # 1

results <- matrix(0, length(Lamb)*length(Epsilon0), 13)
results[, 1] <- rd_seed
results[, 2] <- L
results[, 3] <- p
results[, 4] <- N

count <- 1
for (i in 1:length(Lamb)){
  for (j in 1:length(Epsilon0)){
    print(count)
    lamb1 <- Lamb[i]
    epsilon0 <- Epsilon0[j]
    
    start_time <- proc.time()
    result <- fasjem(data,method="fasjem-g", lambda = lamb1*sqrt(log(L*p)/(N)), epsilon = epsilon0) # 
    end_time <- proc.time()
    
    Errors <- matrix(0,L,5)
    for (l in 1:L){
      error <- error_cal(result[[l]], True_Precission[[l]])
      Errors[l,] <- c(error$recall, error$precision, error$f1, error$max_error, error$norm_error)
    }
    results[count, 5] <- lamb1
    results[count, 6] <- epsilon0
    results[count,7:11] <- colMeans(Errors)
    results[count,12] <- bic_my(data, result)
    results[count,13] <- end_time[3] - start_time[3]
    count <- count + 1
  }
}
write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

