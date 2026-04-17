library(fasjem)

source('Generate_data.R') #setwd('simulation-compare_l1_vs_l0/FASJEM/')
source('funcs.R')

args <- commandArgs(trailingOnly = TRUE)

L <- 10
p <- 200
N <- 2000

# If argument provided, overwrite default
if (length(args) >= 1) {
  N <- as.numeric(args[1])
}
no_regs <- 30
reg_id <- 1:no_regs

file_name <- paste0("FASJEM_",N,".csv")

dir_path <- paste0("Output/AUPRC/")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)

# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("T,p,N,Lamb,Epsilon0,Recall,Precision,F1 score, BIC,runtime(s), total edges(directed)", file_name)
}
cat("Reading data...\n")
data <- load_dataset(no_regs,p,L)          # new data loader (8 clusters)
True_Precission <- read_prec(no_regs,p,L)  # new precision loader
cat("Data read successfully.\n")

Lamb <- 10^seq(3, -5, length.out = 200) # 5
Epsilon0 <- c(0.1) # default is 0.1

results <- matrix(0, length(Lamb) * length(Epsilon0), 11)

results[, 1] <- T
results[, 2] <- p
results[, 3] <- N


count <- 1
for (i in 1:length(Lamb)){
  for (j in 1:length(Epsilon0)){
    print(count)
    lamb1 <- Lamb[i]
    epsilon0 <- Epsilon0[j]
    
    start_time <- proc.time()
    result <- fasjem(data,method="fasjem-g", lambda = lamb1*sqrt(log(L*p)/(N)), epsilon = epsilon0) # 
    end_time <- proc.time()
    
    result_copy <- result
    True_Precission_copy <- True_Precission
    
    for (l in 1:L) {
      for (i in reg_id) {
        result_copy[[l]][i, i] <- 0
        True_Precission_copy[[l]][i, i] <- 0
      }
    }
    
    TP <- 0
    FP <- 0
    FN <- 0
    edges_total <- 0
    
    for (l in 1:L) {
    
      pred <- result_copy[[l]][reg_id, , drop = FALSE]
      truth <- True_Precission_copy[[l]][reg_id, , drop = FALSE]
    
      pred_bin  <- abs(pred) > 1e-5
      truth_bin <- abs(truth) > 1e-5
    
      TP <- TP + sum(pred_bin & truth_bin)
      FP <- FP + sum(pred_bin & !truth_bin)
      FN <- FN + sum(!pred_bin & truth_bin)
    
      edges_total <- edges_total + sum(pred_bin)
    }
    
    # Global metrics
    precision <- TP / (TP + FP )
    recall    <- TP / (TP + FN )
    f1        <- 2 * precision * recall / (precision + recall)
    
    results[count, 4] <- lamb1
    results[count, 5] <- epsilon0
    results[count, 6:8] <- c(recall, precision, f1)
    
    # bic and time
    results[count, 9]  <- bic_my(data, result)
    results[count, 10] <- end_time[3] - start_time[3]
    results[count, 11] <- edges_total
    
    count <- count + 1
  }
}
write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
