library(JGL)
#library(JointNets)

source('Generate_data.R') #setwd('simulation-compare_l1_vs_l0/JGL/')
source('funcs.R')

L <- 10
p <- 100
args <- commandArgs(trailingOnly = TRUE)
no_regs <- 25
reg_id <- 1:no_regs

if(length(args) == 0){
  N <- 2500   
} else {
  N <- as.numeric(args[1])
}


file_name <- paste0("FGL_",N,".csv")

dir_path <- paste0("Output/AUPRC")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)
# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("T,p,N,Lamb1,Lamb2,Recall,Precision,F1 score, BIC,runtime(s), total edges(directed)", file_name)
}
cat("Reading data...\n")
data <- load_dataset(no_regs,p,L)          # new data loader (8 clusters)
True_Precission <- read_prec(no_regs,p,L)  # new precision loader
cat("Data read successfully.\n")

Lamb1 <- 10^seq(2, -6, length.out = 50)
Lamb2 <- c(1) # 0.1

results <- matrix(0, length(Lamb1)*length(Lamb2), 11)
results[, 1] <- L
results[, 2] <- p
results[, 3] <- N

count <- 1
for (i in 1:length(Lamb1)){
  for (j in 1:length(Lamb2)){
    print(count)
    lamb1 <- Lamb1[i]
    lamb2 <- Lamb2[j]
    
    start_time <- proc.time()
    result <- JGL(Y=data,penalty="fused",lamb1, lamb2, return.whole.theta=TRUE) # 
    end_time <- proc.time()
    
    result_copy <- result$theta
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
    results[count, 5] <- lamb2
    results[count, 6:8] <- c(recall, precision, f1)
    
    # bic and time
    results[count, 9]  <- bic_my(data, result$theta)
    results[count, 10 ] <- end_time[3] - start_time[3]
    results[count, 11] <- edges_total
    
    count <- count + 1
    line <- c(
      T,
      p,
      N,
      lamb1,
      lamb2,
      recall,
      precision,
      f1,
      bic_my(data, result$theta),
      end_time[3] - start_time[3],
      edges_total
        )
        
    write.table(
      t(line),
      file = file_name,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
)
  }
}
#write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
