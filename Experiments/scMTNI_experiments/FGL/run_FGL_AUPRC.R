library(JGL)
#library(JointNets)

source('Generate_data_scMTNI.R') #setwd('simulation-compare_l1_vs_l0/JGL/')
source('funcs.R')

L <- 3
p <- 65
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0){
  N <- 1000   
} else {
  N <- as.numeric(args[1])
}

base_path <- "../"
grn_path  <- "../Simulated_GRNs/"
idx <- regulator_idx_R(grn_path)
reg_id <- idx$hsc

file_name <- paste0("FGL_",N,".csv")

dir_path <- paste0("Output/AUPRC")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)
# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("T,p,N,Lamb1,Lamb2,Recall_1,Precision_1,F1 score-1,Recall_2,Precision_2,F1 score-2,Recall_3,Precision_3,F1 score-3,Recall_avg,Precision_avg,F1 score-avg, BIC,runtime(s), total edges(directed)", file_name)
}
cat("Reading data...\n")
data <- load_scMTNI(N,base_path,TRUE)
True_Precission <- read_scMTNI_prec(grn_path)
cat("Data read successfully.\n")

Lamb1 <- 10^seq(2, -6, length.out = 200)
Lamb2 <- c(1) # 0.1

results <- matrix(0, length(Lamb1)*length(Lamb2), 11+3*L)
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
    
    masked <- mask_diag(True_Precission, result$theta)
    True_Precission_copy <- masked$true_prec
    result_copy <- masked$inferred
    edges_total <- 0
    
    Scores <- matrix(0, L, 3)
    for (l in 1:L){
      edges_total <- edges_total + sum(abs(result_copy[[l]]) > 1e-5)
      error <- error_cal(
        result_copy[[l]][reg_id, , drop = FALSE],
        True_Precission_copy[[l]][reg_id, , drop = FALSE]
      )
      Scores[l,] <- c(error$recall, error$precision, error$f1)
    }
    
    results[count, 4] <- lamb1
    results[count, 5] <- lamb2
    
    # per-dataset scores start at col 6
    results[count, 6:(5 + 3*L)] <- as.vector(t(Scores))
    
    # averages next
    results[count, (6 + 3*L):(8 + 3*L)] <- colMeans(Scores)
    
    # bic and time
    results[count, 9 + 3*L]  <- bic_my(data, result$theta)
    results[count, 10 + 3*L] <- end_time[3] - start_time[3]
    results[count, 11 + 3*L] <- edges_total
    
    count <- count + 1
  }
}
write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
