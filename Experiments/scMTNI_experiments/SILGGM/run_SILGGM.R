library(SILGGM)

source('Generate_data_scMTNI.R') 
source('funcs.R')

L <- 3
p <- 65
N <- 200


#W <- diag(1, L) +diag(1, L-1, 1) +diag(1, L-1, -1)


base_path <- "../"
grn_path  <- "../Simulated_GRNs/"
idx <- regulator_idx_R(grn_path)
reg_id <- idx$hsc

file_name <- paste0("SILGGM",N,".csv")

dir_path <- paste0("Output/Table")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)
# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("T,p,N,Lamb1,Lamb2,Recall_1,Precision_1,F1 score-1,AUPRC_1,Recall_2,Precision_2,F1 score-2,AUPRC_2,Recall_3,Precision_3,F1 score-3,AUPRC_3,Recall_avg,Precision_avg,F1 score-avg,AUPRC_avg, BIC,runtime(s),total edges(directed)", file_name)
}
cat("Reading data...\n")
data <- load_scMTNI(N,base_path,TRUE)
True_Precission <- read_scMTNI_prec(grn_path)
cat("Data read successfully.\n")

Lamb1 <- c("NA") # defualt for ackagesqrt(2 * log(p / sqrt(N)) / N)
Lamb2 <- c("NA") # 0.1

results <- matrix(0, length(Lamb1) * length(Lamb2), 12 + 4*L)

results[, 1] <- L
results[, 2] <- p
results[, 3] <- N

masked <- mask_diag(True_Precission)
True_Precission_copy <- masked$true_prec

Scores <- matrix(0, L, 4)


edges_total <- 0
for (l in 1:L){

  start_time <- proc.time()
  result <- SILGGM(data[[l]])
  end_time <- proc.time()
  
  result_copy <- result$precision   # copy
  diag(result_copy) <- 0
  
  edges_total <- edges_total + sum(abs(result_copy) > 1e-5)
  
  error <- error_cal(
    result_copy[reg_id, , drop = FALSE],
    True_Precission_copy[[l]][reg_id, , drop = FALSE]
  )
  auprc <- auprc_from_networks(True_Precission_copy[[l]],
                               result_copy,
                               reg_id, FALSE)
  Scores[l,] <- c(error$recall, error$precision, error$f1, auprc)

}
results[, 4] <- "NA"
results[, 5] <- "NA"

# per-dataset scores start at col 6
results[, 6:(5 + 4*L)] <- as.vector(t(Scores))

# averages next
results[, (6 + 4*L):(9 + 4*L)] <- colMeans(Scores)

# bic and time
results[, 10 + 4*L]  <- "NA"
results[, 11 + 4*L] <- "NA"
results[, 12 + 4*L] <- edges_total


write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

