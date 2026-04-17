library(fasjem)

source('Generate_data_scMTNI.R') #setwd('simulation-compare_l1_vs_l0/FASJEM/')
source('funcs.R')

args <- commandArgs(trailingOnly = TRUE)

L <- 3
p <- 65
N <- 1000

# If argument provided, overwrite default
if (length(args) >= 1) {
  N <- as.numeric(args[1])
}

base_path <- "../"
grn_path  <- "../Simulated_GRNs/"

idx <- regulator_idx_R(grn_path)
reg_id <- idx$hsc

file_name <- paste0("FASJEM_",N,".csv")

dir_path <- paste0("Output/AUPRC/")
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, file_name)

# Check if the file exists
if (!file.exists(file_name)) {
  writeLines("T,p,N,Lamb,Epsilon0,Recall_1,Precision_1,F1 score-1,Recall_2,Precision_2,F1 score-2,Recall_3,Precision_3,F1 score-3,Recall_avg,Precision_avg,F1 score-avg, BIC,runtime(s), total edges(directed)", file_name)
}
cat("Reading data...\n")
data <- load_scMTNI(N,base_path,TRUE)
True_Precission <- read_scMTNI_prec(grn_path)
cat("Data read successfully.\n")

Lamb <- 10^seq(3, -5, length.out = 200) # 5
Epsilon0 <- c(0.1) # default is 0.1

results <- matrix(0, length(Lamb) * length(Epsilon0), 11 + 3*L)

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
    
    masked <- mask_diag(True_Precission, result)
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
    results[count, 5] <- epsilon0
    
    # per-dataset scores start at col 6
    results[count, 6:(5 + 3*L)] <- as.vector(t(Scores))
    
    # averages next
    results[count, (6 + 3*L):(8 + 3*L)] <- colMeans(Scores)
    
    # bic and time
    results[count, 9 + 3*L]  <- bic_my(data, result)
    results[count, 10 + 3*L] <- end_time[3] - start_time[3]
    results[count, 11 + 3*L] <- edges_total
    
    count <- count + 1
  }
}
write.table(results, file = file_name, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
