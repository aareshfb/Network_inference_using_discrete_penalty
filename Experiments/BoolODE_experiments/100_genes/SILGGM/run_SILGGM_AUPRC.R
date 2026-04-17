library(SILGGM)
library(PRROC)

source("Generate_data.R")
source("funcs.R")

L <- 10
p <- 100
no_regs <- 25
reg_id <- 1:no_regs
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  N <- 2000
} else {
  N <- as.numeric(args[1])
}

dir_path <- "Output/AUPRC"
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

file_name <- file.path(dir_path, paste0("SILGGM_", N, ".csv"))

# header
if (!file.exists(file_name)) {
  writeLines(
    "population,p,N,Recall,Precision,F1,AUPRC,runtime_s,total_edges_directed",
    file_name
  )
}

cat("Reading data...\n")
data <- load_dataset(no_regs,p,L)          # list of length L
True_Precission <- read_prec(no_regs,p,L)  # list of length L
cat("Data read successfully.\n")

# store per-population results
results <- data.frame(
  population = as.character(1:L),
  p = rep(p, L),
  N = rep(N, L),
  Recall = rep(NA_real_, L),
  Precision = rep(NA_real_, L),
  F1 = rep(NA_real_, L),
  AUPRC = rep(NA_real_, L),
  runtime_s = rep(NA_real_, L),
  total_edges_directed = rep(NA_real_, L),
  stringsAsFactors = FALSE
)

# for global AUPRC
all_scores <- c()
all_truth  <- c()

# for global recall/precision/f1 if you want them too
TP_global <- 0
FP_global <- 0
FN_global <- 0

for (l in 1:L) {
  cat("Processing population", l, "\n")

  start_time <- proc.time()
  result <- SILGGM(data[[l]])
  end_time <- proc.time()

  runtime_l <- (end_time - start_time)[3]

  # estimated precision
  est <- result$precision
  diag(est) <- 0

  # true precision
  truth <- True_Precission[[l]]
  diag(truth) <- 0

  # restrict to regulator rows only
  est_sub <- est[reg_id, , drop = FALSE]
  truth_sub <- truth[reg_id, , drop = FALSE]

  # total predicted directed edges
  edges_total <- sum(abs(est_sub) > 1e-5)

  # per-population recall / precision / f1
  error <- error_cal(est_sub, truth_sub)

  # per-population auprc
  auprc_l <- auprc_from_networks(truth_sub, est_sub, reg_id, FALSE)

  results$Recall[l] <- error$recall
  results$Precision[l] <- error$precision
  results$F1[l] <- error$f1
  results$AUPRC[l] <- auprc_l
  results$runtime_s[l] <- runtime_l
  results$total_edges_directed[l] <- edges_total

  # write row immediately
  write.table(
    results[l, , drop = FALSE],
    file = file_name,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )

  # accumulate global TP / FP / FN
  pred_bin  <- abs(est_sub) > 1e-5
  truth_bin <- abs(truth_sub) > 1e-5

  TP_global <- TP_global + sum(pred_bin & truth_bin)
  FP_global <- FP_global + sum(pred_bin & !truth_bin)
  FN_global <- FN_global + sum(!pred_bin & truth_bin)

  # accumulate for global AUPRC
  all_scores <- c(all_scores, as.vector(abs(est_sub)))
  all_truth  <- c(all_truth,  as.vector(abs(truth_sub) > 1e-5))
}

# global precision / recall / f1
precision_global <- if ((TP_global + FP_global) == 0) NA_real_ else TP_global / (TP_global + FP_global)
recall_global    <- if ((TP_global + FN_global) == 0) NA_real_ else TP_global / (TP_global + FN_global)
f1_global <- if (is.na(precision_global) || is.na(recall_global) || (precision_global + recall_global) == 0) {
  NA_real_
} else {
  2 * precision_global * recall_global / (precision_global + recall_global)
}

# global AUPRC
if (sum(all_truth == 1) == 0 || sum(all_truth == 0) == 0) {
  global_auprc <- NA_real_
} else {
  pr <- PRROC::pr.curve(
    scores.class0 = all_scores[all_truth == 1],
    scores.class1 = all_scores[all_truth == 0],
    curve = TRUE
  )
  global_auprc <- pr$auc.integral
}

# summary/global row
summary_row <- data.frame(
  population = "global",
  p = p,
  N = N,
  Recall = recall_global,
  Precision = precision_global,
  F1 = f1_global,
  AUPRC = global_auprc,
  runtime_s = sum(results$runtime_s, na.rm = TRUE),
  total_edges_directed = sum(results$total_edges_directed, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.table(
  summary_row,
  file = file_name,
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE
)

cat("\nPer-population results:\n")
print(results)

cat("\nGlobal summary:\n")
print(summary_row)