library(psych)
error_cal <- function(A, invSigma){
  Sparse_estimation <- abs(A)>1e-5
  Sparse_true <- abs(invSigma)>1e-5
  
  TP = sum(Sparse_true*Sparse_estimation)
  FP = sum(Sparse_estimation)-TP
  FN = sum(Sparse_true) - TP
  recall = TP/(TP+FN)
  precision = TP/(TP+FP)
  f1 = 2*recall*precision/(recall+precision)
  
  E = A-invSigma
  nnz_E = E[E!=0]
  nnz_invSigma = invSigma[invSigma!=0]
  max_error = max(abs(nnz_E))
  norm_error = norm(nnz_E,'2')/norm(nnz_invSigma,'2')
  
  result <- list(recall, precision, f1, max_error, norm_error)
  names(result) <- c("recall","precision","f1","max_error","norm_error")
  return(result)
}

mask_diag <- function(true_prec) {
  
  L <- length(true_prec)
  
  true_out <- vector("list", L)

  
  for (l in 1:L) {
    
    Tmat <- true_prec[[l]]

    
    diag(Tmat) <- 0

    
    true_out[[l]] <- Tmat

  }
  
  return(list(true_prec = true_out))
}

error_cal_diff <- function(theta, true){
  n = length(theta)
  results <- matrix(nrow=n*(n-1)/2, ncol=5)
  count <- 1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      theta_diff <- theta[[i]]-theta[[j]]
      true_diff <- true[[i]]-true[[j]]
      result_ij <- error_cal(theta_diff, true_diff)
      results[count,1:5] <- c(result_ij$recall, result_ij$precision, result_ij$f1, result_ij$max_error, result_ij$norm_error)
      count = count + 1
    }
  }
  names(results) <- c("recall","precision","f1","max_error","norm_error")
  return(results)
}


bic_my <- function(sample, theta){
  bic_total <- 0
  total_nnz <- 0
  for (i in 1:length(sample)) {
    temp_t <- theta[[i]]
    x <- sample[[i]]
    nnz <- sum(temp_t[upper.tri(temp_t)] != 0)
    total_nnz <- total_nnz + nnz
    bic_total <- bic_total + dim(sample[[i]])[1]*(-determinant(temp_t)$modulus[1]+tr(cov(x)%*%temp_t)) + nnz*log(dim(sample[[i]])[1])+4*nnz*dim(sample[[i]])[2]
  }
  #print(bic_total)
  #print(total_nnz)
  return(bic_total)
}

bic_jgl <- function(sample, theta){
  bic_total <- 0
  total_nnz <- 0
  for (i in 1:length(sample)) {
    temp_t <- theta[[i]]
    x <- sample[[i]]
    nnz <- sum(temp_t[upper.tri(temp_t)] != 0)
    total_nnz <- total_nnz + nnz
    bic_total <- bic_total + dim(sample[[i]])[1]*(-determinant(temp_t)$modulus[1]+tr(cov(x)%*%temp_t)) + nnz*2
  }
  print(bic_total)
  print(total_nnz)
  return(bic_total)
}