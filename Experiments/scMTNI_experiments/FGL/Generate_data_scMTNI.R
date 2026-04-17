load_scMTNI <- function(n, base_path, transpose = FALSE) {
  if (n == 2000) {
    files <- c(paste0(base_path, "hsc.table"),
               paste0(base_path, "cmp.table"),
               paste0(base_path, "gmp.table"))
  } else if (n == 1000) {
    files <- c(paste0(base_path, "hsc_n1000.table"),
               paste0(base_path, "cmp_n1000.table"),
               paste0(base_path, "gmp_n1000.table"))
  } else if (n == 200) {
    files <- c(paste0(base_path, "hsc_n200.table"),
               paste0(base_path, "cmp_n200.table"),
               paste0(base_path, "gmp_n200.table"))
  } else {
    stop("n must be one of {200, 1000, 2000}")
  }

  data <- vector("list", 3)

  for (i in 1:3) {
    # Read everything as raw text table (no headers interpreted)
    M <- read.table(files[i], sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE, check.names = FALSE)

    # Drop first row and first column (Python: df.iloc[1:,1:])
    A <- as.matrix(M[-1, -1, drop = FALSE])

    # Convert to numeric safely
    A <- apply(A, 2, as.numeric)

    if (transpose) A <- t(A)

    data[[i]] <- A

    cat(sprintf("Loaded %s: %dx%d (%s)\n",
                files[i], nrow(A), ncol(A),
                if (transpose) "cells x genes" else "genes x cells"))
  }

  return(data)
}

read_scMTNI_prec <- function(grn_path) {
  t=3
  names <- c("hsc", "cmp", "gmp")
  true_networks <- vector("list", t)
  
  for (i in 1:t) {
    
    path <- file.path(grn_path, paste0(names[i], "_precision_matrix.csv"))
    # Skip first row (NumHeaderLines = 1)
    M <- read.csv(path, skip = 1, header = FALSE)
    
    # Drop first column (MATLAB: M(:,2:end))
    A <- as.matrix(M[, -1])
    
    true_networks[[i]] <- A
  }
  
  return(true_networks)
}

regulator_idx_R <- function(grn_path) {
  
  # Cell types
  cell_types <- c("hsc", "cmp", "gmp")
  
  idx <- list()
  
  for (ct in cell_types) {
    
    infile <- paste0(grn_path, ct, "_regulators_idx.txt")
    
    v <- scan(infile, what = numeric(), quiet = TRUE)
    
    v <- as.vector(v)
    
    # Convert 0-based -> 1-based indexing
    v <- v + 1
    
    idx[[ct]] <- v
  }
  
  return(idx)
}
