load_dataset <- function(no_regs,p,L,transpose = TRUE) {

  base_path <- "../Expression_data"
  set.seed(0)

  expected_regs <- paste0("reg", 1:no_regs)
  expected_genes <- paste0("gene", 1:(p-no_regs))

  natural_key <- function(x) {
    m <- regexec("([A-Za-z]+)([0-9]+)$", x)
    parts <- regmatches(x, m)
    sapply(parts, function(p1) {
      if (length(p1) == 3) {
        return(sprintf("%s%03d", p1[2], as.integer(p1[3])))
      } else {
        return(x)
      }
    })
  }

  ordered_regs <- expected_regs[order(natural_key(expected_regs))]
  ordered_genes <- expected_genes[order(natural_key(expected_genes))]
  ordered_rows <- c(ordered_regs, ordered_genes)

  data <- vector("list", L)

  for (cluster_idx in 1:L) {

    path <- file.path(
      base_path,
      sprintf("cluster_%d_sparse_expression.csv", cluster_idx)
    )

    df <- read.csv(path, row.names = 1, check.names = FALSE)
    df[] <- lapply(df, function(col) as.numeric(as.character(col)))

    missing_rows <- setdiff(ordered_rows, rownames(df))

    if (length(missing_rows) > 0) {
      noise <- matrix(
        rnorm(length(missing_rows) * ncol(df), mean = 0, sd = 1e-5),
        nrow = length(missing_rows),
        ncol = ncol(df)
      )
      rownames(noise) <- missing_rows
      colnames(noise) <- colnames(df)

      df <- rbind(df, noise)
    }

    df <- df[ordered_rows, , drop = FALSE]

    arr <- as.matrix(df)

    if (transpose) arr <- t(arr)

    data[[cluster_idx]] <- arr

    cat(sprintf(
      "Loaded %s: added missing=%d, final shape=%dx%d\n",
      path, length(missing_rows), nrow(arr), ncol(arr)
    ))
    cat(sprintf("dim = %d x %d\n", nrow( arr), ncol( arr)))
  }

  return(data)
}

read_prec <- function(no_regs,p,L) {

  base_path <- "../Ground_truth"

  expected_regs <- paste0("reg", 1:no_regs)
  expected_genes <- paste0("gene", 1:(p-no_regs))

  natural_key <- function(x) {
    m <- regexec("([A-Za-z]+)([0-9]+)$", x)
    parts <- regmatches(x, m)
    sapply(parts, function(p1) {
      if (length(p1) == 3) {
        return(sprintf("%s%03d", p1[2], as.integer(p1[3])))
      } else {
        return(x)
      }
    })
  }

  ordered_regs <- expected_regs[order(natural_key(expected_regs))]
  ordered_genes <- expected_genes[order(natural_key(expected_genes))]
  ordered_cols <- c(ordered_regs, ordered_genes)

  prec <- vector("list", L)

  for (i in 1:L) {

    path <- file.path(base_path, sprintf("network-%d.csv", i))

    df <- read.csv(path, row.names = 1, check.names = FALSE)

    row_names <- rownames(df)
    col_names <- colnames(df)

    # Check row order
    if (!identical(row_names, ordered_regs)) {
      stop(sprintf(
        "%s: row order mismatch.\nExpected: %s\nFound: %s",
        path,
        paste(ordered_regs, collapse = ", "),
        paste(row_names, collapse = ", ")
      ))
    }

    # Check column order
    if (!identical(col_names, ordered_cols)) {
      stop(sprintf(
        "%s: column order mismatch.\nExpected: %s\nFound: %s",
        path,
        paste(ordered_cols, collapse = ", "),
        paste(col_names, collapse = ", ")
      ))
    }

    arr <- as.matrix(df)

    # Check shape
    if (!all(dim(arr) == c(no_regs, p))) {
      stop(sprintf(
        "%s: expected shape (no_regs, p), found (%d, %d)",
        path, nrow(arr), ncol(arr)
      ))
    }

    prec[[i]] <- arr

    cat(sprintf("Loaded %s: shape=%dx%d\n",
                path, nrow(arr), ncol(arr)))
  }

  return(prec)
}

