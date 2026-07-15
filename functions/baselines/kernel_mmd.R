.mmd_scan_stat <- function(distmat, c = 0.1) {
  distmat <- .run_all_validate_distmat(distmat)
  n <- nrow(distmat)
  cp_min <- ceiling(n * c)
  cp_max <- n - ceiling(n * c)
  if (cp_min < 2) cp_min <- 2
  if (cp_max > n - 1) cp_max <- n - 1

  off_diag <- distmat[upper.tri(distmat)]
  bw <- stats::median(off_diag[off_diag > 0], na.rm = TRUE)
  if (!is.finite(bw) || bw <= 0) {
    bw <- 1
  }
  gram_matrix <- exp(-(distmat^2) / (2 * bw^2))

  obj_vals <- rep(NA_real_, n)
  for (t in cp_min:cp_max) {
    left <- seq_len(t)
    right <- (t + 1):n
    term1 <- length(right) / (length(left) * n) * sum(gram_matrix[left, left, drop = FALSE])
    term2 <- -2 / n * sum(gram_matrix[left, right, drop = FALSE])
    term3 <- length(left) / (length(right) * n) * sum(gram_matrix[right, right, drop = FALSE])
    obj_vals[t] <- term1 + term2 + term3
  }

  loc <- which.max(obj_vals)
  c(stat = unname(obj_vals[loc]), loc = loc)
}

run_kernel_mmd <- function(distmat,
                           c = 0.1,
                           num_permut = 500,
                           alpha = 0.05,
                           seed = NULL,
                           progress = FALSE) {
  distmat <- .run_all_validate_distmat(distmat)
  n <- nrow(distmat)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  runtime <- system.time({
    stats <- matrix(NA_real_, nrow = num_permut + 1, ncol = 2)
    colnames(stats) <- c("stat", "loc")
    for (b in 0:num_permut) {
      if (isTRUE(progress)) {
        message("kernel permutation ", b, "/", num_permut)
      }
      current_distmat <- if (b == 0) {
        distmat
      } else {
        ind <- sample.int(n)
        distmat[ind, ind, drop = FALSE]
      }
      stats[b + 1, ] <- .mmd_scan_stat(current_distmat, c = c)
    }
  })[["elapsed"]]

  test_stat <- unname(stats[1, "stat"])
  permuted_test_stat <- if (num_permut > 0) unname(stats[-1, "stat"]) else numeric(0)
  p_val <- if (length(permuted_test_stat) == 0) {
    NA_real_
  } else {
    (1 + sum(permuted_test_stat >= test_stat)) / (1 + length(permuted_test_stat))
  }
  critical_value <- if (length(permuted_test_stat) == 0) {
    NA_real_
  } else {
    unname(stats::quantile(permuted_test_stat, probs = 1 - alpha, type = 1, names = FALSE))
  }
  reject <- if (is.na(p_val)) NA else p_val <= alpha
  candidate_loc <- unname(stats[1, "loc"])
  loc <- if (isTRUE(reject)) candidate_loc else NA_integer_

  list(
    method = "kernel_mmd",
    test_stat = test_stat,
    permuted_test_stat = permuted_test_stat,
    critical_value = critical_value,
    p_val = p_val,
    reject = reject,
    loc = loc,
    candidate_loc = candidate_loc,
    observed_stat = test_stat,
    permuted_stats = permuted_test_stat,
    all_stats = unname(stats[, "stat"]),
    detected = reject,
    c = c,
    num_permut = num_permut,
    alpha = alpha,
    runtime_sec = runtime
  )
}
