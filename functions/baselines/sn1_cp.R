# SN1 (D_n,1) from Jiang, Zhu, and Shao (2024), "Two-sample and change-point
# inference for non-Euclidean valued time series", EJS, Section 4.1.
#
# D_n,1(k) = n * [Tn(k/n;0,1)]^2 /
#            ( sum_l [Tn(l/n;0,k/n)]^2 + sum_l [Tn(l/n;k/n,1)]^2 )
#
# Uses only the plain Fréchet-variance contrast Tn (not the "contaminated
# variance" cross term TnC that SN2 additionally uses) - the paper notes
# this means SN1 targets changes in Fréchet variance only, not changes that
# are purely in the mean. See sn2_cp.R for the mean-and-variance-sensitive
# sibling statistic, and ss_sn_cp.R for the separate, mean-shift-focused
# SS-SN statistic from Zhang, Zhu, and Shao (2026).
#
# source("functions/baselines/wbs_common.R") first for .wbs_intervals and
# .frechet_tn.

.sn1_interval_stat <- function(data, grid_fraction = 0.15, trim_fraction = 0.05) {
  data <- as.matrix(data)
  n <- nrow(data)
  grid <- floor(n * grid_fraction)
  trim <- floor(n * trim_fraction)
  if (grid < 1 || trim < 1 || (grid + 1) > (n - grid)) {
    return(c(stat = -Inf, loc = NA_real_))
  }

  stats <- sapply((grid + 1):(n - grid), function(k) {
    d_num <- .frechet_tn(k / n, 0, 1, n, data)[1]^2
    left_grid <- ((trim + 1):(k - trim)) / n
    right_grid <- ((k + trim + 1):(n - trim)) / n
    if (length(left_grid) == 0 || length(right_grid) == 0) {
      return(NA_real_)
    }
    left <- vapply(left_grid, function(l) .frechet_tn(l, 0, k / n, n, data)[1], numeric(1))
    right <- vapply(right_grid, function(r) .frechet_tn(r, k / n, 1, n, data)[1], numeric(1))
    denom <- sum(left^2) + sum(right^2)
    if (!is.finite(denom) || denom <= 0) {
      return(NA_real_)
    }
    n * d_num / denom
  })

  c(stat = max(stats, na.rm = TRUE), loc = which.max(stats) + grid)
}

.sn1_recursive <- function(start, end, data, intervals, threshold, min_size,
                           grid_fraction, trim_fraction) {
  if ((end - start) < min_size) {
    return(integer(0))
  }
  candidates <- intervals[intervals[, 1] >= start & intervals[, 2] <= end, , drop = FALSE]
  if (nrow(candidates) == 0) {
    return(integer(0))
  }

  values <- matrix(NA_real_, nrow(candidates), 2)
  for (i in seq_len(nrow(candidates))) {
    segment <- data[candidates[i, 1]:candidates[i, 2], , drop = FALSE]
    values[i, ] <- .sn1_interval_stat(segment, grid_fraction, trim_fraction)
  }

  best <- which.max(values[, 1])
  if (!is.finite(values[best, 1]) || values[best, 1] <= threshold) {
    return(integer(0))
  }
  cp <- as.integer(values[best, 2] + candidates[best, 1] - 1L)
  c(
    .sn1_recursive(start, cp, data, intervals, threshold, min_size, grid_fraction, trim_fraction),
    cp,
    .sn1_recursive(cp + 1L, end, data, intervals, threshold, min_size, grid_fraction, trim_fraction)
  )
}

# Default threshold calibrated for n=300, M=50, p=10 (alpha=0.05, iid
# Gaussian null, 1000 replicates - see simulation/wbs_threshold_calibration.R
# and simulation/wbs_threshold_calibration_n300.csv). Recalibrate for other
# n/M/p.
run_sn1_cp <- function(data,
                       threshold = NULL,
                       M = 50,
                       min_size = 20,
                       seed = 2,
                       grid_fraction = 0.15,
                       trim_fraction = 0.05,
                       alpha = 0.05,
                       num_permut = 0) {
  if (is.null(threshold)) threshold <- 15.19  # placeholder until calibration run completes
  data <- as.matrix(data)
  n <- nrow(data)
  intervals <- .wbs_intervals(n, M = M, min_size = min_size, seed = seed)
  if (is.null(intervals) || nrow(intervals) == 0) {
    stop("No valid WBS intervals generated; increase n or decrease min_size.")
  }

  scan_stats <- matrix(NA_real_, nrow(intervals), 2)
  for (i in seq_len(nrow(intervals))) {
    segment <- data[intervals[i, 1]:intervals[i, 2], , drop = FALSE]
    scan_stats[i, ] <- .sn1_interval_stat(segment, grid_fraction, trim_fraction)
  }
  observed_stat <- max(scan_stats[, 1], na.rm = TRUE)
  candidate_interval <- which.max(scan_stats[, 1])
  candidate_loc <- as.integer(scan_stats[candidate_interval, 2] + intervals[candidate_interval, 1] - 1L)

  permuted_stats <- numeric(0)
  critical_value <- threshold
  p_val <- NA_real_
  if (num_permut > 0) {
    set.seed(seed)
    permuted_stats <- replicate(num_permut, {
      ind <- sample.int(n)
      perm_data <- data[ind, , drop = FALSE]
      max(vapply(seq_len(nrow(intervals)), function(i) {
        segment <- perm_data[intervals[i, 1]:intervals[i, 2], , drop = FALSE]
        .sn1_interval_stat(segment, grid_fraction, trim_fraction)[["stat"]]
      }, numeric(1)), na.rm = TRUE)
    })
    critical_value <- unname(stats::quantile(permuted_stats, probs = 1 - alpha, type = 1, names = FALSE))
    p_val <- (1 + sum(permuted_stats >= observed_stat)) / (1 + num_permut)
  }

  reject <- observed_stat > critical_value
  locs <- if (isTRUE(reject)) {
    sort(unique(.sn1_recursive(
      1L, n, data, intervals, critical_value, min_size,
      grid_fraction, trim_fraction
    )))
  } else {
    integer(0)
  }

  list(
    method = "sn1",
    test_stat = observed_stat,
    permuted_test_stat = permuted_stats,
    critical_value = critical_value,
    p_val = p_val,
    reject = reject,
    loc = if (length(locs) > 0) locs[1] else NA_integer_,
    candidate_loc = candidate_loc,
    all_locs = locs,
    intervals = intervals,
    alpha = alpha,
    num_permut = num_permut
  )
}
