# SN2 (D_n,2) from Jiang, Zhu, and Shao (2024), "Two-sample and change-point
# inference for non-Euclidean valued time series", EJS, Section 4.1.
#
# D_n,2(k) = n * { [Tn(k/n;0,1)]^2 + [TnC(k/n;0,1)]^2 } / (Ln(k) + Rn(k))
#
# Uses both the Fréchet-variance contrast Tn and the "contaminated variance"
# cross term TnC, so unlike SN1 it is sensitive to both mean and variance
# changes. This file previously lived as hdcp_wbs.R / .hdcp_tn /
# .hdcp_interval_stat / run_hdcp_wbs; the math is unchanged, only the name
# is corrected to match its actual source (this was verified term-by-term
# against the paper's Tn/TnC/Dn,2 equations, not against the "MultCpt"/hdcp
# code the original file's comment pointed to - that comment was misleading
# about provenance, not about the math itself, which is correct).
#
# source("functions/baselines/wbs_common.R") first for .wbs_intervals and
# .frechet_tn.

.sn2_interval_stat <- function(data, grid_fraction = 0.15, trim_fraction = 0.05) {
  data <- as.matrix(data)
  n <- nrow(data)
  grid <- floor(n * grid_fraction)
  trim <- floor(n * trim_fraction)
  if (grid < 1 || trim < 1 || (grid + 1) > (n - grid)) {
    return(c(stat = -Inf, loc = NA_real_))
  }

  stats <- sapply((grid + 1):(n - grid), function(k) {
    ds <- .frechet_tn(k / n, 0, 1, n, data)
    d_num <- ds[1]^2
    s_num <- ds[2]^2
    left_grid <- ((trim + 1):(k - trim)) / n
    right_grid <- ((k + trim + 1):(n - trim)) / n
    if (length(left_grid) == 0 || length(right_grid) == 0) {
      return(NA_real_)
    }
    left <- sapply(left_grid, .frechet_tn, 0, k / n, n, data)
    right <- sapply(right_grid, .frechet_tn, k / n, 1, n, data)
    d_denom <- sum(left[1, ]^2) + sum(right[1, ]^2)
    s_denom <- sum(left[2, ]^2) + sum(right[2, ]^2)
    denom <- d_denom + s_denom
    if (!is.finite(denom) || denom <= 0) {
      return(NA_real_)
    }
    (n * d_num + n * s_num) / denom
  })

  c(stat = max(stats, na.rm = TRUE), loc = which.max(stats) + grid)
}

.sn2_recursive <- function(start, end, data, intervals, threshold, min_size,
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
    values[i, ] <- .sn2_interval_stat(segment, grid_fraction, trim_fraction)
  }

  best <- which.max(values[, 1])
  if (!is.finite(values[best, 1]) || values[best, 1] <= threshold) {
    return(integer(0))
  }
  cp <- as.integer(values[best, 2] + candidates[best, 1] - 1L)
  c(
    .sn2_recursive(start, cp, data, intervals, threshold, min_size, grid_fraction, trim_fraction),
    cp,
    .sn2_recursive(cp + 1L, end, data, intervals, threshold, min_size, grid_fraction, trim_fraction)
  )
}

# Default threshold recalibrated for n=300, M=50, p=10 (alpha=0.05, iid
# Gaussian null, 1000 replicates - see simulation/wbs_threshold_calibration.R
# and simulation/wbs_threshold_calibration_n300.csv, computed under this
# file's former name hdcp_wbs.R; the math is identical so the value carries
# over unchanged). Recalibrate for other n/M/p.
run_sn2_cp <- function(data,
                       threshold = 26.49393,
                       M = 50,
                       min_size = 20,
                       seed = 1,
                       grid_fraction = 0.15,
                       trim_fraction = 0.05,
                       alpha = 0.05,
                       num_permut = 0) {
  data <- as.matrix(data)
  n <- nrow(data)
  intervals <- .wbs_intervals(n, M = M, min_size = min_size, seed = seed)
  if (is.null(intervals) || nrow(intervals) == 0) {
    stop("No valid WBS intervals generated; increase n or decrease min_size.")
  }

  scan_stats <- matrix(NA_real_, nrow(intervals), 2)
  for (i in seq_len(nrow(intervals))) {
    segment <- data[intervals[i, 1]:intervals[i, 2], , drop = FALSE]
    scan_stats[i, ] <- .sn2_interval_stat(segment, grid_fraction, trim_fraction)
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
        .sn2_interval_stat(segment, grid_fraction, trim_fraction)[["stat"]]
      }, numeric(1)), na.rm = TRUE)
    })
    critical_value <- unname(stats::quantile(permuted_stats, probs = 1 - alpha, type = 1, names = FALSE))
    p_val <- (1 + sum(permuted_stats >= observed_stat)) / (1 + num_permut)
  }

  reject <- observed_stat > critical_value
  locs <- if (isTRUE(reject)) {
    sort(unique(.sn2_recursive(
      1L, n, data, intervals, critical_value, min_size,
      grid_fraction, trim_fraction
    )))
  } else {
    integer(0)
  }

  list(
    method = "sn2",
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
