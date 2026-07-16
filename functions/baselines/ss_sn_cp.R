# SS-SN from Zhang, Zhu, and Shao (2026), "Change-Point Detection for
# Object-Valued Time Series", JBES, Section 2.2 (equations 3-5), embedded in
# their WBS-SN algorithm (Section 3) for multiple change-point search.
#
# Splits the sequence into head/middle/tail (fraction b each end), projects
# each middle point onto the head-mean-minus-tail-mean direction to get a
# univariate contrast Z_t, then runs a self-normalized CUSUM (Shao and Zhang
# 2010) on {Z_t}. Per the paper's Remark 2, Z_t can be computed from
# pairwise distances alone:
#   Zhat_t = sum_{j in tail} d(X_j, X_t) - sum_{j in head} d(X_j, X_t)
# which is exactly what .ss_sn_interval_stat computes below.
#
# This file previously lived as wbs_sn_cp.R / .wbs_sn_stat / .wbs_sn_norm /
# .wbs_sn_interval_stat / run_wbs_sn_cp, with a header comment attributing
# it to "the WBS-SN code supplied in code_object_valued_CP.zip" (the SN1/SN2
# paper's own zip). That attribution was checked term-by-term against both
# papers' equations and was wrong: the math here matches Zhang, Zhu, and
# Shao (2026)'s SS-SN exactly (see Remark 2's distance-only Zhat_t), not
# Jiang, Zhu, and Shao (2024)'s SN1/SN2 (Dn,1/Dn,2, which use Fréchet
# variance contrasts - see sn1_cp.R / sn2_cp.R for those). Only the name is
# corrected here; the statistic itself was already a faithful
# implementation of SS-SN.
#
# source("functions/baselines/wbs_common.R") first for .wbs_intervals.

.ss_sn_stat <- function(z, k) {
  n <- length(z)
  mean_left <- mean(z[seq_len(k)])
  mean_right <- mean(z[(k + 1):n])
  inter_left <- cumsum(z[seq_len(k)]) - seq_len(k) * mean_left
  inter_right <- cumsum(z[n:(k + 1)]) - seq_len(n - k) * mean_right
  denom <- sqrt(sum(inter_left^2) / n^2 + sum(inter_right^2) / n^2)
  if (!is.finite(denom) || denom <= 0) {
    return(0)
  }
  sqrt(n) * (((n - k) * k / n^2) * (mean_left - mean_right)) / denom
}

.ss_sn_norm <- function(x) sqrt(sum(abs(x)^2))

.ss_sn_interval_stat <- function(data, trim_fraction = 0.15) {
  data <- t(as.matrix(data))
  n <- ncol(data)
  trim <- floor(n * trim_fraction)
  if (trim < 1 || n - 2 * trim < 2) {
    return(c(stat = -Inf, loc = NA_real_))
  }

  head_data <- data[, seq_len(trim), drop = FALSE]
  middle_data <- data[, (trim + 1):(n - trim), drop = FALSE]
  tail_data <- data[, (n - trim + 1):n, drop = FALSE]

  z <- numeric(ncol(middle_data))
  for (t in seq_len(ncol(middle_data))) {
    head_z <- head_data - middle_data[, t]
    tail_z <- tail_data - middle_data[, t]
    z[t] <- sum(apply(tail_z, 2, .ss_sn_norm)) -
      sum(apply(head_z, 2, .ss_sn_norm))
  }

  stats <- rep(0, length(z))
  if (length(z) > 1) {
    for (k in seq_len(length(z) - 1)) {
      stats[k] <- .ss_sn_stat(z, k)
    }
  }
  c(stat = max(stats), loc = which.max(stats) + trim)
}

.ss_sn_recursive <- function(start, end, data, intervals, threshold, min_size) {
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
    values[i, ] <- .ss_sn_interval_stat(segment)
  }

  best <- which.max(values[, 1])
  if (!is.finite(values[best, 1]) || values[best, 1] <= threshold) {
    return(integer(0))
  }
  cp <- as.integer(values[best, 2] + candidates[best, 1] - 1L)
  c(
    .ss_sn_recursive(start, cp, data, intervals, threshold, min_size),
    cp,
    .ss_sn_recursive(cp + 1L, end, data, intervals, threshold, min_size)
  )
}

# Default threshold recalibrated for n=300, M=50, p=10 (alpha=0.05, iid
# Gaussian null, 1000 replicates - see simulation/wbs_threshold_calibration.R
# and simulation/wbs_threshold_calibration_n300.csv). Reconfirmed after the
# floor/round and degenerate-segment fixes in wbs_common.R (SS-SN's own
# statistic was never affected by either bug, so the value is unchanged to
# 5 decimal places from its pre-fix calibration). The original 9.756101 was
# calibrated for n=96 (see code_object_valued_CP.zip, WBS_critic_value.R)
# and corresponds to only the ~92nd percentile at n=300 (true alpha approx
# 0.08, not 0.05). Recalibrate for other n/M/p.
run_ss_sn_cp <- function(data,
                         threshold = 10.288973,
                         M = 50,
                         min_size = 20,
                         seed = 100,
                         trim_fraction = 0.15,
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
    scan_stats[i, ] <- .ss_sn_interval_stat(segment, trim_fraction = trim_fraction)
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
        .ss_sn_interval_stat(segment, trim_fraction = trim_fraction)[["stat"]]
      }, numeric(1)), na.rm = TRUE)
    })
    critical_value <- unname(stats::quantile(permuted_stats, probs = 1 - alpha, type = 1, names = FALSE))
    p_val <- (1 + sum(permuted_stats >= observed_stat)) / (1 + num_permut)
  }

  reject <- observed_stat > critical_value
  locs <- if (isTRUE(reject)) {
    sort(unique(.ss_sn_recursive(1L, n, data, intervals, critical_value, min_size)))
  } else {
    integer(0)
  }

  list(
    method = "ss_sn",
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
