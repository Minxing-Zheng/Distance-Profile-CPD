# High-dimensional/object-valued WBS baseline.
#
# This is a cleaned wrapper around the WBS code supplied in code.zip. It uses
# the self-normalized D2 statistic from the supplied MultCpt/wbs scripts.

.hdcp_tn <- function(r, a, b, n, data) {
  l <- round((b - a) * n, 0)
  i <- round((r - a) * n, 0)
  left <- data[(round(a * n, 0) + 1):(r * n), , drop = FALSE]
  right <- data[(round(r * n, 0) + 1):(b * n), , drop = FALSE]
  mu_left <- colMeans(matrix(left, nrow = i))
  mu_right <- colMeans(matrix(right, nrow = l - i))
  v_left <- sum((left - rep(1, i) %*% t(mu_left))^2) / i
  v_left_cross <- sum((left - rep(1, i) %*% t(mu_right))^2) / i
  v_right <- sum((right - rep(1, l - i) %*% t(mu_right))^2) / (l - i)
  v_right_cross <- sum((right - rep(1, l - i) %*% t(mu_left))^2) / (l - i)
  (r - a) * (b - r) / (b - a) *
    c(v_left - v_right, v_left_cross - v_left + v_right_cross - v_right)
}

.hdcp_interval_stat <- function(data, grid_fraction = 0.15, trim_fraction = 0.05) {
  data <- as.matrix(data)
  n <- nrow(data)
  grid <- round(n * grid_fraction, 0)
  trim <- round(n * trim_fraction, 0)
  if (grid < 1 || trim < 1 || (grid + 1) > (n - grid)) {
    return(c(stat = -Inf, loc = NA_real_))
  }

  stats <- sapply((grid + 1):(n - grid), function(k) {
    ds <- .hdcp_tn(k / n, 0, 1, n, data)
    d_num <- ds[1]^2
    s_num <- ds[2]^2
    left_grid <- ((trim + 1):(k - trim)) / n
    right_grid <- ((k + trim + 1):(n - trim)) / n
    if (length(left_grid) == 0 || length(right_grid) == 0) {
      return(NA_real_)
    }
    left <- sapply(left_grid, .hdcp_tn, 0, k / n, n, data)
    right <- sapply(right_grid, .hdcp_tn, k / n, 1, n, data)
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

.hdcp_recursive <- function(start, end, data, intervals, threshold, min_size,
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
    values[i, ] <- .hdcp_interval_stat(segment, grid_fraction, trim_fraction)
  }

  best <- which.max(values[, 1])
  if (!is.finite(values[best, 1]) || values[best, 1] <= threshold) {
    return(integer(0))
  }
  cp <- as.integer(values[best, 2] + candidates[best, 1] - 1L)
  c(
    .hdcp_recursive(start, cp, data, intervals, threshold, min_size, grid_fraction, trim_fraction),
    cp,
    .hdcp_recursive(cp + 1L, end, data, intervals, threshold, min_size, grid_fraction, trim_fraction)
  )
}

run_hdcp_wbs <- function(data,
                         threshold = 145.4887,
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
    scan_stats[i, ] <- .hdcp_interval_stat(segment, grid_fraction, trim_fraction)
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
        .hdcp_interval_stat(segment, grid_fraction, trim_fraction)[["stat"]]
      }, numeric(1)), na.rm = TRUE)
    })
    critical_value <- unname(stats::quantile(permuted_stats, probs = 1 - alpha, type = 1, names = FALSE))
    p_val <- (1 + sum(permuted_stats >= observed_stat)) / (1 + num_permut)
  }

  reject <- observed_stat > critical_value
  locs <- if (isTRUE(reject)) {
    sort(unique(.hdcp_recursive(
      1L, n, data, intervals, critical_value, min_size,
      grid_fraction, trim_fraction
    )))
  } else {
    integer(0)
  }

  list(
    method = "hdcp_wbs",
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
