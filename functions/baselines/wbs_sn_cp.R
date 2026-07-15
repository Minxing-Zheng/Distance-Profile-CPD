# Wild binary segmentation with a self-normalized statistic.
#
# This is a cleaned, sourceable wrapper around the WBS-SN code supplied in
# code_object_valued_CP.zip. It accepts an n by p matrix with observations in
# rows and returns estimated change-point locations.

.wbs_sn_stat <- function(z, k) {
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

.wbs_sn_norm <- function(x) sqrt(sum(abs(x)^2))

.wbs_sn_interval_stat <- function(data, trim_fraction = 0.15) {
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
    z[t] <- sum(apply(tail_z, 2, .wbs_sn_norm)) -
      sum(apply(head_z, 2, .wbs_sn_norm))
  }

  stats <- rep(0, length(z))
  if (length(z) > 1) {
    for (k in seq_len(length(z) - 1)) {
      stats[k] <- .wbs_sn_stat(z, k)
    }
  }
  c(stat = max(stats), loc = which.max(stats) + trim)
}

.wbs_intervals <- function(n, M = 50, min_size = 20, seed = 100) {
  if (n < min_size) {
    return(NULL)
  }
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
  on.exit({
    if (is.null(old_seed)) {
      rm(.Random.seed, envir = .GlobalEnv)
    } else {
      .Random.seed <<- old_seed
    }
  }, add = TRUE)

  set.seed(seed)
  intervals <- matrix(0L, M, 2)
  counter <- 1L
  attempts <- 0L
  while (counter <= M && attempts < 1000L * M) {
    attempts <- attempts + 1L
    s <- sort(sample.int(n, 2, replace = TRUE))
    if (s[2] - s[1] >= min_size) {
      intervals[counter, ] <- s
      counter <- counter + 1L
    }
  }
  intervals[seq_len(counter - 1L), , drop = FALSE]
}

.wbs_sn_recursive <- function(start, end, data, intervals, threshold, min_size) {
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
    values[i, ] <- .wbs_sn_interval_stat(segment)
  }

  best <- which.max(values[, 1])
  if (!is.finite(values[best, 1]) || values[best, 1] <= threshold) {
    return(integer(0))
  }
  cp <- as.integer(values[best, 2] + candidates[best, 1] - 1L)
  c(
    .wbs_sn_recursive(start, cp, data, intervals, threshold, min_size),
    cp,
    .wbs_sn_recursive(cp + 1L, end, data, intervals, threshold, min_size)
  )
}

run_wbs_sn_cp <- function(data,
                          threshold = 9.756101,
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
    scan_stats[i, ] <- .wbs_sn_interval_stat(segment, trim_fraction = trim_fraction)
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
        .wbs_sn_interval_stat(segment, trim_fraction = trim_fraction)[["stat"]]
      }, numeric(1)), na.rm = TRUE)
    })
    critical_value <- unname(stats::quantile(permuted_stats, probs = 1 - alpha, type = 1, names = FALSE))
    p_val <- (1 + sum(permuted_stats >= observed_stat)) / (1 + num_permut)
  }

  reject <- observed_stat > critical_value
  locs <- if (isTRUE(reject)) {
    sort(unique(.wbs_sn_recursive(1L, n, data, intervals, critical_value, min_size)))
  } else {
    integer(0)
  }

  list(
    method = "wbs_sn_cp",
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
