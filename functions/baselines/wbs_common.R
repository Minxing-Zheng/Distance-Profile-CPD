# Shared utilities for the SN1/SN2 (Jiang, Zhu, and Shao 2024, "Two-sample
# and change-point inference for non-Euclidean valued time series", EJS) and
# SS-SN (Zhang, Zhu, and Shao 2026, "Change-Point Detection for Object-Valued
# Time Series", JBES) baseline wrappers in sn1_cp.R, sn2_cp.R, and
# ss_sn_cp.R.

# Random wild-binary-segmentation intervals: M pairs (s, e) with s < e and
# e - s >= min_size, used identically by all three methods' WBS recursion.
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

# Fréchet subsample-variance contrast at split r within (a, b), for
# Euclidean data. Returns c(Tn, TnC):
#   Tn(r;a,b)  = (r-a)(b-r)/(b-a) * (Vhat[a,r] - Vhat[r,b])
#   TnC(r;a,b) = (r-a)(b-r)/(b-a) * (VhatC[r;a,b] - Vhat[a,r] - Vhat[r,b])
# matching Jiang, Zhu, and Shao (2024), Section 4.1 exactly (their Tn and
# TnC). SN1 (Dn,1) uses only Tn; SN2 (Dn,2) uses both Tn and TnC - see
# sn1_cp.R / sn2_cp.R.
.frechet_tn <- function(r, a, b, n, data) {
  na <- floor(a * n)
  nr <- floor(r * n)
  nb <- floor(b * n)
  i <- nr - na
  l <- nb - na
  # Compound floating-point rounding (b, r are often k/n recomputed from an
  # already-divided fraction) can occasionally push na/nr/nb to a degenerate,
  # zero-or-negative-width segment even though the intended integer split is
  # well-defined. Rather than change floor()'s semantics (which would break
  # the exact match to the reference implementation - see
  # simulation/verify_sn_baselines.R), return a neutral zero contribution for
  # these rare edge cases instead of crashing on the matrix reshape below.
  if (i <= 0 || (l - i) <= 0) {
    return(c(0, 0))
  }
  left <- data[(na + 1):nr, , drop = FALSE]
  right <- data[(nr + 1):nb, , drop = FALSE]
  mu_left <- colMeans(matrix(left, nrow = i))
  mu_right <- colMeans(matrix(right, nrow = l - i))
  v_left <- sum((left - rep(1, i) %*% t(mu_left))^2) / i
  v_left_cross <- sum((left - rep(1, i) %*% t(mu_right))^2) / i
  v_right <- sum((right - rep(1, l - i) %*% t(mu_right))^2) / (l - i)
  v_right_cross <- sum((right - rep(1, l - i) %*% t(mu_left))^2) / (l - i)
  (r - a) * (b - r) / (b - a) *
    c(v_left - v_right, v_left_cross - v_left + v_right_cross - v_right)
}
