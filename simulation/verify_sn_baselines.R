setwd("/home1/mzheng97/dist_profile_cpd/github_repo/Distance-Profile-CPD")

## ---- Reference implementations, extracted verbatim from the original
## authors' own code (functions/baselines/code_object_valued_CP.zip,
## "github update/Section_4_2/DGP3_size_and_power.R") ----

ref_compute_Tn <- function(r, a, b, n, X, X.cumsum) {
  nr <- floor(r * n)
  na <- floor(a * n)
  nb <- floor(b * n)
  p <- ncol(X)
  X0 <- matrix(X[(na + 1):nr, ], nr - na, p)
  X1 <- matrix(X[(nr + 1):nb, ], nb - nr, p)
  if (na == 0) {
    mu0 <- X.cumsum[nr, ] / nr
  } else {
    mu0 <- (X.cumsum[nr, ] - X.cumsum[na, ]) / (nr - na)
  }
  mu1 <- (X.cumsum[nb, ] - X.cumsum[nr, ]) / (nb - nr)
  V0 <- sum((X0 - cbind(rep(1, nr - na)) %*% rbind(mu0))^2) / (nr - na)
  V0c <- sum((X0 - cbind(rep(1, nr - na)) %*% rbind(mu1))^2) / (nr - na)
  V1 <- sum((X1 - cbind(rep(1, nb - nr)) %*% rbind(mu1))^2) / (nb - nr)
  V1c <- sum((X1 - cbind(rep(1, nb - nr)) %*% rbind(mu0))^2) / (nb - nr)
  (r - a) * (b - r) / (b - a) * c((V0 - V1), (V0c - V0 + V1c - V1))
}

ref_Jiang_test <- function(X, eta1 = 0.15, eta2 = 0.05) {
  res_Jiang <- matrix(0, nrow = 1, ncol = 10)
  X.cumsum <- apply(X, 2, cumsum)
  n <- nrow(X)
  grid <- floor(n * eta1)
  trim <- floor(n * eta2)
  Test_res <- NULL
  for (k in (grid + 1):(n - grid)) {
    DS_res <- ref_compute_Tn(k / n, 0, 1, n, X, X.cumsum)
    Dn_num <- (DS_res[1])^2
    Sn_num <- (DS_res[2])^2
    Left <- sapply(((trim + 1):(k - trim)) / n, ref_compute_Tn, 0, k / n, n, X, X.cumsum)
    Right <- sapply(((k + trim + 1):(n - trim)) / n, ref_compute_Tn, k / n, 1, n, X, X.cumsum)
    Dn_denorm <- sum(Left[1, ]^2) + sum(Right[1, ]^2)
    Sn_denorm <- sum(Left[2, ]^2) + sum(Right[2, ]^2)
    Tn <- c(n * Dn_num / Dn_denorm, n * (Dn_num + Sn_num) / (Dn_denorm + Sn_denorm))
    Test_res <- cbind(Test_res, Tn)
  }
  res_Jiang[1:2] <- apply(Test_res, 1, max)
  res_Jiang
}

ref_SN_test <- function(ts, k) {
  n <- length(ts)
  mean1 <- mean(ts[1:k])
  mean2 <- mean(ts[(k + 1):n])
  inter1 <- cumsum(ts[1:k]) - (1:k) * mean1
  inter2 <- cumsum(ts[n:(k + 1)]) - (1:(n - k)) * mean2
  M1 <- sum(inter1^2) / n^2
  M2 <- sum(inter2^2) / n^2
  sqrt(n) * (((n - k) * k / n^2) * (mean1 - mean2)) / sqrt(M1 + M2)
}

ref_sum_squared_vals <- function(x) (sum(abs(x)^2))^(1 / 2)

ref_ss_sn_stat <- function(seqi, b = 0.15) {
  # seqi: p x n (matches the original script's `t(seqi)` layout)
  n <- ncol(seqi)
  seqihead <- seqi[, 1:floor(n * b)]
  seqimiddle <- seqi[, (floor(n * b) + 1):(n - floor(n * b))]
  seqitail <- seqi[, (n - floor(n * b) + 1):n]
  Zt <- rep(0, ncol(seqimiddle))
  for (t in seq_len(ncol(seqimiddle))) {
    headZ <- seqihead - seqimiddle[, t]
    tailZ <- seqitail - seqimiddle[, t]
    Zt[t] <- sum(apply(tailZ, 2, ref_sum_squared_vals)) - sum(apply(headZ, 2, ref_sum_squared_vals))
  }
  TNpl <- rep(0, ncol(seqimiddle))
  for (k in seq_len(ncol(seqimiddle) - 1)) TNpl[k] <- ref_SN_test(Zt, k)
  max(TNpl)
}

## ---- Our implementation ----
source("functions/baselines/wbs_common.R")
source("functions/baselines/sn1_cp.R")
source("functions/baselines/sn2_cp.R")
source("functions/baselines/ss_sn_cp.R")

## ---- Compare on several independent random datasets ----
set.seed(7)
n <- 100
p <- 3
n_checks <- 5
cat(sprintf("%-8s %14s %14s %14s\n", "check", "SN1(diff)", "SN2(diff)", "SS-SN(diff)"))
max_abs_diff <- c(sn1 = 0, sn2 = 0, ss_sn = 0)
for (i in seq_len(n_checks)) {
  X <- matrix(rnorm(n * p), n, p)
  if (i > 1) {
    # inject an obvious mean+variance change on later checks
    X[(n / 2 + 1):n, ] <- X[(n / 2 + 1):n, ] * 1.8 + 1.5
  }

  ref_res <- ref_Jiang_test(X, eta1 = 0.15, eta2 = 0.05)
  ref_sn1 <- ref_res[1, 1]
  ref_sn2 <- ref_res[1, 2]
  ref_ss <- ref_ss_sn_stat(t(X), b = 0.15)

  our_sn1 <- .sn1_interval_stat(X, grid_fraction = 0.15, trim_fraction = 0.05)[["stat"]]
  our_sn2 <- .sn2_interval_stat(X, grid_fraction = 0.15, trim_fraction = 0.05)[["stat"]]
  our_ss <- .ss_sn_interval_stat(X, trim_fraction = 0.15)[["stat"]]

  d1 <- abs(our_sn1 - ref_sn1)
  d2 <- abs(our_sn2 - ref_sn2)
  d3 <- abs(our_ss - ref_ss)
  max_abs_diff["sn1"] <- max(max_abs_diff["sn1"], d1)
  max_abs_diff["sn2"] <- max(max_abs_diff["sn2"], d2)
  max_abs_diff["ss_sn"] <- max(max_abs_diff["ss_sn"], d3)
  cat(sprintf("%-8d %14.10f %14.10f %14.10f\n", i, d1, d2, d3))
}

cat("\nMax absolute difference across all checks (should be ~0, floating point only):\n")
print(max_abs_diff)
cat(ifelse(all(max_abs_diff < 1e-8), "PASS: all three match the original reference code to floating-point precision.\n",
           "FAIL: at least one statistic disagrees with the reference implementation.\n"))
