setwd("/home1/mzheng97/dist_profile_cpd/github_repo/Distance-Profile-CPD")
source("functions/baselines/wbs_common.R")
source("functions/baselines/sn1_cp.R")
source("functions/baselines/sn2_cp.R")
source("functions/baselines/ss_sn_cp.R")

n <- 300
p <- 10
M <- 50
min_size <- 20
reps <- 1000
quantile_levels <- c(0.90, 0.95, 0.975, 0.99)

intervals <- .wbs_intervals(n, M = M, min_size = min_size, seed = 100)
cat("n_blocks intervals:", nrow(intervals), "\n")

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("using", n_cores, "cores\n")

t0 <- Sys.time()
one_rep <- function(k) {
  set.seed(20300000 + k)
  data <- matrix(rnorm(n * p), n, p)
  stat_sn1 <- max(vapply(seq_len(nrow(intervals)), function(i) {
    seg <- data[intervals[i, 1]:intervals[i, 2], , drop = FALSE]
    .sn1_interval_stat(seg)[["stat"]]
  }, numeric(1)), na.rm = TRUE)
  stat_sn2 <- max(vapply(seq_len(nrow(intervals)), function(i) {
    seg <- data[intervals[i, 1]:intervals[i, 2], , drop = FALSE]
    .sn2_interval_stat(seg)[["stat"]]
  }, numeric(1)), na.rm = TRUE)
  stat_ss_sn <- max(vapply(seq_len(nrow(intervals)), function(i) {
    seg <- data[intervals[i, 1]:intervals[i, 2], , drop = FALSE]
    .ss_sn_interval_stat(seg)[["stat"]]
  }, numeric(1)), na.rm = TRUE)
  c(sn1 = stat_sn1, sn2 = stat_sn2, ss_sn = stat_ss_sn)
}

results <- parallel::mclapply(seq_len(reps), one_rep, mc.cores = n_cores)
results <- do.call(rbind, results)
cat("done, elapsed:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

out_rows <- list()
for (m in c("sn1", "sn2", "ss_sn")) {
  q <- quantile(results[, m], quantile_levels)
  cat("\n===", m, "max-statistic quantiles at n=300, M=50, p=10 (iid Gaussian null) ===\n")
  print(q)
  out_rows[[m]] <- data.frame(
    method = m, quantile = quantile_levels, alpha = 1 - quantile_levels,
    threshold_n300 = as.numeric(q)
  )
}
out <- do.call(rbind, out_rows)
write.csv(out, "simulation/wbs_threshold_calibration_n300.csv", row.names = FALSE)
cat("\nWrote simulation/wbs_threshold_calibration_n300.csv\n")
