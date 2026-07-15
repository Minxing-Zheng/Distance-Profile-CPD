setwd("/home1/mzheng97/dist_profile_cpd/github_repo/Distance-Profile-CPD")
source("functions/run_benchmark_methods.R")
source_run_distCPD()
source("simulation/common/generators.R")

n <- 100
mc <- 25
num_permut <- 79
schemes <- c("iid", "block", "circular")
results <- data.frame()

for (rep in seq_len(mc)) {
  cfg <- list(
    generator = "dependent_change", n = n, tau = 1/3, p = 10,
    change_kind = "mean", dependence = "arma11", phi = 0.3, theta = 0.3,
    null = TRUE, signal = 0
  )
  gen <- generate_simulation_data(cfg, seed = 70000000 + rep)
  D <- as.matrix(dist(gen$data))

  for (sch in schemes) {
    res <- run_distCPD(D, c = 0.1, num_permut = num_permut, variants = "dist_cpd",
                        ndSup = 200, seed = 80000000 + rep, permutation_scheme = sch,
                        progress = FALSE)
    results <- rbind(results, data.frame(rep = rep, scheme = sch, reject = res$dist_cpd$reject))
  }
  if (rep %% 5 == 0) cat("done replicate", rep, "/", mc, "\n")
}

cat("\n=== Empirical size (arma_pos null, dist_cpd) by permutation scheme ===\n")
print(aggregate(reject ~ scheme, data = results, FUN = mean))
write.csv(results, "simulation/toy_block_perm_results.csv", row.names = FALSE)
