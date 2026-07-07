# Distance-Profile Change-Point Detection

This repository contains code for distance-profile change-point detection, plus baseline wrappers for energy, graph, and kernel/MMD change-point methods.

The easiest entry points are:

- `functions/run_distCPD.R`: run the proposed distance-profile variants from a distance matrix.
- `functions/run_benchmark_methods.R`: run distance-profile variants and baseline methods from a common input.

## Quick Start

Start R from the repository root, or set the working directory to this folder:

```r
setwd("/path/to/Distance-Profile-CPD")
```

Create a simple mean-change example:

```r
set.seed(1)

n1 <- 100
n2 <- 200
p <- 10

X <- rbind(
  matrix(rnorm(n1 * p), nrow = n1),
  matrix(rnorm(n2 * p, mean = 0.5), nrow = n2)
)

D <- as.matrix(dist(X))
true_cp <- n1
```

## Run Distance-Profile Only

Use `run_distCPD()` when you already have a distance matrix.

```r
source("functions/run_distCPD.R")

res <- run_distCPD(
  distmat = D,
  c = 0.1,
  num_permut = 300,
  alpha = 0.05,
  seed = 1
)
```

The default runs all four variants:

```r
names(res)
```

```text
dist_cpd_uniform
dist_cpd
dist_cpd_AD
dist_cpd_W
```

Each variant returns the same core fields:

```r
res$dist_cpd$p_val
res$dist_cpd$reject
res$dist_cpd$loc
res$dist_cpd$candidate_loc
res$dist_cpd$test_stat
res$dist_cpd$critical_value
```

Use `loc` as the accepted change-point estimate. It is `NA` if the method does not reject. Use `candidate_loc` for the scan maximizer, even when there is no rejection.

## Run Selected Variants

Running fewer variants can save time.

```r
# Basic dF variant only
res <- run_distCPD(D, variants = "dist_cpd", num_permut = 300)

# Uniform variant only
res <- run_distCPD(D, variants = "dist_cpd_uniform", num_permut = 300)

# Non-uniform variants together
res <- run_distCPD(
  D,
  variants = c("dist_cpd", "dist_cpd_AD", "dist_cpd_W"),
  num_permut = 300
)
```

Short aliases also work:

```r
res <- run_distCPD(D, variants = c("uniform", "AD"), num_permut = 300)
```

## The `ndSup` Argument

`ndSup` controls the integration grid size for `dist_cpd_uniform` only.

```r
res <- run_distCPD(
  D,
  variants = "dist_cpd_uniform",
  ndSup = 200,
  num_permut = 300
)
```

Larger `ndSup` is more accurate for the uniform-grid integral but slower. It does not affect:

- `dist_cpd`
- `dist_cpd_AD`
- `dist_cpd_W`

For large simulations, a useful sensitivity check is:

```r
for (g in c(200, 500, 1000)) {
  res <- run_distCPD(D, variants = "dist_cpd_uniform", ndSup = g)
  print(c(ndSup = g, p_val = res$dist_cpd_uniform$p_val,
          loc = res$dist_cpd_uniform$candidate_loc))
}
```

## Run Distance-Profile and Baselines

Use `run_benchmark_methods()` to run distance-profile plus energy, graph, and kernel/MMD baselines.

```r
source("functions/run_benchmark_methods.R")

res <- run_benchmark_methods(
  data = X,
  distmat = D,
  c = 0.1,
  num_permut = 300,
  alpha = 0.05,
  min_size = 50,
  methods = c("dist_profile", "energy", "graph", "kernel"),
  dist_profile_variants = c("dist_cpd", "dist_cpd_AD"),
  dist_profile_ndSup = 200,
  seed = 1
)
```

Access results:

```r
res$dist_profile$dist_cpd
res$dist_profile$dist_cpd_AD
res$energy
res$graph
res$kernel
```

The baseline outputs use the same conventions:

```r
res$energy$p_val
res$energy$reject
res$energy$loc
res$energy$candidate_loc
```

For energy-CP, `candidate_loc` is the last proposed split location from the energy maximization step. If energy-CP does not reject, `loc` is `NA`.

## Summarize One Simulation Run

Example for power/localization bookkeeping:

```r
methods <- c(
  "dist_cpd" = res$dist_profile$dist_cpd,
  "dist_cpd_AD" = res$dist_profile$dist_cpd_AD,
  "energy" = res$energy,
  "graph" = res$graph,
  "kernel" = res$kernel
)

summary <- do.call(rbind, lapply(names(methods), function(name) {
  x <- methods[[name]]
  data.frame(
    method = name,
    reject = x$reject,
    loc = x$loc,
    candidate_loc = x$candidate_loc,
    abs_error = if (isTRUE(x$reject)) abs(x$loc - true_cp) else NA_real_,
    candidate_abs_error = abs(x$candidate_loc - true_cp),
    p_val = x$p_val,
    runtime_sec = x$runtime_sec
  )
}))

summary
```

For simulation tables, report detection and localization separately:

- rejection rate / power: `mean(reject)`
- conditional MAE: `mean(abs_error, na.rm = TRUE)`
- runtime: `mean(runtime_sec)`

For null simulations, report type-I error and runtime; localization error is not meaningful when there is no true change point.

## Batch Example

The script below runs a multivariate mean-change experiment and writes detail, summary, and RDS outputs:

```bash
Rscript simulation/mul_vec/run_mean_p10_summary.R \
  --mc 100 \
  --num_permut 300 \
  --p 10 \
  --n1 100 \
  --n2 200 \
  --delta 0.5 \
  --dist_profile_variants all \
  --dist_profile_ndSup 200 \
  --cores 8 \
  --seed 2026 \
  --out_dir simulation/mul_vec/results_mean_p10_mc100_B300
```

Parallelism is over Monte Carlo replications. A single replication does not become faster just because `--cores` is larger than one.

## Repository Layout

- `functions/`: core methods and wrappers.
- `simulation/`: synthetic experiments.
- `real_data/`: real-data analysis scripts.

Important files:

- `functions/distCPD_combined.cpp`: optimized C++ implementation for distance-profile variants.
- `functions/baselines/`: energy, graph, and kernel baseline helpers.
- `functions/data_generation/gen_data.R`: covariance-generation helper used by simulation scripts.
- `functions/run_distCPD.R`: tutorial-friendly wrapper for distance-profile variants.
- `functions/run_benchmark_methods.R`: common wrapper for distance-profile and baselines.
- `simulation/mul_vec/run_mean_p10_summary.R`: batch runner for a multivariate mean-change experiment.

## Dependencies

Core distance-profile code needs:

```r
install.packages("Rcpp")
```

Baseline and simulation scripts may also use:

```r
install.packages(c("MASS", "pracma", "ade4"))
```

The graph baseline uses `gSeg`. Install it before running graph-CP experiments.
