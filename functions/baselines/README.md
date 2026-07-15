# Baseline helpers

This folder contains sourceable helper functions for benchmark methods used in
the simulation comparisons.

## Existing wrappers

- `energy_cp.R`, `ecp_distmat_input.R`, `energyChangePoint.cpp`: energy-CP
  baseline for a distance matrix.
- `graph_cp.R`: graph-CP baseline via `gSeg`.
- `kernel_mmd.R`: kernel/MMD scan baseline.

## Added external WBS baselines

- `wbs_sn_cp.R`: cleaned wrapper for the WBS-SN object-valued change-point code
  supplied in `code_object_valued_CP.zip`.
- `hdcp_wbs.R`: cleaned wrapper for the high-dimensional/object-valued WBS code
  supplied in `code.zip`.

The original zip files include many simulation scripts, outputs, and macOS
metadata. Only the reusable method code was wrapped here.

Example:

```r
source("functions/run_benchmark_methods.R")

set.seed(1)
X <- rbind(
  matrix(rnorm(30 * 3), 30, 3),
  matrix(rnorm(30 * 3, mean = 1), 30, 3)
)
D <- as.matrix(dist(X))

res <- run_benchmark_methods(
  data = X,
  distmat = D,
  methods = c("wbs_sn", "hdcp_wbs"),
  num_permut = 0,
  progress = FALSE
)

res$wbs_sn$candidate_loc
res$hdcp_wbs$candidate_loc
```
