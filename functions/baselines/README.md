# Baseline helpers

This folder contains sourceable helper functions for benchmark methods used in
the simulation comparisons.

## Existing wrappers

- `energy_cp.R`, `ecp_distmat_input.R`, `energyChangePoint.cpp`: energy-CP
  baseline for a distance matrix.
- `graph_cp.R`: graph-CP baseline via `gSeg`.
- `kernel_mmd.R`: kernel/MMD scan baseline.

## Self-normalized WBS baselines (`wbs_common.R`, `sn1_cp.R`, `sn2_cp.R`, `ss_sn_cp.R`)

Three related but distinct methods, sharing `wbs_common.R`'s random-interval
generator and (for SN1/SN2) Fréchet-variance-contrast helper:

- `sn1_cp.R` (method name `sn1`): SN1 / $D_{n,1}$ from Jiang, Zhu, and Shao
  (2024), "Two-sample and change-point inference for non-Euclidean valued
  time series", EJS, Section 4.1. Uses only the Fréchet-variance contrast
  $T_n$ - sensitive to variance changes only, not pure mean shifts.
- `sn2_cp.R` (method name `sn2`): SN2 / $D_{n,2}$, same paper. Adds the
  "contaminated variance" cross term $T_n^C$ on top of $T_n$ - sensitive to
  both mean and variance changes.
- `ss_sn_cp.R` (method name `ss_sn`): SS-SN from Zhang, Zhu, and Shao (2026),
  "Change-Point Detection for Object-Valued Time Series", JBES, Section 2.2
  (their own separate paper and method, built on a head/tail-projected,
  distance-only self-normalized CUSUM - see the file header for the exact
  correspondence to their equations 3-5 and Remark 2).

These three were reconstructed and verified term-by-term against both
papers' published equations (not just against the supplied code archives,
`code_object_valued_CP.zip` and `code.zip`, whose internal file/variable
naming turned out to be misleading about which paper each statistic
actually came from - `wbs_sn_cp.R`/`hdcp_wbs.R`, the files these replaced,
had the SS-SN/SN2 correspondence backwards relative to their names).

Default thresholds (`sn1`: n=300-calibrated but unverified against a paper
table; `sn2`/`ss_sn`: n=300-calibrated, see
`simulation/wbs_threshold_calibration.R` and
`simulation/wbs_threshold_calibration_n300.csv`) assume `n=300, M=50, p=10`
and should be recalibrated for other settings.

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
  methods = c("sn1", "sn2", "ss_sn"),
  num_permut = 0,
  progress = FALSE
)

res$sn1$candidate_loc
res$sn2$candidate_loc
res$ss_sn$candidate_loc
```
