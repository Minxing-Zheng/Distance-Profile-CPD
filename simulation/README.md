# Clean simulation pipeline

`run_simulations.R` is the production entry point for new single-change-point
experiments. The older scenario scripts remain available for provenance, but
new cluster runs should use this configuration-driven pipeline.

## Design

- One configuration row defines one scenario.
- `mc` gives the number of Monte Carlo replicates for that scenario.
- One Slurm array index maps to one scenario/replicate pair.
- Every replicate is checkpointed independently, so rerunning an array skips
  completed results.
- A small tidy CSV is written for analysis and an optional RDS preserves raw
  permutation statistics and method output.
- A location error is recorded only when the method rejects. This makes
  `conditional_mae` genuinely conditional on detection.

Supported generators are:

- `gaussian_mean`
- `gaussian_scale`
- `gaussian_mixture`
- `gaussian_tail`
- `gaussian_covariance`
- `ar_mean`
- `distributional_location`
- `distributional_scale`
- `network_pa`

Supported methods are `dist_profile`, `energy`, `graph`, `kernel`, `sn1`,
`sn2`, and `ss_sn`. Separate multiple-change-point experiments should use a
dedicated runner because their output and performance metrics differ.

## Local smoke test

From the repository root:

```bash
Rscript simulation/run_simulations.R \
  --config simulation/configs/smoke.csv \
  --output_dir simulation/results/smoke \
  --array_index 1

Rscript simulation/run_simulations.R \
  --config simulation/configs/smoke.csv \
  --output_dir simulation/results/smoke \
  --array_index 2
```

Run all generator smoke tests with indices 1 through 9 using
`all_generators_smoke.csv`.

## Cluster pilot

`revision_pilot.csv` contains 20 null and 20 alternative replicates, so it has
40 array tasks:

```bash
mkdir -p simulation/logs
sbatch --array=1-40 simulation/slurm_array.sh \
  simulation/configs/revision_pilot.csv \
  simulation/results/revision_pilot
```

The scripts assume the same GCC, OpenBLAS, and R modules used by the existing
cluster jobs. Adjust module versions or the Slurm partition at submission time
if needed.

## Summaries

```bash
Rscript simulation/summarize_results.R \
  --input_dir simulation/results/revision_pilot \
  --output simulation/results/revision_pilot_summary.csv
```

The summary reports empirical rejection rate, conditional mean and median
absolute error, candidate-location MAE, failures, and runtime. Under a null
scenario, rejection rate is empirical size and localization metrics are left
undefined.

## Important configuration columns

- `scenario_id`: unique output directory name.
- `generator`: one of the supported generators above.
- `mc`: number of Monte Carlo replicates.
- `n`, `tau`, `p`, `signal`: scenario parameters.
- `null`: if true, generate no change and save `true_cp = NA`.
- `c`: scan trimming fraction.
- `num_permut`, `alpha`: calibration settings.
- `methods`: semicolon-separated method names.
- `variants`: semicolon-separated distance-profile variants.
- `ndsup`: integration-grid size for `dist_cpd_uniform`.
- `base_seed`: replicate `r` uses `base_seed + r - 1` for data generation.

Generator-specific optional columns include `covariance`, `rho`, `grid_size`,
`object_size`, `base_variance`, `spike_count`, `active_fraction`, and the
outlier/transform columns documented in `common/generators.R`.
