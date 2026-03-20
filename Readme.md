# Distance-Profile-CPD: File Guide

This short guide explains what each major file does.

## Repository layout

- `functions/`: Core methods (R + C++) used by simulations and real-data scripts.
- `simulation/`: Synthetic experiments for different shift settings.
- `real_data/`: Real-data pipelines and analysis scripts.

## `functions/` (core methods)

### Main R and C++ functions
- `functions/depth_CPDcpp.cpp`: C++ implementation of the change-point scan based on a single test statistic, used for the `dist-CP-U` method.
- `functions/depth_CPDcpp_ALL.cpp`: C++ implementation of the change-point scan that returns multiple test-statistic variants, used for `dist-CP-F`, `dist-CP-AD`, and `dist-CP`.

### Baseline wrappers / helpers
- `functions/ecp_distmat_input.R`: R wrapper for ECP-style methods with distance matrix input.
- `functions/kcp_distmat_input.R`: R wrapper for kernel CP methods with distance/kernel matrix input.

## `real_data/` 

### Electricity
- `real_data/electricity/elec_preprocess.R`: Preprocessing pipeline for raw electricity dataset.
- `real_data/electricity/electricity_compositional.R`: CPD script for electricity dataset.

### MIT Reality Mining
- `real_data/mit_reality_mining/reality_mining_1392.RData`: MIT Reality Mining data object from https://github.com/cran/GreedySBTM/blob/master/data/reality_mining_1392.RData
- `real_data/mit_reality_mining/mit_reality_mining.R`: CPD script for MIT Reality Mining data.

### Human vs Machine Text
- `real_data/mix_human_machine_text/mix_text_detection.R`: CPD script for mixed human/machine text dataset.
