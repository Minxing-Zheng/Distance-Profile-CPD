#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=simulation/logs/wbs_calibration_%j.out
#SBATCH --error=simulation/logs/wbs_calibration_%j.err

set -euo pipefail
module purge
module load legacy/CentOS7
module load gcc/11.3.0
module load openblas/0.3.20
module load r/4.2.3

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript simulation/wbs_threshold_calibration.R
