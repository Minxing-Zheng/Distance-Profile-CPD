#!/bin/bash
# Submit with, for example:
#   sbatch --array=1-40 simulation/slurm_array.sh simulation/configs/revision_pilot.csv simulation/results/revision_pilot

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=simulation/logs/%A_%a.out
#SBATCH --error=simulation/logs/%A_%a.err

set -euo pipefail

CONFIG=${1:?Pass the configuration CSV as argument 1}
OUTPUT_DIR=${2:-simulation/results}

module purge
module load legacy/CentOS7
module load gcc/11.3.0
module load openblas/0.3.20
module load r/4.2.3

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p simulation/logs "$OUTPUT_DIR"
Rscript simulation/run_simulations.R \
  --config "$CONFIG" \
  --output_dir "$OUTPUT_DIR" \
  --array_index "$SLURM_ARRAY_TASK_ID" \
  --save_raw true \
  --overwrite false
