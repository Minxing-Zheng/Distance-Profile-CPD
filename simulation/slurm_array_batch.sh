#!/bin/bash
# Batched variant of slurm_array.sh: each array task claims a contiguous
# block of BATCH_SIZE global replicate indices (as defined by the config's
# mc column and run_simulations.R's map_array_index) and runs them
# concurrently across the job's cpus-per-task.
#
# Submit with, for example (100 total replicates, batches of 10 -> 10 tasks):
#   BATCH_SIZE=10 sbatch --array=1-10 simulation/slurm_array_batch.sh \
#     simulation/configs/toy_run.csv simulation/results/toy_run

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --output=simulation/logs/%A_%a.out
#SBATCH --error=simulation/logs/%A_%a.err

set -euo pipefail

CONFIG=${1:?Pass the configuration CSV as argument 1}
OUTPUT_DIR=${2:-simulation/results}
BATCH_SIZE=${BATCH_SIZE:-10}
PARALLEL_JOBS=${PARALLEL_JOBS:-${SLURM_CPUS_PER_TASK:-4}}

module purge
module load legacy/CentOS7
module load gcc/11.3.0
module load openblas/0.3.20
module load r/4.2.3

# Each replicate runs single-threaded so PARALLEL_JOBS concurrent
# replicates saturate the job's cpus-per-task without oversubscribing.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p simulation/logs "$OUTPUT_DIR"

START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))

seq "$START" "$END" | xargs -P "$PARALLEL_JOBS" -I{} \
  Rscript simulation/run_simulations.R \
    --config "$CONFIG" \
    --output_dir "$OUTPUT_DIR" \
    --array_index {} \
    --save_raw true \
    --overwrite false
