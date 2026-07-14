#!/bin/bash
# Submit after all 40,000 phase-1 ABMs succeed. Each element handles 100 fits.
#SBATCH --job-name=alfak-inference
#SBATCH --array=1-400%40
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=8G

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"

cd "$PROJECT_DIR"
fits_per_task=100
fit_start=$(( (SLURM_ARRAY_TASK_ID - 1) * fits_per_task + 1 ))
Rscript 04_alfak_inference/04_fit_alfak_landscapes.R \
  "$PROJECT_DIR" 1 45 "$fit_start" "$fits_per_task"
