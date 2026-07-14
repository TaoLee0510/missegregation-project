#!/bin/bash
# Submit after phase-2 ABM succeeds. Each element fits 100 trajectories.
# Full design: 800,000 fits / 100 = 8,000 array elements.
#SBATCH --job-name=phase2-alfak-inference
#SBATCH --array=1-8000%40
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=8G

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"

fits_per_task=100
fit_start=$(( (SLURM_ARRAY_TASK_ID - 1) * fits_per_task + 1 ))
cd "$PROJECT_DIR"
Rscript 06_phase2_alfak_inference/06_fit_phase2_landscapes.R \
  "$PROJECT_DIR" 1 45 "$fit_start" "$fits_per_task"
