#!/bin/bash
# Submit after phase-2 ABM succeeds. Each element fits 100 trajectories.
#SBATCH --job-name=phase2-alfak-inference
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
#SBATCH --qos=large

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"

if [[ -n "${R_MODULE:-}" ]]; then
  module load "$R_MODULE"
fi
fits_per_task=${FITS_PER_TASK:-100}
cd "$PROJECT_DIR"
manifest=${ALFAK_MANIFEST:-}
if [[ -n "$manifest" ]]; then
  fit_start=$(( (SLURM_ARRAY_TASK_ID - 1) * fits_per_task + 1 ))
  Rscript 06_phase2_alfak_inference/06_fit_phase2_landscapes.R \
    "$PROJECT_DIR" 1 45 "$fit_start" "$fits_per_task" "$manifest"
else
  fit_start=$(( (SLURM_ARRAY_TASK_ID - 1) * fits_per_task + 1 ))
  Rscript 06_phase2_alfak_inference/06_fit_phase2_landscapes.R \
    "$PROJECT_DIR" 1 45 "$fit_start" "$fits_per_task"
fi
