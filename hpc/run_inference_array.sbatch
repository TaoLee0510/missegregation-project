#!/bin/bash
# Submit after all phase-1 ABMs succeed. Each element handles FITS_PER_TASK fits.
#SBATCH --job-name=alfak-inference
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
cd "$PROJECT_DIR"
fits_per_task=${FITS_PER_TASK:-100}
manifest=${ALFAK_MANIFEST:-}
if [[ -n "$manifest" ]]; then
  fit_start=$(( (SLURM_ARRAY_TASK_ID - 1) * fits_per_task + 1 ))
  Rscript 04_alfak_inference/04_fit_alfak_landscapes.R \
    "$PROJECT_DIR" 1 45 "$fit_start" "$fits_per_task" "$manifest"
else
  fit_start=$(( (SLURM_ARRAY_TASK_ID - 1) * fits_per_task + 1 ))
  Rscript 04_alfak_inference/04_fit_alfak_landscapes.R \
    "$PROJECT_DIR" 1 45 "$fit_start" "$fits_per_task"
fi
