#!/bin/bash
# Set PROJECT_DIR and R_LIBS_USER before submitting.
#SBATCH --job-name=bounded-grf-abm
#SBATCH --array=1-200%20
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=16G

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"

cd "$PROJECT_DIR"
Rscript 03_bounded_grf_abm/03_run_bounded_grf_abm.R \
  "$PROJECT_DIR" "${SLURM_CPUS_PER_TASK:-1}" "$SLURM_ARRAY_TASK_ID"
