#!/bin/bash
# Each element owns one (landscape, p_mis_1) block: 10 inherited replicates x
# 20 second-stage rates = 200 phase-2 trajectories.
#SBATCH --job-name=phase2-inferred-abm
#SBATCH --array=1-4000%20
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=16G

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"

phase1_index=$(( (SLURM_ARRAY_TASK_ID - 1) % 20 + 1 ))
landscape_index=$(( (SLURM_ARRAY_TASK_ID - 1) / 20 + 1 ))
cd "$PROJECT_DIR"
Rscript 05_phase2_inferred_abm/05_run_phase2_inferred_abm.R \
  "$PROJECT_DIR" "${SLURM_CPUS_PER_TASK:-1}" "$landscape_index" "$phase1_index"
