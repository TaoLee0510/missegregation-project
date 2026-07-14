#!/bin/bash
# One bounded topology-analysis task per original landscape.
#SBATCH --job-name=phase2-topology
#SBATCH --array=1-200%20
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"
cd "$PROJECT_DIR"
Rscript 07_phase2_topology_analysis/07_compare_phase2_topologies.R "$PROJECT_DIR" "$SLURM_ARRAY_TASK_ID"
