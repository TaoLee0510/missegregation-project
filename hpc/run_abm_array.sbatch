#!/bin/bash
# One array element runs one phase-1 ABM simulation:
# one landscape x one missegregation rate x one replicate.
#SBATCH --job-name=bounded-grf-abm
#SBATCH --array=1-40000
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
#SBATCH --qos=large

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"

task_id=$(( ${TASK_OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
if (( task_id < 1 || task_id > 40000 )); then
  echo "Phase-1 ABM task id out of range: $task_id" >&2
  exit 1
fi
zero_based=$(( task_id - 1 ))
landscape_index=$(( zero_based / 200 + 1 ))
within_landscape=$(( zero_based % 200 ))
p_index=$(( within_landscape / 10 + 1 ))
replicate_id=$(( within_landscape % 10 + 1 ))

cd "$PROJECT_DIR"
Rscript 03_bounded_grf_abm/03_run_bounded_grf_abm.R \
  "$PROJECT_DIR" "${SLURM_CPUS_PER_TASK:-1}" "$landscape_index" "$p_index" "$replicate_id"
