#!/bin/bash
# One array element runs one phase-2 ABM simulation:
# one landscape x one phase-1 rate x one inherited replicate x one phase-2 rate.
# Submit eight chunks with TASK_OFFSET=0,100000,...,700000 because the cluster
# MaxArraySize is 100001.
#SBATCH --job-name=phase2-inferred-abm
#SBATCH --array=1-100000
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
#SBATCH --qos=large

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"

task_id=$(( ${TASK_OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
if (( task_id < 1 || task_id > 800000 )); then
  echo "Phase-2 ABM task id out of range: $task_id" >&2
  exit 1
fi
zero_based=$(( task_id - 1 ))
p2_index=$(( zero_based % 20 + 1 ))
replicate_id=$(( (zero_based / 20) % 10 + 1 ))
phase1_index=$(( (zero_based / 200) % 20 + 1 ))
landscape_index=$(( zero_based / 4000 + 1 ))

cd "$PROJECT_DIR"
Rscript 05_phase2_inferred_abm/05_run_phase2_inferred_abm.R \
  "$PROJECT_DIR" "${SLURM_CPUS_PER_TASK:-1}" "$landscape_index" "$phase1_index" "$replicate_id" "$p2_index"
