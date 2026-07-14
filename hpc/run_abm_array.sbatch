#!/bin/bash
# One array element runs one phase-1 ABM simulation:
# one landscape x one missegregation rate x one replicate.
# Submit with --array=1-${N}, where N = n_landscapes * 20 * 10.
#SBATCH --job-name=bounded-grf-abm
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
#SBATCH --qos=large

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"

n_landscapes=$(Rscript -e 'cat(nrow(read.csv(file.path(Sys.getenv("PROJECT_DIR"), "data", "landscapes", "manifest.csv"), stringsAsFactors = FALSE)))')
n_p_mis=${N_P_MIS:-20}
n_replicates=${N_REPLICATES:-10}
per_landscape=$(( n_p_mis * n_replicates ))
total_tasks=$(( n_landscapes * per_landscape ))
task_id=$(( ${TASK_OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
if (( task_id < 1 || task_id > total_tasks )); then
  echo "Phase-1 ABM task id out of range: $task_id (valid 1-${total_tasks})" >&2
  exit 1
fi
zero_based=$(( task_id - 1 ))
landscape_index=$(( zero_based / per_landscape + 1 ))
within_landscape=$(( zero_based % per_landscape ))
p_index=$(( within_landscape / n_replicates + 1 ))
replicate_id=$(( within_landscape % n_replicates + 1 ))

cd "$PROJECT_DIR"
Rscript 03_bounded_grf_abm/03_run_bounded_grf_abm.R \
  "$PROJECT_DIR" "${SLURM_CPUS_PER_TASK:-1}" "$landscape_index" "$p_index" "$replicate_id"
