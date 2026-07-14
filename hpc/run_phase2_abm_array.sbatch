#!/bin/bash
# One array element runs one phase-2 ABM simulation:
# one landscape x one phase-1 rate x one inherited replicate x one phase-2 rate.
# Submit with --array=1-${N}, chunked with TASK_OFFSET when N exceeds the
# cluster MaxArraySize. N = n_landscapes * 20 * 10 * 20.
#SBATCH --job-name=phase2-inferred-abm
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=4G
#SBATCH --qos=large

set -euo pipefail
: "${PROJECT_DIR:?Export PROJECT_DIR before submitting.}"
: "${R_LIBS_USER:?Export R_LIBS_USER before submitting.}"
: "${SLURM_ARRAY_TASK_ID:?Submit through Slurm.}"

n_landscapes=$(Rscript -e 'cat(nrow(read.csv(file.path(Sys.getenv("PROJECT_DIR"), "data", "landscapes", "manifest.csv"), stringsAsFactors = FALSE)))')
n_p_mis=$(Rscript -e 'cat(nrow(read.csv(file.path(Sys.getenv("PROJECT_DIR"), "outputs", "bounded_grf", "p_mis_lhs.csv"), stringsAsFactors = FALSE)))')
n_replicates=${N_REPLICATES:-10}
per_phase1_rate=$(( n_replicates * n_p_mis ))
per_landscape=$(( n_p_mis * per_phase1_rate ))
total_tasks=$(( n_landscapes * per_landscape ))
task_id=$(( ${TASK_OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
if (( task_id < 1 || task_id > total_tasks )); then
  echo "Phase-2 ABM task id out of range: $task_id (valid 1-${total_tasks})" >&2
  exit 1
fi
zero_based=$(( task_id - 1 ))
p2_index=$(( zero_based % n_p_mis + 1 ))
replicate_id=$(( (zero_based / n_p_mis) % n_replicates + 1 ))
phase1_index=$(( (zero_based / per_phase1_rate) % n_p_mis + 1 ))
landscape_index=$(( zero_based / per_landscape + 1 ))

cd "$PROJECT_DIR"
Rscript 05_phase2_inferred_abm/05_run_phase2_inferred_abm.R \
  "$PROJECT_DIR" "${SLURM_CPUS_PER_TASK:-1}" "$landscape_index" "$phase1_index" "$replicate_id" "$p2_index"
