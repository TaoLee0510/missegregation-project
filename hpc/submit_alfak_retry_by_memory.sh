#!/bin/bash
# Submit ALFAK inference in memory tiers, retrying only unfinished fits.
#
# Usage:
#   PROJECT_DIR=/path/to/repo R_LIBS_USER=/path/to/Rlib \
#     bash hpc/submit_alfak_retry_by_memory.sh [phase1|phase2|both]

set -euo pipefail

phase_mode=${1:-both}
case "$phase_mode" in
  phase1|phase2|both) ;;
  *) echo "Usage: $0 [phase1|phase2|both]" >&2; exit 2 ;;
esac

PROJECT_DIR=${PROJECT_DIR:-$(pwd)}
RESULTS_DIR=${RESULTS_DIR:-"$PROJECT_DIR/results"}
LOG_DIR=${LOG_DIR:-"$RESULTS_DIR/logs"}
MANIFEST_DIR=${MANIFEST_DIR:-"$RESULTS_DIR/manifests/alfak_retry_$(date +%Y%m%d_%H%M%S)"}
R_LIBS_USER=${R_LIBS_USER:?Export R_LIBS_USER before submitting.}
ALFAKR_SOURCE=${ALFAKR_SOURCE:-"$PROJECT_DIR/packages/alfakR"}
ALFAKR_COMPILE=${ALFAKR_COMPILE:-false}
R_MODULE=${R_MODULE:-}
NBOOT=${NBOOT:-45}
FITS_PER_TASK=${FITS_PER_TASK:-1}
QOS=${QOS:-large}
TIME_LIMIT=${TIME_LIMIT:-1-00:00:00}
MEMORY_TIERS=${MEMORY_TIERS:-"32G 64G 128G 256G 512G"}
POLL_SECONDS=${POLL_SECONDS:-60}
TOPOLOGY_MEM=${TOPOLOGY_MEM:-4G}
FIGURES_MEM=${FIGURES_MEM:-4G}
SUBMIT_DOWNSTREAM=${SUBMIT_DOWNSTREAM:-true}

mkdir -p "$LOG_DIR" "$MANIFEST_DIR"

if [[ "$FITS_PER_TASK" -lt 1 ]]; then
  echo "FITS_PER_TASK must be a positive integer." >&2
  exit 2
fi

run_rscript() {
  if [[ -n "$R_MODULE" ]]; then
    module load "$R_MODULE"
  fi
  Rscript "$@"
}

manifest_path_for() {
  local phase=$1
  local label=$2
  echo "$MANIFEST_DIR/${phase}_pending_${label}.csv"
}

pending_count() {
  local manifest=$1
  awk 'NR > 1 { n++ } END { print n + 0 }' "$manifest"
}

array_count_for() {
  local pending=$1
  echo $(( (pending + FITS_PER_TASK - 1) / FITS_PER_TASK ))
}

build_manifest() {
  local phase=$1
  local label=$2
  local manifest
  manifest=$(manifest_path_for "$phase" "$label")
  echo "Building $phase pending manifest: $manifest"
  run_rscript "$PROJECT_DIR/hpc/build_alfak_pending_manifest.R" "$PROJECT_DIR" "$phase" "$manifest" "$NBOOT"
}

phase_script() {
  case "$1" in
    phase1) echo "$PROJECT_DIR/hpc/run_inference_array.sbatch" ;;
    phase2) echo "$PROJECT_DIR/hpc/run_phase2_inference_array.sbatch" ;;
  esac
}

phase_job_prefix() {
  case "$1" in
    phase1) echo "mis-p1-alfak" ;;
    phase2) echo "mis-p2-alfak" ;;
  esac
}

submit_phase_job() {
  local phase=$1
  local mem=$2
  local manifest=$3
  local pending=$4
  local array_count
  local mem_label
  local job_name
  array_count=$(array_count_for "$pending")
  mem_label=${mem//[^[:alnum:]]/}
  job_name="$(phase_job_prefix "$phase")-${mem_label}"
  sbatch --parsable \
    --job-name="$job_name" \
    --qos="$QOS" \
    --time="$TIME_LIMIT" \
    --cpus-per-task=1 \
    --mem="$mem" \
    --array=1-"$array_count" \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    --export=ALL,PROJECT_DIR="$PROJECT_DIR",R_LIBS_USER="$R_LIBS_USER",ALFAKR_SOURCE="$ALFAKR_SOURCE",ALFAKR_COMPILE="$ALFAKR_COMPILE",R_MODULE="$R_MODULE",FITS_PER_TASK="$FITS_PER_TASK",ALFAK_MANIFEST="$manifest" \
    "$(phase_script "$phase")"
}

wait_for_job() {
  local job_id=$1
  echo "Waiting for Slurm job $job_id"
  while squeue -h -j "$job_id" | grep -q .; do
    sleep "$POLL_SECONDS"
  done
  sleep 15
}

print_sacct_counts() {
  local job_id=$1
  echo "Slurm accounting for $job_id:"
  sacct -j "$job_id" -X -n -P --format=State,ExitCode |
    awk -F'|' 'NF >= 2 { key = $1 "|" $2; count[key]++ } END { for (key in count) print count[key] "|" key }' |
    sort || true
}

run_phase_with_memory_tiers() {
  local phase=$1
  local attempt=0
  local manifest
  local pending
  local mem
  local job_id

  for mem in $MEMORY_TIERS; do
    attempt=$((attempt + 1))
    build_manifest "$phase" "before_${attempt}_${mem}"
    manifest=$(manifest_path_for "$phase" "before_${attempt}_${mem}")
    pending=$(pending_count "$manifest")
    echo "$phase pending before ${mem}: $pending"
    if [[ "$pending" -eq 0 ]]; then
      echo "$phase is already complete."
      return 0
    fi

    job_id=$(submit_phase_job "$phase" "$mem" "$manifest" "$pending")
    echo "$phase submitted at $mem: $job_id"
    wait_for_job "$job_id"
    print_sacct_counts "$job_id"

    build_manifest "$phase" "after_${attempt}_${mem}"
    manifest=$(manifest_path_for "$phase" "after_${attempt}_${mem}")
    pending=$(pending_count "$manifest")
    echo "$phase pending after ${mem}: $pending"
    if [[ "$pending" -eq 0 ]]; then
      echo "$phase completed at $mem."
      return 0
    fi
  done

  echo "$phase still has $pending pending fit(s) after memory tiers: $MEMORY_TIERS" >&2
  echo "Last pending manifest: $manifest" >&2
  return 1
}

submit_downstream_jobs() {
  if [[ "$SUBMIT_DOWNSTREAM" != "true" ]]; then
    echo "SUBMIT_DOWNSTREAM=$SUBMIT_DOWNSTREAM; not submitting topology/figures."
    return 0
  fi

  local n_landscapes
  local topo_job
  local fig_job
  n_landscapes=$(run_rscript -e 'args <- commandArgs(trailingOnly = TRUE); cat(nrow(read.csv(file.path(args[[1]], "data", "landscapes", "manifest.csv"), stringsAsFactors = FALSE)))' "$PROJECT_DIR")
  if [[ "$n_landscapes" -lt 1 ]]; then
    echo "Could not determine landscape count." >&2
    return 1
  fi

  topo_job=$(sbatch --parsable \
    --job-name=mis-p2-topology \
    --qos="$QOS" \
    --time="$TIME_LIMIT" \
    --cpus-per-task=1 \
    --mem="$TOPOLOGY_MEM" \
    --array=1-"$n_landscapes" \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    --export=ALL,PROJECT_DIR="$PROJECT_DIR",R_LIBS_USER="$R_LIBS_USER",R_MODULE="$R_MODULE" \
    "$PROJECT_DIR/hpc/run_phase2_topology_array.sbatch")
  echo "Submitted topology: $topo_job"

  fig_cmd='if [ -n "${R_MODULE:-}" ]; then module load "$R_MODULE"; fi; cd "$PROJECT_DIR"; Rscript 08_phase2_visualization/08_make_phase2_figures.R "$PROJECT_DIR"'
  fig_job=$(sbatch --parsable \
    --job-name=mis-p2-figures \
    --qos="$QOS" \
    --time="$TIME_LIMIT" \
    --cpus-per-task=1 \
    --mem="$FIGURES_MEM" \
    --dependency=afterok:"$topo_job" \
    --output="$LOG_DIR/%x_%j.out" \
    --error="$LOG_DIR/%x_%j.err" \
    --export=ALL,PROJECT_DIR="$PROJECT_DIR",R_LIBS_USER="$R_LIBS_USER",R_MODULE="$R_MODULE" \
    --wrap="$fig_cmd")
  echo "Submitted figures: $fig_job"
}

case "$phase_mode" in
  phase1)
    run_phase_with_memory_tiers phase1
    ;;
  phase2)
    run_phase_with_memory_tiers phase2
    ;;
  both)
    run_phase_with_memory_tiers phase1
    run_phase_with_memory_tiers phase2
    submit_downstream_jobs
    ;;
esac
