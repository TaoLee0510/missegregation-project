#!/bin/bash
# SLURM submission script for 02_run_alfak.R (Slim Version - Corrected R call)

# --- Configuration ---
PROJECT_ROOT_DIR=$(pwd)
LOG_SUBDIR="logs"
R_SCRIPT_RELATIVE_PATH="scripts/02_run_alfak.R" # Path to your R script from project root
INPUT_FILES_RELATIVE_DIR="data/processed/salehi/alfak_inputs" # Path to input files from project root

LOG_DIR="${PROJECT_ROOT_DIR}/${LOG_SUBDIR}"
R_SCRIPT_FULL_PATH="${PROJECT_ROOT_DIR}/${R_SCRIPT_RELATIVE_PATH}"
INPUT_FILES_FULL_DIR="${PROJECT_ROOT_DIR}/${INPUT_FILES_RELATIVE_DIR}"

mkdir -p "${LOG_DIR}"

declare -a MINOBS_VALUES_ARRAY=(5 10 20) # Define your minobs values
INPUT_FILES_BASENAME_ARRAY=()

if [ ! -d "${INPUT_FILES_FULL_DIR}" ] || [ ! -r "${INPUT_FILES_FULL_DIR}" ]; then
    echo "Error: Input directory ${INPUT_FILES_FULL_DIR} not found or not readable. Exiting."
    exit 1
fi

while IFS= read -r -d $'\0' file_path; do
    INPUT_FILES_BASENAME_ARRAY+=("$(basename "$file_path")")
done < <(find "${INPUT_FILES_FULL_DIR}" -maxdepth 1 -name "*.Rds" -print0)

if [ ${#INPUT_FILES_BASENAME_ARRAY[@]} -eq 0 ]; then
    echo "Error: No .Rds files found in ${INPUT_FILES_FULL_DIR}. Exiting."
    exit 1
fi

declare -a ALL_TASK_MINOBS_ARGS=()
declare -a ALL_TASK_INPUT_FILE_BASENAMES=()

for minobs_val in "${MINOBS_VALUES_ARRAY[@]}"; do
    for file_basename in "${INPUT_FILES_BASENAME_ARRAY[@]}"; do
        ALL_TASK_MINOBS_ARGS+=("${minobs_val}")
        ALL_TASK_INPUT_FILE_BASENAMES+=("${file_basename}")
    done
done

NUM_TOTAL_TASKS="${#ALL_TASK_MINOBS_ARGS[@]}"

if [ "$NUM_TOTAL_TASKS" -eq 0 ]; then
  echo "Error: No tasks generated. Exiting."
  exit 1
fi

SLURM_JOB_NAME_PREFIX=${0##*/} && SLURM_JOB_NAME_PREFIX=${SLURM_JOB_NAME_PREFIX%%.*}
CONCURRENCY_JOB_LIMIT="" # Optional: e.g., %10

echo "Submitting ${NUM_TOTAL_TASKS} tasks for ${SLURM_JOB_NAME_PREFIX} to SLURM."
echo "R_SCRIPT: ${R_SCRIPT_FULL_PATH}"
echo "LOGS_DIR: ${LOG_DIR}"

sbatch <<EOSUBMIT
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-06:00:00
#SBATCH --qos=small
#SBATCH --array=0-$((${NUM_TOTAL_TASKS}-1))${CONCURRENCY_JOB_LIMIT}
#SBATCH --job-name=${SLURM_JOB_NAME_PREFIX}
#SBATCH --output=${LOG_DIR}/${SLURM_JOB_NAME_PREFIX}_%A_%a.out
#SBATCH --error=${LOG_DIR}/${SLURM_JOB_NAME_PREFIX}_%A_%a.err

# --- SLURM Task ---
# Ensure your R environment is set up (e.g., module load R/version) if Rscript is not in default PATH
ml R/4.2.0

TASK_MINOBS_LIST_EXPANDED=(${ALL_TASK_MINOBS_ARGS[*]})
TASK_FILES_LIST_EXPANDED=(${ALL_TASK_INPUT_FILE_BASENAMES[*]})

CURRENT_MINOBS_ARG="\${TASK_MINOBS_LIST_EXPANDED[\$SLURM_ARRAY_TASK_ID]}"
CURRENT_INPUT_FILE_BASENAME="\${TASK_FILES_LIST_EXPANDED[\$SLURM_ARRAY_TASK_ID]}"

R_SCRIPT_TO_EXECUTE="${R_SCRIPT_FULL_PATH}"
INPUT_FILES_PARENT_DIR="${INPUT_FILES_FULL_DIR}"
FULL_PATH_TO_CURRENT_INPUT_FILE="\${INPUT_FILES_PARENT_DIR}/\${CURRENT_INPUT_FILE_BASENAME}"

echo "SLURM Array Task ID: \${SLURM_ARRAY_TASK_ID} | MinObs: \${CURRENT_MINOBS_ARG} | File: \${CURRENT_INPUT_FILE_BASENAME}"

# R_COMMAND to use positional arguments
R_COMMAND="Rscript \${R_SCRIPT_TO_EXECUTE} \${CURRENT_MINOBS_ARG} \"\${FULL_PATH_TO_CURRENT_INPUT_FILE}\""
eval "\${R_COMMAND}"
EXIT_STATUS=\$?

if [ \$EXIT_STATUS -ne 0 ]; then
    echo "Error in SLURM Array Task ID: \${SLURM_ARRAY_TASK_ID} | MinObs: \${CURRENT_MINOBS_ARG} | File: \${CURRENT_INPUT_FILE_BASENAME} | Exit Status: \${EXIT_STATUS}"
fi
exit \${EXIT_STATUS}
EOSUBMIT

echo "SLURM job array submitted. Check status with squeue -u \$USER and logs in ${LOG_DIR}"


