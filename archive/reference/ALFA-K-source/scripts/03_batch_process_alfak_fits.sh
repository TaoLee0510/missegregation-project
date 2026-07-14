#! /bin/bash
#submit to Slurm
sbatch <<EOSUBMIT
#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --job-name=03_process_alfak_fits
#SBATCH --cpus-per-task=51
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1-6
#SBATCH --qos=small
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err


# Load R module and check if Rscript is available
ml R/4.2.0
which Rscript  # Verify Rscript is available

Rscript scripts/03_process_alfak_fits.R

EOSUBMIT


