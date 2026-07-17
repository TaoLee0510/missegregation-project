# HPC execution (Slurm)

Load an R module with C++17 support, then use a persistent library location;
do not use node-local temporary storage for `R_LIBS_USER`.

```sh
export PROJECT_DIR=/path/to/Missegregation_Project
export R_LIBS_USER=/path/to/persistent/Rlib
cd "$PROJECT_DIR"
Rscript hpc/install_r_environment.R "$PROJECT_DIR"
```

The setup step installs all CRAN dependencies and recompiles `alfakR` from
source for the cluster's Linux architecture. It requires network access to
CRAN (or set `R_REPOSITORY` to the institution's CRAN mirror) and a C++17
compiler. Run it on a login or build node, not in every array task.

Generate the two shared prerequisites once:

```sh
Rscript 01_reference_ploidy/01_derive_reference_ploidy.R "$PROJECT_DIR"
Rscript 02_landscape_generation/02_generate_landscapes.R "$PROJECT_DIR" 200
```

The second command accepts the landscape count as its second argument. It
also writes `data/landscapes/grf_generation_parameters.csv`, which records
the exact GRF lambda, seeds, fixed settings, file digests, and centroid
coordinates for every generated landscape.

Submit the phase-1 ABM as one independent array element per ABM simulation:
one landscape, one phase-1 missegregation rate, and one replicate.

```sh
n_landscapes=$(Rscript -e 'cat(nrow(read.csv("data/landscapes/manifest.csv", stringsAsFactors = FALSE)))')
phase1_tasks=$(( n_landscapes * 20 * 10 ))
sbatch --array=1-${phase1_tasks} hpc/run_abm_array.sbatch
```

Wait for every ABM element to succeed. A failed simulation or fit returns a
non-zero exit code, so Slurm records the element as failed. Per-landscape
status files are written under `outputs/bounded_grf/`.

Afterward, submit phase-1 ALFA-K fits with the retry driver. These fits
estimate the phase-1 inferred landscapes for downstream topology analysis, but
they are not the fitness source for phase-2 ABM. The driver builds a pending
manifest before every submission, so Slurm arrays contain only fits without a
complete, provenance-matching `fit_metadata.rds` and `landscape.Rds`.

```sh
MEMORY_TIERS="32G 64G 128G 256G 512G" FITS_PER_TASK=1 \
  bash hpc/submit_alfak_retry_by_memory.sh phase1
```

Update the `#SBATCH` time and memory directives after a one-landscape ABM
smoke run and one-fit smoke run on the target cluster. The templates avoid
shared output writes: each array element owns a landscape or inference result.

## Phase 2

For every original landscape, the second phase forms a 20 × 20 grid of
`p_mis_phase1` and `p_mis_phase2`. The 10 phase-1 replicate lineages are
continued on the same original bounded-GRF landscape, yielding
`n_landscapes × 20 × 20 × 10` phase-2 trajectories. Each phase-2 array element
owns one ABM simulation: one landscape, one phase-1 rate, one inherited
replicate, and one phase-2 rate. Phase-2 ABM only requires the completed
phase-1 ABM endpoints and the original GRF inputs; it does not depend on
phase-1 ALFA-K. If the resulting task count exceeds the cluster MaxArraySize,
submit chunks with `TASK_OFFSET`.

```sh
n_rates=$(Rscript -e 'cat(nrow(read.csv("outputs/bounded_grf/p_mis_lhs.csv", stringsAsFactors = FALSE)))')
phase2_tasks=$(( n_landscapes * n_rates * 10 * n_rates ))
max_array_tasks=100000
offset=0
while [ "$offset" -lt "$phase2_tasks" ]; do
  chunk=$(( phase2_tasks - offset ))
  if [ "$chunk" -gt "$max_array_tasks" ]; then chunk=$max_array_tasks; fi
  sbatch --array=1-${chunk} --export=ALL,TASK_OFFSET=${offset} hpc/run_phase2_abm_array.sbatch
  offset=$(( offset + chunk ))
done
```

Phase-2 inference uses the same pending-manifest retry driver. It retries only
the incomplete phase-2 fits at each memory tier.

```sh
MEMORY_TIERS="32G 64G 128G 256G 512G" FITS_PER_TASK=1 \
  bash hpc/submit_alfak_retry_by_memory.sh phase2
```

After all phase-2 fits succeed and the phase-1 ALFA-K branch is complete,
compute topology in one task per landscape, then build the figures from their
compact CSV shards.

```sh
sbatch --array=1-${n_landscapes} hpc/run_phase2_topology_array.sbatch
Rscript 08_phase2_visualization/08_make_phase2_figures.R "$PROJECT_DIR"
```

To run both ALFA-K phases and then submit topology/figures automatically:

```sh
MEMORY_TIERS="32G 64G 128G 256G 512G" FITS_PER_TASK=1 \
  bash hpc/submit_alfak_retry_by_memory.sh both
```
