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
Rscript 02_landscape_generation/02_generate_landscapes.R "$PROJECT_DIR"
```

Submit the ABM as one independent array element per landscape. Each element
runs 200 rate/replicate combinations (20 rates × 10 replicates) using the
CPUs allocated to it.

```sh
sbatch hpc/run_abm_array.sbatch
```

Wait for every ABM element to succeed. A failed simulation or fit returns a
non-zero exit code, so Slurm records the element as failed. Per-landscape
status files are written under `outputs/bounded_grf/`.

Afterward, count the observations and submit 100 fits per array element. The
full phase-1 experiment has 40,000 fits, therefore 400 inference jobs.

```sh
n=$(find outputs/bounded_grf -name abm_observations.rds | wc -l)
tasks=$(( (n + 99) / 100 ))
test "$n" -eq 40000 || { echo "Expected 40000 phase-1 observations, found $n"; exit 1; }
test "$tasks" -eq 400 && sbatch hpc/run_inference_array.sbatch
```

Update the `#SBATCH` time and memory directives after a one-landscape ABM
smoke run and one-fit smoke run on the target cluster. The templates avoid
shared output writes: each array element owns a landscape or inference result.

## Phase 2

For every original landscape, the second phase forms a 20 × 20 grid of
`p_mis_phase1` and `p_mis_phase2`. The 10 phase-1 replicate lineages are
continued, yielding 200 × 20 × 20 × 10 = 800,000 phase-2 trajectories. Each
phase-2 array element owns one landscape and one phase-1 rate, so it runs 200
trajectories (10 replicates × 20 phase-2 rates).

```sh
sbatch hpc/run_phase2_abm_array.sbatch
```

Phase-2 inference is batched at 100 fits per element to avoid scheduler array
limits. After the phase-2 ABM completes, count observations and submit the
corresponding number of batches.

```sh
n=$(find outputs/phase2_abm -name abm_observations.rds | wc -l)
tasks=$(( (n + 99) / 100 ))
test "$n" -eq 800000 || { echo "Expected 800000 phase-2 observations, found $n"; exit 1; }
test "$tasks" -eq 8000 && sbatch hpc/run_phase2_inference_array.sbatch
```

After all 800,000 phase-2 fits succeed, compute topology in 200 bounded tasks,
then build the figures from their compact CSV shards.

```sh
sbatch hpc/run_phase2_topology_array.sbatch
Rscript 08_phase2_visualization/08_make_phase2_figures.R "$PROJECT_DIR"
```
