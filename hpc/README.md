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

Afterward, count the observations and submit 100 fits per array element.

```sh
n=$(find outputs/bounded_grf -name abm_observations.rds | wc -l)
tasks=$(( (n + 99) / 100 ))
test "$n" -eq "$phase1_tasks" || { echo "Expected ${phase1_tasks} phase-1 observations, found $n"; exit 1; }
sbatch --array=1-${tasks} hpc/run_inference_array.sbatch
```

Update the `#SBATCH` time and memory directives after a one-landscape ABM
smoke run and one-fit smoke run on the target cluster. The templates avoid
shared output writes: each array element owns a landscape or inference result.

## Phase 2

For every original landscape, the second phase forms a 20 × 20 grid of
`p_mis_phase1` and `p_mis_phase2`. The 10 phase-1 replicate lineages are
continued, yielding `n_landscapes × 20 × 20 × 10` phase-2 trajectories. Each
phase-2 array element owns one ABM simulation: one landscape, one phase-1
rate, one inherited replicate, and one phase-2 rate. If the resulting task
count exceeds the cluster MaxArraySize, submit chunks with `TASK_OFFSET`.

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

Phase-2 inference is batched at 100 fits per element to avoid scheduler array
limits. After the phase-2 ABM completes, count observations and submit the
corresponding number of batches.

```sh
n=$(find outputs/phase2_abm -name abm_observations.rds | wc -l)
tasks=$(( (n + 99) / 100 ))
test "$n" -eq "$phase2_tasks" || { echo "Expected ${phase2_tasks} phase-2 observations, found $n"; exit 1; }
sbatch --array=1-${tasks} hpc/run_phase2_inference_array.sbatch
```

After all phase-2 fits succeed, compute topology in one task per landscape,
then build the figures from their compact CSV shards.

```sh
sbatch --array=1-${n_landscapes} hpc/run_phase2_topology_array.sbatch
Rscript 08_phase2_visualization/08_make_phase2_figures.R "$PROJECT_DIR"
```
