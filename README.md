# Two-stage bounded-GRF CIN / ALFA-K validation project

This project generates bounded synthetic karyotype-evolution data calibrated
to untreated PDX ploidy, then applies a second CIN phase on each phase-1
ALFA-K-inferred fitness support.

## Execution order

1. `01_reference_ploidy/01_derive_reference_ploidy.R` derives PDX-informed
   chromosome-specific copy-number support and the reference mean ploidy.
2. `02_landscape_generation/02_generate_landscapes.R` creates 200 reproducible
   GRF landscapes.
3. `03_bounded_grf_abm/03_run_bounded_grf_abm.R` runs 20 LHS `p_mis` values ×
   10 replicates per landscape (40,000 phase-1 trajectories).
4. `04_alfak_inference/04_fit_alfak_landscapes.R` fits the phase-1 landscapes.
5. `05_phase2_inferred_abm/05_run_phase2_inferred_abm.R` continues every
   phase-1 replicate under each of the 20 second-stage rates on its own
   inferred fitness support (800,000 phase-2 trajectories).
6. `06_phase2_alfak_inference/06_fit_phase2_landscapes.R` fits phase-2
   landscapes.
7. `07_phase2_topology_analysis/07_compare_phase2_topologies.R` compares
   phase-2 fits with their phase-1 inferred source landscapes.
8. `08_phase2_visualization/08_make_phase2_figures.R` produces phase-2
   figures.

Run each script from the project root, passing the root as its first argument:

```sh
Rscript 01_reference_ploidy/01_derive_reference_ploidy.R .
Rscript 02_landscape_generation/02_generate_landscapes.R .
Rscript 03_bounded_grf_abm/03_run_bounded_grf_abm.R . <workers> [landscape_index] [n_steps]
Rscript 05_phase2_inferred_abm/05_run_phase2_inferred_abm.R . <workers> <landscape_index> <p_mis_phase1_index>
```

For local use, install the local package into the R library selected by
`R_LIBS_USER` (or an existing entry in `.libPaths()`):

```sh
R CMD INSTALL --preclean -l "$R_LIBS_USER" packages/alfakR
```

For Slurm execution, dependency installation, task partitioning, restart
behavior, and submission templates, follow [hpc/README.md](hpc/README.md).

## Layout

- `data/`: Salehi reference inputs, derived reference ploidy, and generated GRF inputs.
- `outputs/bounded_grf/`: destination for the HPC ABM screen.
- `outputs/archive/`: prior exploratory and smoke-test outputs, retained but
  excluded from the current workflow.
- `packages/alfakR/`: local package compiled on the target HPC.
- `docs/`: paper and workflow notes.
- `archive/`: upstream ALFA-K reference repository and deferred future-work protocol; neither is executable.
