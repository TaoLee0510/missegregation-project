# Phase-2 original-GRF continuation experiment

For each landscape listed in `data/landscapes/manifest.csv`, phase-1 rate index `i = 1, …, 20`,
phase-2 rate index `j = 1, …, 20`, and matched replicate `r = 1, …, 10`,
the phase-2 trajectory is

`X₂(L, i, j, r) = ABM(X₁_terminal(L, i, r), f_GRF(L), p_mis[j])`.

`f_GRF(L)` is the original bounded-GRF landscape used by the corresponding
phase-1 ABM. The design is therefore a 20 × 20 rate-pair matrix per original
landscape, with 10 matched replicate lineages in every cell:

`n_landscapes × 20 × 20 × 10` phase-2 trajectories.

Phase-1 ALFA-K fitting remains part of the analysis branch: it estimates
`f̂₁(L, i, r)` from the corresponding phase-1 trajectory, but `f̂₁` is not
used as the phase-2 ABM fitness source. The phase-2 starting population
retains every nonzero phase-1 endpoint karyotype and proportionally rescales
their counts to 10,000 cells.

The C++ ABM receives the original bounded-GRF centroids, lambda, and
chromosome-specific PDX bounds. It evaluates every in-bounds descendant
directly on the GRF; descendants outside the bounded domain are removed. The
all-diploid karyotype is excluded from phase 1 and phase 2, including initial
populations and descendants. Phase-2 fitting uses the recorded phase-2
trajectory, matching phase 1.

Topology analysis is three-way: original GRF truth, phase-1 inference, and
phase-2 inference. Full-support node, edge, component, and local-maximum
readouts are reported for each inferred map. The original GRF is evaluated on
the phase-1 support, phase-2 support, and their union; the full bounded
22-dimensional GRF is not enumerated. Shared-node fitness fidelity is reported
separately, alongside support overlap, retention, and phase-specific nodes.
