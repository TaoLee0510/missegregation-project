# Phase-2 inferred-landscape experiment

For landscape `L = 1, …, 200`, phase-1 rate index `i = 1, …, 20`,
phase-2 rate index `j = 1, …, 20`, and matched replicate `r = 1, …, 10`,
the phase-2 trajectory is

`X₂(L, i, j, r) = ABM(X₁_terminal(L, i, r), f̂₁(L, i, r), p_mis[j])`.

`f̂₁(L, i, r)` is the ALFA-K landscape inferred from the corresponding
phase-1 trajectory. The design is therefore a 20 × 20 rate-pair matrix per
original landscape, with 10 matched replicate lineages in every cell:

`200 × 20 × 20 × 10 = 800,000` phase-2 trajectories.

Phase-1 ALFA-K fitting includes every karyotype present at that replicate's
endpoint, using its full recorded trajectory. The phase-2 starting population
then retains every nonzero phase-1 endpoint karyotype and proportionally
rescales their counts to 10,000 cells. Every retained state therefore has a
direct phase-1 inferred mean fitness; no Kriging fill-in or posterior samples
are used.

The C++ ABM receives the original phase-1 ALFA-K support plus the retained
endpoint states. It retains offspring only when their karyotype occurs in this
fixed map; all descendants beyond that existing ALFA-K support are removed.
The all-diploid karyotype is excluded from phase 1 and phase 2, including
initial populations and descendants. Phase-2 fitting uses every nonzero
phase-2 terminal karyotype with its full recorded trajectory, matching phase 1.
Topology analysis is three-way: original GRF truth, phase-1 inference, and
phase-2 inference. Full-support node, edge, component, and local-maximum
readouts are reported for each inferred map. The original GRF is evaluated on
the phase-1 support, phase-2 support, and their union; the full bounded
22-dimensional GRF is not enumerated. Shared-node fitness fidelity is reported
separately, alongside support overlap, retention, and phase-specific nodes.
