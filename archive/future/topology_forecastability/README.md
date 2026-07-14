# Reserved: topology forecastability branch

This directory is deliberately non-executable. It preserves only the proposed
protocol and input-data note for future work; it is not part of the current
pipeline and contains no runnable scripts.

## Scientific question

For a fixed chromosome-missegregation rate (`p_mis`), why do some fitness
landscapes yield reproducible evolutionary trajectories while others yield many
plausible karyotypic futures?

This is a companion analysis to the main `p_mis` screen. It does **not** select
or target individual karyotypes, simulate a second treatment, or use shortest
paths as intervention candidates.

## Primary hypothesis

At a given `p_mis`, landscapes with redundant viable routes between the founder
and high-fitness basins have shorter forecast horizons: independent replicates
diverge earlier and are less likely to occupy the same basin. Funnelled
landscapes should remain more reproducible even when CIN is high.

## Analysis units

* **Landscape topology:** a graph whose nodes are observed karyotypes and whose
  edges represent one-missegregation adjacency.
* **Transition-flow topology:** a directed, weighted graph whose edges are
  realized parent-to-daughter transitions. This is the preferred graph for the
  primary analysis.
* **Forecastability:** agreement among independent ABM replicates through time,
  quantified by distributional divergence and basin concordance.

## Planned endpoints

1. Time-resolved replicate divergence (Jensen-Shannon divergence).
2. Basin concordance: probability that two replicates occupy the same dominant
   high-fitness basin at each time point.
3. Route entropy and effective number of flow routes.
4. Directed-flow concentration, bypass redundancy, and community persistence.
5. Forecast horizon: the earliest time at which replicate agreement crosses a
   pre-specified threshold.

## Data requirements

The current screen supplies final karyotypes and population summaries. That is
enough for an occupancy-topology inventory, implemented in
`01_build_occupancy_graphs.R`.

The primary analysis additionally requires sparse karyotype snapshots and
aggregate parent-to-daughter transition counts. The schema and a focused
follow-up design are recorded in `PROTOCOL.md`. Aggregate counts, not
cell-level lineage records, are sufficient.

## Interpretation guardrails

* Use weighted route ensembles, not only shortest paths.
* Never interpret a high-centrality karyotype as a therapeutic target.
* Treat longer bypass routes as the mechanism of interest.
* Pre-specify forecastability endpoints before selecting landscapes for the
  confirmatory follow-up.
