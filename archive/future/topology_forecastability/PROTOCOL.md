# Pre-specified topology and forecastability protocol

## Discovery phase

Use the completed 200-landscape x 20-rate x 10-replicate phase-1 `p_mis` screen to create
an occupancy-topology inventory. This phase is descriptive and selects no
intervention.

For each landscape-rate pair, construct the induced graph of all karyotypes
observed in the final populations across the three replicates. An undirected
edge denotes Manhattan distance one (one chromosome gain or loss). Summarize:

* occupied-node count and component structure;
* largest-component mass fraction;
* degree and occupancy entropy;
* number of high-fitness components;
* founder-to-peak connectivity when observed nodes support it.

## Confirmatory follow-up

Select 9-12 landscape-rate settings prospectively across low, medium, and high
route-redundancy strata. For each, run 10-20 independent replicates with:

* a sparse state snapshot every 5 days;
* aggregate directed parent-to-daughter edge counts over each 5-day interval;
* unchanged biological parameters from the validated main ABM.

No individual karyotype is killed, penalized, or nominated as a drug target.

## Graph quantities

For a directed transition-flow graph, edge weight `F_ij` is the aggregate number
of daughters of parent state `i` entering daughter state `j` in an interval.

* **Route entropy:** entropy of normalized source-to-basin transition-path
  fluxes; report its exponential as the effective number of routes.
* **Flow concentration:** Herfindahl concentration of edge flux.
* **Bypass redundancy:** number of edge-disjoint positive-flow routes from the
  founder basin to a destination basin, plus the fraction of flux outside the
  highest-flow route.
* **Community persistence:** adjusted agreement between flow communities in
  adjacent time intervals.

## Forecastability endpoints

For each replicate pair at each recorded time:

1. Jensen-Shannon divergence between their karyotype frequency distributions.
2. Indicator of whether their dominant states belong to the same high-fitness
   graph component.

The forecast horizon is the first recorded time at which the median pairwise
Jensen-Shannon divergence exceeds a pre-specified threshold. Sensitivity to the
threshold is reported rather than optimized.

## Statistical test

Use a mixed model or permutation test with landscape as the resampling unit:

`forecast_horizon ~ p_mis * route_entropy + p_mis * bypass_redundancy`

Adjust for landscape wavelength and founder fitness. The claim is supported only
if topology explains forecastability beyond `p_mis` and roughness alone.
