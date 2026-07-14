# Input schemas for the confirmatory topology follow-up

These files are deliberately not synthesized. They are produced only by an
instrumented follow-up ABM after the discovery screen has selected settings.

## `karyotype_snapshots.csv`

One row per occupied karyotype at every recording time.

| column | meaning |
|---|---|
| `landscape_id` | e.g. `landscape_03` |
| `p_mis` | per-chromosome per-division probability |
| `replicate_id` | independent RNG replicate |
| `time_days` | observation time |
| `karyotype` | dot-delimited copy-number vector |
| `count` | cells with this karyotype |
| `fitness` | GRF fitness of this karyotype |

## `transition_flows.csv`

One row per aggregate realized parent-to-daughter transition within a recording
interval. Counts are aggregated; individual cell lineages are never stored.

| column | meaning |
|---|---|
| `landscape_id`, `p_mis`, `replicate_id` | simulation identifiers |
| `interval_start_day`, `interval_end_day` | aggregation interval |
| `parent_karyotype`, `daughter_karyotype` | directed transition |
| `daughter_count` | realized number of daughter cells on that edge |

The primary analysis excludes faithful self-transitions from route metrics, but
they may be retained in the raw file for auditing.
