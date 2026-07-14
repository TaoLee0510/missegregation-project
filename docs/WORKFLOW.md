# Full-GRF validation experiment

`03_bounded_grf_abm/03_run_bounded_grf_abm.R` runs a chromosome-specifically
bounded GRF ABM. It starts 10,000 cells across 50 in-bounds FQ states sampled
near the empirical reference ploidy, then evaluates GRF fitness for every
in-bounds descendant. The bounds are broad PDX-derived support—not an FQ + NN
lookup-table shell—and boundary rejections are reported in run metadata.
The all-diploid state is excluded from both initialization and descendants.

`01_reference_ploidy/01_derive_reference_ploidy.R` derives the empirical
target ploidy from the longest entirely untreated lineage of every PDX in the
bundled Salehi data. It uses chromosome lengths from the bundled GRCh38 loci
and gives every PDX equal weight. The current four-PDX reference is 2.413768.
