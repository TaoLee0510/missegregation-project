#!/usr/bin/env Rscript

# Derive a biological reference ploidy from the longest untreated lineage in
# each available PDX.  PDXs contribute equally, preventing the extensively
# sampled SA609 lineage from dominating the reference.

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]]) else getwd()
salehi_dir <- file.path(project_dir, "data", "salehi_reference")
raw_dir <- file.path(salehi_dir, "raw")
processed_dir <- file.path(salehi_dir, "processed")
input_dir <- file.path(processed_dir, "alfak_inputs")
output_dir <- file.path(project_dir, "data", "reference_ploidy")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(file.path(raw_dir, "metadata.csv"), stringsAsFactors = FALSE)
lineages <- readRDS(file.path(processed_dir, "lineages.Rds"))
arm_loci <- readRDS(file.path(raw_dir, "arm_loci.Rds"))

# The bundled loci are GRCh38 coordinates.  The endpoint of each autosome's
# q arm is its chromosome length; sex chromosomes are intentionally excluded,
# matching the 22-chromosome ALFA-K input profiles.
length_table <- aggregate(end ~ chrom, data = subset(arm_loci, chrom %in% as.character(1:22)), FUN = max)
chromosome_lengths <- length_table$end[match(as.character(seq_len(22L)), length_table$chrom)]
stopifnot(length(chromosome_lengths) == 22L, all(is.finite(chromosome_lengths)))

lineage_table <- data.frame(
  lineage_id = names(lineages),
  n_passages = vapply(lineages, function(x) length(x$ids), integer(1)),
  endpoint_uid = vapply(lineages, function(x) tail(x$ids, 1L), numeric(1)),
  stringsAsFactors = FALSE
)
lineage_table$pdx_id <- metadata$PDX_id[match(lineage_table$endpoint_uid, metadata$uid)]
lineage_table$untreated_entire_chain <- vapply(lineages, function(x) {
  all(metadata$on_treatment[match(x$ids, metadata$uid)] == "n")
}, logical(1))
lineage_table$has_alfak_input <- file.exists(file.path(input_dir, paste0(lineage_table$lineage_id, ".Rds")))

# SA039 and SA906 are the in-vitro controls in this repository, not PDXs.
candidates <- subset(
  lineage_table,
  untreated_entire_chain & has_alfak_input & !(pdx_id %in% c("SA039", "SA906"))
)
if (!nrow(candidates)) stop("No untreated PDX lineages with ALFA-K inputs were found.", call. = FALSE)

# Retain every distinct PDX's maximal untreated chain.  With the bundled data
# this yields four PDXs: SA609 (10 passages), SA532 (8), SA1035 (5), and SA535
# (5).  This avoids an arbitrary tie break between SA1035 and SA535.
selected <- do.call(rbind, lapply(split(candidates, candidates$pdx_id), function(x) {
  x[x$n_passages == max(x$n_passages), , drop = FALSE]
}))
selected <- selected[order(-selected$n_passages, selected$pdx_id, selected$lineage_id), , drop = FALSE]
if (anyDuplicated(selected$pdx_id)) {
  stop("A PDX has multiple equally long eligible chains; specify a tie rule before deriving a reference.", call. = FALSE)
}

weighted_ploidy <- function(karyotype_matrix) {
  drop(karyotype_matrix %*% chromosome_lengths / sum(chromosome_lengths))
}

summarise_lineage <- function(lineage_id) {
  x <- readRDS(file.path(input_dir, paste0(lineage_id, ".Rds")))$x
  karyotypes <- do.call(rbind, strsplit(rownames(x), ".", fixed = TRUE))
  storage.mode(karyotypes) <- "numeric"
  if (ncol(karyotypes) != 22L) stop(sprintf("%s does not contain 22 autosomes.", lineage_id), call. = FALSE)
  k_ploidy <- weighted_ploidy(karyotypes)
  passage_ploidy <- drop(crossprod(k_ploidy, x) / colSums(x))
  data.frame(
    lineage_id = lineage_id,
    passage = as.numeric(colnames(x)),
    passage_weighted_ploidy = passage_ploidy,
    passage_cells = colSums(x),
    stringsAsFactors = FALSE
  )
}

passage_summary <- do.call(rbind, lapply(selected$lineage_id, summarise_lineage))
lineage_summary <- do.call(rbind, lapply(split(passage_summary, passage_summary$lineage_id), function(x) {
  data.frame(
    lineage_id = x$lineage_id[[1L]],
    pooled_weighted_ploidy = weighted.mean(x$passage_weighted_ploidy, x$passage_cells),
    initial_weighted_ploidy = x$passage_weighted_ploidy[[1L]],
    final_weighted_ploidy = x$passage_weighted_ploidy[[nrow(x)]],
    total_cells = sum(x$passage_cells),
    stringsAsFactors = FALSE
  )
}))
selected <- merge(selected, lineage_summary, by = "lineage_id", sort = FALSE)
selected <- selected[match(unique(passage_summary$lineage_id), selected$lineage_id), , drop = FALSE]

reference_mean_ploidy <- mean(selected$pooled_weighted_ploidy)

# A broad, chromosome-specific support is safer than one global copy-number
# cap.  The upper bound is the empirical 99.5th percentile plus one copy: the
# extra copy prevents the observed tail itself from becoming the simulation
# boundary.  Lower bounds remain zero because nullisomies occur in the source
# data; their biological consequences are handled by GRF fitness.
cell_copy_numbers <- do.call(rbind, lapply(selected$lineage_id, function(lineage_id) {
  x <- readRDS(file.path(input_dir, paste0(lineage_id, ".Rds")))$x
  k <- do.call(rbind, strsplit(rownames(x), ".", fixed = TRUE))
  storage.mode(k) <- "numeric"
  k[rep(seq_len(nrow(k)), rowSums(x)), , drop = FALSE]
}))
upper_quantile <- 0.995
upper_copy_numbers <- ceiling(apply(cell_copy_numbers, 2, stats::quantile, probs = upper_quantile, names = FALSE)) + 1L
lower_copy_numbers <- rep.int(0L, 22L)
bound_profile <- data.frame(
  chromosome = seq_len(22L),
  lower_copy_number = lower_copy_numbers,
  upper_copy_number = upper_copy_numbers,
  empirical_upper_quantile = upper_quantile,
  stringsAsFactors = FALSE
)
reference <- list(
  reference_mean_ploidy = reference_mean_ploidy,
  aggregation = "unweighted mean of one pooled, cell-count-weighted lineage ploidy per untreated PDX",
  chromosome_lengths_source = "data/salehi_reference/raw/arm_loci.Rds (GRCh38 autosomal q-arm endpoints)",
  selected_lineages = selected,
  bound_profile = bound_profile
)

utils::write.csv(selected, file.path(output_dir, "selected_untreated_pdx_lineages.csv"), row.names = FALSE)
utils::write.csv(passage_summary, file.path(output_dir, "passage_weighted_ploidy.csv"), row.names = FALSE)
utils::write.csv(bound_profile, file.path(output_dir, "empirical_copy_number_bounds.csv"), row.names = FALSE)
saveRDS(reference, file.path(output_dir, "reference_ploidy.rds"))
writeLines(sprintf("%.10f", reference_mean_ploidy), file.path(output_dir, "reference_mean_ploidy.txt"))
message(sprintf("Reference mean ploidy: %.6f", reference_mean_ploidy))
