#!/usr/bin/env Rscript

# Three-way topology analysis: original GRF, phase-1 inference, and phase-2
# inference. Original-GRF topology is evaluated on observed supports only;
# the complete 22-dimensional bounded GRF cannot be enumerated.
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]]) else getwd()
landscape_index <- if (length(args) >= 2L) as.integer(args[[2]]) else NA_integer_
if (is.na(landscape_index) || !is.finite(landscape_index) || landscape_index < 1L || landscape_index > 200L) stop("Run one landscape index (1-200) per task.", call. = FALSE)
if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.", call. = FALSE)

phase1_fit_root <- file.path(project_dir, "outputs", "alfak_inference")
phase2_fit_root <- file.path(project_dir, "outputs", "phase2_alfak_inference")
phase2_abm_root <- file.path(project_dir, "outputs", "phase2_abm")
out_dir <- file.path(project_dir, "outputs", "phase2_topology")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parse_karyotypes <- function(tags) {
  x <- do.call(rbind, strsplit(tags, ".", fixed = TRUE))
  storage.mode(x) <- "numeric"
  x
}
grf_fitness <- function(tags, landscape) {
  karyotypes <- parse_karyotypes(tags)
  vapply(seq_len(nrow(karyotypes)), function(i) {
    d <- sqrt(rowSums((landscape$centroids - rep(karyotypes[i, ], each = nrow(landscape$centroids)))^2))
    sum(sin(d / landscape$lambda)) / (pi * sqrt(nrow(landscape$centroids)))
  }, numeric(1))
}
build_edges <- function(tags) {
  index <- setNames(seq_along(tags), tags)
  states <- parse_karyotypes(tags)
  rows <- lapply(seq_along(tags), function(i) {
    candidates <- character()
    for (j in seq_len(ncol(states))) {
      plus <- states[i, ]; plus[j] <- plus[j] + 1L
      candidates <- c(candidates, paste(plus, collapse = "."))
      if (states[i, j] > 0) { minus <- states[i, ]; minus[j] <- minus[j] - 1L; candidates <- c(candidates, paste(minus, collapse = ".")) }
    }
    child <- unique(candidates[candidates %in% tags])
    child <- child[index[child] > i]
    if (length(child)) data.frame(from = tags[i], to = child, stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}
topology_metrics <- function(tags, fitness) {
  edges <- build_edges(tags)
  if (is.null(edges)) edges <- data.frame(from = character(), to = character())
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = tags))
  adjacent <- igraph::as_adj_list(graph, mode = "all")
  names(fitness) <- tags
  maxima <- vapply(seq_along(tags), function(i) {
    neighbours <- names(adjacent[[i]])
    !length(neighbours) || all(fitness[[tags[i]]] >= fitness[neighbours])
  }, logical(1))
  c(nodes = length(tags), edges = nrow(edges), components = igraph::components(graph)$no,
    local_maxima = sum(maxima), edge_density = if (length(tags) > 1L) 2 * nrow(edges) / (length(tags) * (length(tags) - 1L)) else 0)
}
spearman <- function(x, y) if (length(x) < 2L || !all(is.finite(x)) || !all(is.finite(y))) NA_real_ else suppressWarnings(stats::cor(x, y, method = "spearman"))
prefix_metrics <- function(prefix, values) stats::setNames(as.list(values), paste0(prefix, "_", names(values)))
clean_map <- function(x, path) {
  if (!is.data.frame(x) || !all(c("k", "mean") %in% names(x))) stop("Invalid inferred landscape: ", path, call. = FALSE)
  x <- x[is.finite(x$mean) & !duplicated(x$k), c("k", "mean")]
  if (nrow(x) < 2L) stop("Too little inferred support: ", path, call. = FALSE)
  x
}

landscape_id <- sprintf("landscape_%02d", landscape_index)
truth <- readRDS(file.path(project_dir, "data", "landscapes", paste0(landscape_id, ".rds")))
paths <- sort(Sys.glob(file.path(phase2_fit_root, landscape_id, "p_mis_*", "replicate_*", "p_mis_*", "landscape.Rds")))
if (!length(paths)) stop("No phase-2 inferred landscapes found for ", landscape_id, call. = FALSE)
metric_rows <- list(); summary_rows <- list()
for (path in paths) {
  relative <- sub(paste0("^", phase2_fit_root, "/"), "", dirname(path))
  bits <- strsplit(relative, "/", fixed = TRUE)[[1]]
  p1_dir <- bits[[2]]; replicate_id <- bits[[3]]; p2_dir <- bits[[4]]
  phase2 <- clean_map(readRDS(path), path)
  phase1_path <- file.path(phase1_fit_root, landscape_id, p1_dir, replicate_id, "landscape.Rds")
  phase1 <- clean_map(readRDS(phase1_path), phase1_path)
  shared <- intersect(phase1$k, phase2$k)
  union <- union(phase1$k, phase2$k)
  p1_fit_shared <- phase1$mean[match(shared, phase1$k)]
  p2_fit_shared <- phase2$mean[match(shared, phase2$k)]
  truth_p1 <- grf_fitness(phase1$k, truth)
  truth_p2 <- grf_fitness(phase2$k, truth)
  truth_union <- grf_fitness(union, truth)
  support <- list(
    phase1_nodes = nrow(phase1), phase2_nodes = nrow(phase2), shared_nodes = length(shared), union_nodes = length(union),
    phase1_only_nodes = length(setdiff(phase1$k, phase2$k)), phase2_only_nodes = length(setdiff(phase2$k, phase1$k)),
    jaccard = length(shared) / length(union), phase1_retained = length(shared) / nrow(phase1), phase2_inherited = length(shared) / nrow(phase2)
  )
  row <- c(list(landscape_id = landscape_id, p1_dir = p1_dir, replicate_id = replicate_id, p2_dir = p2_dir),
    support,
    prefix_metrics("phase1_inferred", topology_metrics(phase1$k, phase1$mean)),
    prefix_metrics("phase2_inferred", topology_metrics(phase2$k, phase2$mean)),
    prefix_metrics("original_grf_on_phase1_support", topology_metrics(phase1$k, truth_p1)),
    prefix_metrics("original_grf_on_phase2_support", topology_metrics(phase2$k, truth_p2)),
    prefix_metrics("original_grf_on_union_support", topology_metrics(union, truth_union)),
    list(
      phase1_original_grf_spearman = spearman(phase1$mean, truth_p1),
      phase2_original_grf_spearman = spearman(phase2$mean, truth_p2),
      phase1_phase2_shared_spearman = spearman(p1_fit_shared, p2_fit_shared),
      phase1_phase2_shared_rmse = if (length(shared)) sqrt(mean((p1_fit_shared - p2_fit_shared)^2)) else NA_real_
    ))
  metric_rows[[relative]] <- as.data.frame(row, stringsAsFactors = FALSE)
  final_path <- file.path(phase2_abm_root, relative, "final_karyotypes.csv")
  if (file.exists(final_path)) {
    final <- read.csv(final_path, stringsAsFactors = FALSE)
    original_fitness <- grf_fitness(final$karyotype, truth)
    summary_rows[[relative]] <- data.frame(landscape_id, p1_dir, replicate_id, p2_dir,
      terminal_population = sum(final$count), terminal_diversity = sum(final$count > 0),
      mean_source_inferred_fitness = stats::weighted.mean(final$fitness, final$count),
      mean_original_grf_fitness = stats::weighted.mean(original_fitness, final$count))
  }
}
utils::write.csv(do.call(rbind, metric_rows), file.path(out_dir, sprintf("phase2_topology_metrics_%s.csv", landscape_id)), row.names = FALSE)
utils::write.csv(do.call(rbind, summary_rows), file.path(out_dir, sprintf("phase2_terminal_fitness_%s.csv", landscape_id)), row.names = FALSE)
