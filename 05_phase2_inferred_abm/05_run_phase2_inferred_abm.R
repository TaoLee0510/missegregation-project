#!/usr/bin/env Rscript

# Continue each phase-1 replicate on its own ALFA-K-inferred, finite fitness
# support across every second-stage p_mis value.  Replicate IDs are inherited,
# yielding 10 matched replicate trajectories per (p_mis_1, p_mis_2) pair.

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]]) else getwd()
workers <- if (length(args) >= 2L) as.integer(args[[2]]) else 1L
landscape_index <- if (length(args) >= 3L) as.integer(args[[3]]) else NA_integer_
p1_index <- if (length(args) >= 4L) as.integer(args[[4]]) else NA_integer_
replicate_filter <- if (length(args) >= 6L) as.integer(args[[5]]) else NA_integer_
p2_index <- if (length(args) >= 6L) as.integer(args[[6]]) else NA_integer_
n_steps <- if (length(args) >= 7L) {
  as.integer(args[[7]])
} else if (length(args) == 5L) {
  as.integer(args[[5]])
} else {
  2000L
}
if (!is.finite(workers) || workers < 1L) stop("`workers` must be a positive integer.", call. = FALSE)
if (!is.na(landscape_index) && (!is.finite(landscape_index) || landscape_index < 1L)) stop("`landscape_index` must be positive.", call. = FALSE)
if (!is.na(p1_index) && (!is.finite(p1_index) || p1_index < 1L)) stop("`p1_index` must be positive.", call. = FALSE)
if (!is.na(replicate_filter) && (!is.finite(replicate_filter) || replicate_filter < 1L)) stop("`replicate_id` must be positive.", call. = FALSE)
if (!is.na(p2_index) && (!is.finite(p2_index) || p2_index < 1L)) stop("`p2_index` must be positive.", call. = FALSE)
if (!is.finite(n_steps) || n_steps < 0L) stop("`n_steps` must be a non-negative integer.", call. = FALSE)
library(alfakR)
source(file.path(project_dir, "R", "project_helpers.R"))

phase1_root <- file.path(project_dir, "outputs", "bounded_grf")
phase1_fit_root <- file.path(project_dir, "outputs", "alfak_inference")
out_root <- file.path(project_dir, "outputs", "phase2_abm")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
p_mis_path <- file.path(phase1_root, "p_mis_lhs.csv")
parameters <- read.csv(p_mis_path, stringsAsFactors = FALSE)
if (!all(c("p_index", "p_mis") %in% names(parameters))) stop("Invalid phase-1 p_mis_lhs.csv.", call. = FALSE)
parameters <- parameters[order(parameters$p_index), ]
phase1_parameters <- parameters
phase2_parameters <- parameters
p_mis_lhs_digest <- file_digest(p_mis_path)
n_replicates <- 10L
n_cells <- 10000L
dt <- 0.1
record_interval <- 50L
diploid_tag <- paste(rep.int(2L, 22L), collapse = ".")
model_version <- "phase2_full_endpoint_direct_map_v4_no_diploid_state"
allocation_method <- "largest_remainder_hash_ties_v1"

prepare_initial_population <- function(final_path) {
  final <- read.csv(final_path, stringsAsFactors = FALSE)
  required <- c("karyotype", "count")
  if (!all(required %in% names(final))) stop("Invalid phase-1 final population: ", final_path, call. = FALSE)
  endpoint <- final[final$count > 0, required]
  endpoint <- stats::aggregate(count ~ karyotype, endpoint, sum)
  endpoint <- endpoint[endpoint$karyotype != diploid_tag, , drop = FALSE]
  endpoint <- endpoint[order(endpoint$karyotype), , drop = FALSE]
  if (!nrow(endpoint)) stop("Phase-1 endpoint population is empty: ", final_path, call. = FALSE)
  if (nrow(endpoint) > n_cells) stop("Cannot retain every endpoint karyotype in a 10,000-cell phase-2 population.", call. = FALSE)
  weights <- stats::setNames(endpoint$count, endpoint$karyotype)
  counts <- allocate_counts_exact(weights, n_cells, min_count = 1L, tie_breaker = endpoint$karyotype)
  list(
    tags = endpoint$karyotype,
    counts = as.integer(counts[endpoint$karyotype]),
    source = "complete phase-1 terminal population, proportionally rescaled to 10,000 cells",
    provenance = list(
      allocation_method = allocation_method,
      source_final_digest = file_digest(final_path),
      n_cells = n_cells,
      excluded_karyotypes = diploid_tag
    )
  )
}

source_paths <- function(landscape_id, p1, replicate_id) {
  p1_dir <- sprintf("p_mis_%02d_%.8f", p1$p_index, p1$p_mis)
  replicate_dir <- sprintf("replicate_%02d", replicate_id)
  list(
    p1_dir = p1_dir,
    replicate_dir = replicate_dir,
    final_path = file.path(phase1_root, landscape_id, p1_dir, replicate_dir, "final_karyotypes.csv"),
    inferred_path = file.path(phase1_fit_root, landscape_id, p1_dir, replicate_dir, "landscape.Rds")
  )
}

prepare_source <- function(landscape_id, p1, replicate_id) {
  paths <- source_paths(landscape_id, p1, replicate_id)
  final_path <- paths$final_path
  inferred_path <- paths$inferred_path
  if (!file.exists(final_path) || !file.exists(inferred_path)) stop("Phase-1 endpoint and inferred landscape are required for ", landscape_id, "/", paths$p1_dir, "/", paths$replicate_dir, call. = FALSE)
  inferred <- readRDS(inferred_path)
  if (!is.data.frame(inferred) || !all(c("k", "mean") %in% names(inferred))) stop("Invalid phase-1 inferred landscape: ", inferred_path, call. = FALSE)
  inferred <- inferred[is.finite(inferred$mean) & !duplicated(inferred$k), c("k", "mean")]
  inferred <- inferred[inferred$k != diploid_tag, , drop = FALSE]
  initial <- prepare_initial_population(final_path)
  missing <- setdiff(initial$tags, inferred$k)
  if (length(missing)) stop("Phase-1 inference did not directly estimate every endpoint karyotype: ", paste(missing, collapse = ", "), call. = FALSE)
  fitness_map <- stats::setNames(as.list(inferred$mean), inferred$k)
  list(initial = initial, fitness_map = fitness_map, final_path = final_path, inferred_path = inferred_path)
}

manifest <- read.csv(file.path(project_dir, "data", "landscapes", "manifest.csv"), stringsAsFactors = FALSE)
if (!is.na(landscape_index)) {
  if (landscape_index > nrow(manifest)) stop("`landscape_index` exceeds the manifest.", call. = FALSE)
  manifest <- manifest[landscape_index, , drop = FALSE]
}
if (!is.na(p1_index)) {
  if (p1_index > nrow(phase1_parameters)) stop("`p1_index` exceeds the p_mis table.", call. = FALSE)
  phase1_parameters <- phase1_parameters[p1_index, , drop = FALSE]
}
if (!is.na(p2_index)) {
  if (p2_index > nrow(phase2_parameters)) stop("`p2_index` exceeds the p_mis table.", call. = FALSE)
  phase2_parameters <- phase2_parameters[p2_index, , drop = FALSE]
}
if (!is.na(replicate_filter) && replicate_filter > n_replicates) stop("`replicate_id` exceeds the replicate count.", call. = FALSE)
replicate_ids <- if (is.na(replicate_filter)) seq_len(n_replicates) else replicate_filter

tasks <- list()
for (landscape_id in manifest$landscape_id) for (i in seq_len(nrow(phase1_parameters))) for (replicate_id in replicate_ids) for (j in seq_len(nrow(phase2_parameters))) {
  tasks[[length(tasks) + 1L]] <- list(landscape_id = landscape_id, p1 = phase1_parameters[i, ], replicate_id = replicate_id, p2 = phase2_parameters[j, ])
}

run_task <- function(task) {
  p1_dir <- p_mis_dir(task$p1$p_index, task$p1$p_mis)
  p2_dir <- p_mis_dir(task$p2$p_index, task$p2$p_mis)
  source_dir <- file.path(out_root, task$landscape_id, p1_dir, sprintf("replicate_%02d", task$replicate_id))
  out_dir <- file.path(source_dir, p2_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  metadata_path <- file.path(out_dir, "run_metadata.rds")
  paths <- source_paths(task$landscape_id, task$p1, task$replicate_id)
  stop_if_missing_files(c(paths$final_path, paths$inferred_path), "phase-1 source")
  expected_provenance <- list(
    model_version = model_version,
    allocation_method = allocation_method,
    task = list(
      landscape_id = task$landscape_id,
      p1_index = task$p1$p_index,
      p_mis_phase1 = task$p1$p_mis,
      p2_index = task$p2$p_index,
      p_mis_phase2 = task$p2$p_mis,
      replicate_id = task$replicate_id
    ),
    p_mis_lhs_digest = p_mis_lhs_digest,
    phase2_parameter_count = nrow(parameters),
    source_final_digest = file_digest(paths$final_path),
    source_inferred_digest = file_digest(paths$inferred_path),
    n_cells = n_cells,
    n_steps = n_steps,
    dt = dt,
    record_interval = record_interval,
    culling_survival_fraction = 0.999,
    excluded_karyotypes = diploid_tag
  )
  if (file.exists(metadata_path)) {
    previous <- readRDS(metadata_path)
    if (provenance_matches(previous, expected_provenance)) return(data.frame(status = "skipped", landscape_id = task$landscape_id, p1_index = task$p1$p_index, p2_index = task$p2$p_index, replicate_id = task$replicate_id))
  }
  source <- prepare_source(task$landscape_id, task$p1, task$replicate_id)
  initial_path <- file.path(source_dir, "phase2_initialization.rds")
  if (!file.exists(initial_path) || !provenance_matches(readRDS(initial_path), source$initial$provenance)) {
    write_rds_atomic(source$initial, initial_path)
  }
  seed <- 1200000L + as.integer(sub(".*_", "", task$landscape_id)) * 100000L + task$p1$p_index * 1000L + task$replicate_id * 100L + task$p2$p_index
  raw <- run_alfak_abm(
    initial_population_r = stats::setNames(as.list(source$initial$counts), source$initial$tags),
    fitness_map_r = source$fitness_map, p_missegregation = task$p2$p_mis,
    dt = dt, n_steps = n_steps, max_population_size = n_cells,
    culling_survival_fraction = 0.999, record_interval = record_interval, seed = seed,
    excluded_karyotypes = diploid_tag
  )
  saveRDS(raw_to_observations(raw, dt), file.path(out_dir, "abm_observations.rds"))
  final_counts <- raw[[tail(names(raw), 1L)]]
  final_fitness <- unlist(source$fitness_map[names(final_counts)], use.names = FALSE)
  utils::write.csv(data.frame(karyotype = names(final_counts), count = as.numeric(final_counts), fitness = final_fitness), file.path(out_dir, "final_karyotypes.csv"), row.names = FALSE)
  saveRDS(list(landscape_id = task$landscape_id, p_mis_phase1 = task$p1$p_mis, p_mis_phase2 = task$p2$p_mis,
               replicate_id = task$replicate_id, seed = seed, n_steps = n_steps, source_final_path = source$final_path,
               source_inferred_path = source$inferred_path,
               phase2_support = "phase-1 ALFA-K support plus every retained endpoint karyotype; diploid state excluded", excluded_karyotypes = diploid_tag,
               model_version = model_version, provenance = expected_provenance), metadata_path)
  data.frame(status = "completed", landscape_id = task$landscape_id, p1_index = task$p1$p_index, p2_index = task$p2$p_index, replicate_id = task$replicate_id)
}

safe_run_task <- function(task) tryCatch(run_task(task), error = function(e) data.frame(status = "failed", landscape_id = task$landscape_id, p1_index = task$p1$p_index, p2_index = task$p2$p_index, replicate_id = task$replicate_id, message = conditionMessage(e)))
status <- if (.Platform$OS.type == "windows" || workers <= 1L) lapply(tasks, safe_run_task) else parallel::mclapply(tasks, safe_run_task, mc.cores = workers, mc.preschedule = FALSE)
status <- bind_status(status)
suffix <- paste(
  if (is.na(landscape_index)) "all_landscapes" else sprintf("landscape_%02d", landscape_index),
  if (is.na(p1_index)) "all_p1" else sprintf("p1_%02d", p1_index),
  if (is.na(replicate_filter)) "all_replicates" else sprintf("replicate_%02d", replicate_filter),
  if (is.na(p2_index)) "all_p2" else sprintf("p2_%02d", p2_index),
  sep = "_"
)
utils::write.csv(status, file.path(out_root, paste0("run_status_", suffix, ".csv")), row.names = FALSE)
if (any(status$status == "failed")) quit(status = 1L)
