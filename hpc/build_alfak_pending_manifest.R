#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]]) else getwd()
phase <- if (length(args) >= 2L) args[[2L]] else "phase1"
manifest_path <- if (length(args) >= 3L) args[[3L]] else file.path(project_dir, "results", "manifests", paste0(phase, "_alfak_pending.csv"))
nboot <- if (length(args) >= 4L) as.integer(args[[4L]]) else 45L
if (!phase %in% c("phase1", "phase2")) stop("`phase` must be 'phase1' or 'phase2'.", call. = FALSE)
if (!is.finite(nboot) || nboot < 1L) stop("`nboot` must be a positive integer.", call. = FALSE)

source(file.path(project_dir, "R", "project_helpers.R"))

if (phase == "phase1") {
  abm_root <- file.path(project_dir, "outputs", "bounded_grf")
  out_root <- file.path(project_dir, "outputs", "alfak_inference")
  pattern <- file.path(abm_root, "landscape_*", "p_mis_*", "replicate_*", "abm_observations.rds")
} else {
  abm_root <- file.path(project_dir, "outputs", "phase2_abm")
  out_root <- file.path(project_dir, "outputs", "phase2_alfak_inference")
  pattern <- file.path(abm_root, "landscape_*", "p_mis_*", "replicate_*", "p_mis_*", "abm_observations.rds")
}

observation_paths <- sort(Sys.glob(pattern))
if (!length(observation_paths)) stop("No ", phase, " ABM observation files found under ", abm_root, call. = FALSE)

empty_pending <- data.frame(
  task_index = integer(),
  phase = character(),
  source_index = integer(),
  relative = character(),
  observation_path = character(),
  out_dir = character(),
  metadata_path = character(),
  landscape_path = character(),
  reason = character(),
  stringsAsFactors = FALSE
)

rows <- vector("list", length(observation_paths))
complete <- logical(length(observation_paths))
reasons <- character(length(observation_paths))
for (i in seq_along(observation_paths)) {
  status <- alfak_fit_completion_status(observation_paths[[i]], abm_root, out_root, nboot)
  complete[[i]] <- isTRUE(status$complete)
  reasons[[i]] <- status$reason
  if (!complete[[i]]) {
    rows[[i]] <- data.frame(
      task_index = NA_integer_,
      phase = phase,
      source_index = i,
      relative = status$relative,
      observation_path = observation_paths[[i]],
      out_dir = status$out_dir,
      metadata_path = status$metadata_path,
      landscape_path = status$landscape_path,
      reason = status$reason,
      stringsAsFactors = FALSE
    )
  }
}

pending <- Filter(Negate(is.null), rows)
pending <- if (length(pending)) do.call(rbind, pending) else empty_pending
if (nrow(pending)) pending$task_index <- seq_len(nrow(pending))
write_csv_atomic(pending, manifest_path, row.names = FALSE)

reason_counts <- as.data.frame(table(reason = reasons), stringsAsFactors = FALSE)
summary <- data.frame(
  phase = phase,
  total_observations = length(observation_paths),
  completed = sum(complete),
  pending = nrow(pending),
  manifest_path = manifest_path,
  stringsAsFactors = FALSE
)
summary_path <- sub("\\.csv$", "_summary.csv", manifest_path)
write_csv_atomic(summary, summary_path, row.names = FALSE)
write_csv_atomic(reason_counts, sub("\\.csv$", "_reason_counts.csv", manifest_path), row.names = FALSE)

cat(sprintf("phase=%s\n", phase))
cat(sprintf("total_observations=%d\n", length(observation_paths)))
cat(sprintf("completed=%d\n", sum(complete)))
cat(sprintf("pending=%d\n", nrow(pending)))
cat(sprintf("manifest=%s\n", manifest_path))
