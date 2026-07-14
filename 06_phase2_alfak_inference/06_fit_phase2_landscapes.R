#!/usr/bin/env Rscript

# Fit an ALFA-K landscape to every phase-2 ABM trajectory using every nonzero
# terminal karyotype and its full recorded trajectory.
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]]) else getwd()
workers <- if (length(args) >= 2L) as.integer(args[[2]]) else 1L
nboot <- if (length(args) >= 3L) as.integer(args[[3]]) else 45L
fit_start <- if (length(args) >= 4L) as.integer(args[[4]]) else NA_integer_
fit_count <- if (length(args) >= 5L) as.integer(args[[5]]) else 1L
if (!is.finite(workers) || workers < 1L || !is.finite(nboot) || nboot < 1L) stop("Invalid workers or nboot.", call. = FALSE)
library(alfakR)
source(file.path(project_dir, "R", "project_helpers.R"))

abm_root <- file.path(project_dir, "outputs", "phase2_abm")
out_root <- file.path(project_dir, "outputs", "phase2_alfak_inference")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
paths <- sort(Sys.glob(file.path(abm_root, "landscape_*", "p_mis_*", "replicate_*", "p_mis_*", "abm_observations.rds")))
if (!length(paths)) stop("No phase-2 observations found. Run phase-2 ABM first.", call. = FALSE)
if (!is.na(fit_start)) {
  if (!is.finite(fit_start) || fit_start < 1L || fit_start > length(paths) || !is.finite(fit_count) || fit_count < 1L) stop("Invalid phase-2 fit range.", call. = FALSE)
  paths <- paths[seq.int(fit_start, min(length(paths), fit_start + fit_count - 1L))]
}

fit_one <- function(path) {
  relative <- sub(paste0("^", abm_root, "/"), "", dirname(path))
  out_dir <- file.path(out_root, relative)
  metadata_path <- file.path(out_dir, "fit_metadata.rds")
  if (file.exists(metadata_path)) {
    previous <- readRDS(metadata_path)
    if (identical(previous$fit_mode, "all_terminal_karyotypes_v1") && identical(previous$nboot, nboot)) return(data.frame(status = "skipped", observation = relative))
  }
  prepared <- prepare_terminal_input(path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  alfak(prepared$yi, outdir = out_dir, passage_times = prepared$passage_times, minobs = 1L, nboot = nboot, n0 = 1e4, nb = 1e4, pm = 5e-05, landscape_data_output = FALSE)
  saveRDS(list(observation_path = path, fit_mode = "all_terminal_karyotypes_v1", n_terminal_karyotypes = length(prepared$terminal), terminal_karyotypes = prepared$terminal, nboot = nboot), metadata_path)
  data.frame(status = "completed", observation = relative)
}
safe_fit <- function(path) tryCatch(fit_one(path), error = function(e) data.frame(status = "failed", observation = sub(paste0("^", abm_root, "/"), "", dirname(path)), message = conditionMessage(e)))
status <- if (.Platform$OS.type == "windows" || workers <= 1L) lapply(paths, safe_fit) else parallel::mclapply(paths, safe_fit, mc.cores = workers, mc.preschedule = FALSE)
status <- bind_status(status)
name <- if (is.na(fit_start)) "fit_status.csv" else sprintf("fit_status_%06d_%06d.csv", fit_start, fit_start + length(paths) - 1L)
utils::write.csv(status, file.path(out_root, name), row.names = FALSE)
if (any(status$status == "failed")) quit(status = 1L)
