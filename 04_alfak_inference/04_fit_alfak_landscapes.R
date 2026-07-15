#!/usr/bin/env Rscript

# Fit one ALFA-K local landscape for every bounded-GRF landscape, p_mis value,
# and ABM replicate. Each fit uses the union of karyotypes observed at ABM
# steps 0, 1000, and 2000 in that replicate's phase-1 trajectory.

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]]) else getwd()
workers <- if (length(args) >= 2L) as.integer(args[[2]]) else 1L
nboot <- if (length(args) >= 3L) as.integer(args[[3]]) else 45L
fit_start <- if (length(args) >= 4L) as.integer(args[[4]]) else NA_integer_
fit_count <- if (length(args) >= 5L) as.integer(args[[5]]) else 1L
if (!is.finite(workers) || workers < 1L) stop("`workers` must be a positive integer.", call. = FALSE)
if (!is.na(fit_start) && (!is.finite(fit_start) || fit_start < 1L)) {
  stop("`fit_start` must be a positive integer.", call. = FALSE)
}
if (!is.finite(fit_count) || fit_count < 1L) {
  stop("`fit_count` must be a positive integer.", call. = FALSE)
}
source(file.path(project_dir, "R", "project_helpers.R"))
load_project_alfak(project_dir)

abm_dir <- file.path(project_dir, "outputs", "bounded_grf")
out_root <- file.path(project_dir, "outputs", "alfak_inference")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
observation_paths <- sort(Sys.glob(file.path(abm_dir, "landscape_*", "p_mis_*", "replicate_*", "abm_observations.rds")))
if (!length(observation_paths)) stop("No ABM observation files found. Run step 03 first.", call. = FALSE)
if (!is.na(fit_start)) {
  if (fit_start > length(observation_paths)) stop("`fit_start` exceeds the number of ABM observations.", call. = FALSE)
  fit_end <- min(length(observation_paths), fit_start + fit_count - 1L)
  observation_paths <- observation_paths[seq.int(fit_start, fit_end)]
}

fit_one <- function(observation_path) {
  replicate_dir <- dirname(observation_path)
  relative <- sub(paste0("^", abm_dir, "/"), "", replicate_dir)
  out_dir <- file.path(out_root, relative)
  metadata_path <- file.path(out_dir, "fit_metadata.rds")
  abm_pm <- read_abm_missegregation_rate(observation_path)
  expected_provenance <- list(
    fit_mode = alfak_fit_mode,
    observation_digest = file_digest(observation_path),
    observation_steps = alfak_observation_steps,
    karyotype_selection = alfak_karyotype_selection,
    nboot = nboot,
    minobs = alfak_minobs,
    n0 = 1e4,
    nb = 1e4,
    pm = abm_pm
  )
  if (file.exists(metadata_path)) {
    previous <- readRDS(metadata_path)
    if (provenance_matches(previous, expected_provenance)) {
      return(data.frame(status = "skipped", observation = relative))
    }
  }
  prepared <- prepare_observed_input(observation_path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    alfak(prepared$yi, outdir = out_dir, passage_times = prepared$passage_times,
          minobs = alfak_minobs, nboot = nboot, n0 = 1e4, nb = 1e4, pm = abm_pm,
          landscape_data_output = FALSE)
    saveRDS(list(observation_path = observation_path, fit_mode = prepared$fit_mode,
                 n_observed_karyotypes = length(prepared$observed_karyotypes),
                 observed_karyotypes = prepared$observed_karyotypes, nboot = nboot,
                 observation_steps = prepared$observation_steps,
                 karyotype_selection = prepared$karyotype_selection,
                 passage_times = prepared$passage_times,
                 minobs = alfak_minobs, n0 = 1e4, nb = 1e4, pm = abm_pm,
                 provenance = expected_provenance),
            file.path(out_dir, "fit_metadata.rds"))
    data.frame(status = "completed", observation = relative)
  }, error = function(e) data.frame(status = "failed", observation = relative, message = conditionMessage(e)))
}

fit_one_safely <- function(observation_path) {
  tryCatch(fit_one(observation_path), error = function(e) {
    relative <- sub(paste0("^", abm_dir, "/"), "", dirname(observation_path))
    data.frame(status = "failed", observation = relative, message = conditionMessage(e))
  })
}
status <- if (.Platform$OS.type == "windows" || workers <= 1L) {
  lapply(observation_paths, fit_one_safely)
} else {
  parallel::mclapply(observation_paths, fit_one_safely, mc.cores = workers, mc.preschedule = FALSE)
}
status <- bind_status(status)
status_path <- if (is.na(fit_start)) {
  file.path(out_root, "fit_status.csv")
} else {
  file.path(out_root, sprintf("fit_status_%05d_%05d.csv", fit_start, fit_start + length(observation_paths) - 1L))
}
utils::write.csv(status, status_path, row.names = FALSE)
if (any(status$status == "failed")) quit(status = 1L)
