#!/usr/bin/env Rscript

# Fit an ALFA-K landscape to every phase-2 ABM trajectory using the union of
# karyotypes observed at ABM steps 0, 1000, and 2000.
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]]) else getwd()
workers <- if (length(args) >= 2L) as.integer(args[[2]]) else 1L
nboot <- if (length(args) >= 3L) as.integer(args[[3]]) else 45L
fit_start <- if (length(args) >= 4L) as.integer(args[[4]]) else NA_integer_
fit_count <- if (length(args) >= 5L) as.integer(args[[5]]) else 1L
manifest_path <- if (length(args) >= 6L && nzchar(args[[6]])) normalizePath(args[[6]], mustWork = FALSE) else NA_character_
if (!is.finite(workers) || workers < 1L || !is.finite(nboot) || nboot < 1L) stop("Invalid workers or nboot.", call. = FALSE)
if (!is.na(manifest_path) && !file.exists(manifest_path)) {
  stop("ALFAK manifest does not exist: ", manifest_path, call. = FALSE)
}
source(file.path(project_dir, "R", "project_helpers.R"))
load_project_alfak(project_dir)

abm_root <- file.path(project_dir, "outputs", "phase2_abm")
out_root <- file.path(project_dir, "outputs", "phase2_alfak_inference")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
paths <- sort(Sys.glob(file.path(abm_root, "landscape_*", "p_mis_*", "replicate_*", "p_mis_*", "abm_observations.rds")))
if (!length(paths)) stop("No phase-2 observations found. Run phase-2 ABM first.", call. = FALSE)
if (!is.na(manifest_path)) {
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  if (!all(c("observation_path", "relative") %in% names(manifest))) {
    stop("Invalid ALFAK manifest: ", manifest_path, call. = FALSE)
  }
  if (!nrow(manifest)) stop("ALFAK manifest has no pending fits: ", manifest_path, call. = FALSE)
  paths <- manifest$observation_path
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("ALFAK manifest references missing observation file(s); first missing: ",
         missing[[1L]], call. = FALSE)
  }
  if (is.na(fit_start)) fit_start <- 1L
  if (!is.finite(fit_start) || fit_start < 1L || fit_start > length(paths) || !is.finite(fit_count) || fit_count < 1L) stop("Invalid phase-2 manifest fit range.", call. = FALSE)
  paths <- paths[seq.int(fit_start, min(length(paths), fit_start + fit_count - 1L))]
} else if (!is.na(fit_start)) {
  if (!is.finite(fit_start) || fit_start < 1L || fit_start > length(paths) || !is.finite(fit_count) || fit_count < 1L) stop("Invalid phase-2 fit range.", call. = FALSE)
  paths <- paths[seq.int(fit_start, min(length(paths), fit_start + fit_count - 1L))]
}

fit_one <- function(path) {
  relative <- sub(paste0("^", abm_root, "/"), "", dirname(path))
  out_dir <- file.path(out_root, relative)
  metadata_path <- file.path(out_dir, "fit_metadata.rds")
  abm_pm <- read_abm_missegregation_rate(path)
  expected_provenance <- alfak_expected_provenance(path, nboot)
  completion <- alfak_fit_completion_status(path, abm_root, out_root, nboot)
  if (isTRUE(completion$complete)) {
    previous <- readRDS(completion$metadata_path)
    attempted <- if (!is.null(previous$attempted_minobs)) previous$attempted_minobs else previous$minobs
    selected <- if (!is.null(previous$selected_minobs)) previous$selected_minobs else previous$minobs
    return(data.frame(status = "skipped", observation = relative,
                      selected_minobs = selected,
                      attempted_minobs = collapse_minobs_attempts(attempted),
                      fallback_used = isTRUE(previous$fallback_used)))
  }
  prepared <- prepare_observed_input(path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fit_result <- fit_alfak_with_minobs_fallback(
    prepared$yi, outdir = out_dir, passage_times = prepared$passage_times,
    minobs_candidates = alfak_minobs_candidates, nboot = nboot, n0 = 1e4,
    nb = 1e4, pm = abm_pm, landscape_data_output = FALSE
  )
  if (!fit_result$success) {
    return(data.frame(status = "failed", observation = relative,
                      selected_minobs = NA_integer_,
                      attempted_minobs = collapse_minobs_attempts(fit_result$attempted_minobs),
                      fallback_used = fit_result$fallback_used,
                      message = fit_result$message))
  }
  saveRDS(list(observation_path = path, fit_mode = prepared$fit_mode,
               n_observed_karyotypes = length(prepared$observed_karyotypes),
               observed_karyotypes = prepared$observed_karyotypes,
               nboot = nboot, observation_steps = prepared$observation_steps,
               karyotype_selection = prepared$karyotype_selection,
               passage_times = prepared$passage_times,
               minobs = fit_result$selected_minobs,
               selected_minobs = fit_result$selected_minobs,
               default_minobs = alfak_minobs,
               attempted_minobs = fit_result$attempted_minobs,
               minobs_candidates = alfak_minobs_candidates,
               minobs_strategy = alfak_minobs_strategy,
               minobs_fallback_policy = alfak_minobs_fallback_policy,
               fallback_used = fit_result$fallback_used,
               fallback_attempts = fit_result$attempts,
               n0 = 1e4, nb = 1e4, pm = abm_pm, provenance = expected_provenance),
          metadata_path)
  data.frame(status = "completed", observation = relative,
             selected_minobs = fit_result$selected_minobs,
             attempted_minobs = collapse_minobs_attempts(fit_result$attempted_minobs),
             fallback_used = fit_result$fallback_used)
}
safe_fit <- function(path) tryCatch(fit_one(path), error = function(e) data.frame(status = "failed", observation = sub(paste0("^", abm_root, "/"), "", dirname(path)), message = conditionMessage(e)))
status <- if (.Platform$OS.type == "windows" || workers <= 1L) lapply(paths, safe_fit) else parallel::mclapply(paths, safe_fit, mc.cores = workers, mc.preschedule = FALSE)
status <- bind_status(status)
name <- if (!is.na(manifest_path)) {
  sprintf("fit_status_manifest_%06d_%06d.csv", fit_start, fit_start + length(paths) - 1L)
} else if (is.na(fit_start)) {
  "fit_status.csv"
} else {
  sprintf("fit_status_%06d_%06d.csv", fit_start, fit_start + length(paths) - 1L)
}
utils::write.csv(status, file.path(out_root, name), row.names = FALSE)
if (any(status$status == "failed")) quit(status = 1L)
