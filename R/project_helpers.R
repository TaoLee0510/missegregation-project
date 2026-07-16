raw_to_observations <- function(raw, dt) {
  steps <- names(raw)
  tags <- sort(unique(unlist(lapply(raw, names), use.names = FALSE)))
  x <- matrix(0L, nrow = length(tags), ncol = length(steps), dimnames = list(tags, steps))
  for (step in steps) x[names(raw[[step]]), step] <- as.integer(raw[[step]])
  list(x = x, dt = dt, passage_times = as.numeric(steps) * dt)
}

file_digest <- function(path) {
  if (!file.exists(path)) stop("Required file is missing: ", path, call. = FALSE)
  unname(tools::md5sum(path))
}

object_digest <- function(x) {
  path <- tempfile("object_digest_")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 2)
  file_digest(path)
}

write_rds_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(basename(path), "_"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Could not write RDS atomically: ", path, call. = FALSE)
  invisible(path)
}

write_csv_atomic <- function(object, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(basename(path), "_"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(object, tmp, ...)
  if (!file.rename(tmp, path)) stop("Could not write CSV atomically: ", path, call. = FALSE)
  invisible(path)
}

load_project_alfak <- function(project_dir) {
  source_path <- Sys.getenv("ALFAKR_SOURCE", unset = "")
  if (!nzchar(source_path)) {
    suppressPackageStartupMessages(library(alfakR))
    return(invisible("installed"))
  }
  source_path <- normalizePath(source_path, mustWork = TRUE)
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required when ALFAKR_SOURCE is set.", call. = FALSE)
  }
  compile <- tolower(Sys.getenv("ALFAKR_COMPILE", unset = "false")) %in% c("1", "true", "yes")
  if (!compile) {
    dll_path <- file.path(source_path, "src", paste0("alfakR", .Platform$dynlib.ext))
    if (!file.exists(dll_path)) {
      stop("Compiled repo alfakR DLL is missing: ", dll_path,
           "\nRun once before submitting: ALFAKR_SOURCE=/path/to/packages/alfakR ALFAKR_COMPILE=true Rscript -e 'source(\"R/project_helpers.R\"); load_project_alfak(getwd())'",
           call. = FALSE)
    }
  }
  pkgload::load_all(source_path, compile = compile, quiet = TRUE, export_all = FALSE,
                    helpers = FALSE, attach_testthat = FALSE)
  invisible(source_path)
}

provenance_matches <- function(previous, expected) {
  is.list(previous) && is.list(previous$provenance) && identical(previous$provenance, expected)
}

stable_string_hash <- function(x) {
  vapply(x, function(one) {
    ints <- utf8ToInt(one)
    if (!length(ints)) return(0)
    sum(as.double(ints) * seq_along(ints)) %% .Machine$integer.max
  }, numeric(1), USE.NAMES = FALSE)
}

allocate_counts_exact <- function(weights, total, min_count = 0L, tie_breaker = names(weights)) {
  if (!length(weights)) stop("`weights` must not be empty.", call. = FALSE)
  if (!is.finite(total) || total < 0 || total != floor(total)) stop("`total` must be a non-negative integer.", call. = FALSE)
  if (!is.finite(min_count) || min_count < 0 || min_count != floor(min_count)) stop("`min_count` must be a non-negative integer.", call. = FALSE)
  if (any(!is.finite(weights)) || any(weights < 0)) stop("`weights` must be finite and non-negative.", call. = FALSE)
  if (sum(weights) <= 0) stop("At least one weight must be positive.", call. = FALSE)
  base <- rep.int(as.numeric(min_count), length(weights))
  remaining <- as.numeric(total) - sum(base)
  if (remaining < 0) stop("`total` is smaller than the requested minimum counts.", call. = FALSE)
  raw_extra <- remaining * as.numeric(weights) / sum(weights)
  extra <- floor(raw_extra)
  remainder <- as.integer(remaining - sum(extra))
  if (remainder > 0L) {
    fractional <- raw_extra - extra
    if (is.null(tie_breaker) || length(tie_breaker) != length(weights)) tie_breaker <- as.character(seq_along(weights))
    order_index <- order(-fractional, stable_string_hash(as.character(tie_breaker)), seq_along(weights))
    extra[order_index[seq_len(remainder)]] <- extra[order_index[seq_len(remainder)]] + 1
  }
  counts <- base + extra
  names(counts) <- names(weights)
  counts
}

alfak_observation_steps <- c(0, 1000, 2000)
alfak_karyotype_selection <- "union_of_selected_timepoints"
alfak_minobs_candidates <- c(20L, 10L, 5L, 3L, 1L)
alfak_minobs <- alfak_minobs_candidates[[1L]]
alfak_minobs_fallback_policy <- "retry_all_alfak_fit_errors"
alfak_minobs_strategy <- paste0("fallback_", paste(alfak_minobs_candidates, collapse = "_"),
                                "_", alfak_minobs_fallback_policy)
alfak_fit_mode <- sprintf("timepoint_union_karyotypes_steps_0_1000_2000_minobs_%s_v2", alfak_minobs_strategy)

collapse_minobs_attempts <- function(x) paste(as.integer(x), collapse = ";")

fit_alfak_with_minobs_fallback <- function(yi, outdir, passage_times, nboot, n0, nb, pm,
                                           landscape_data_output = FALSE,
                                           minobs_candidates = alfak_minobs_candidates) {
  if (!length(minobs_candidates)) stop("`minobs_candidates` must not be empty.", call. = FALSE)
  minobs_candidates <- as.integer(minobs_candidates)
  if (any(is.na(minobs_candidates)) || any(minobs_candidates < 1L)) {
    stop("`minobs_candidates` must be positive integers.", call. = FALSE)
  }
  attempts <- list()
  for (i in seq_along(minobs_candidates)) {
    current_minobs <- minobs_candidates[[i]]
    if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    result <- tryCatch({
      alfak(yi, outdir = outdir, passage_times = passage_times,
            minobs = current_minobs, nboot = nboot, n0 = n0, nb = nb, pm = pm,
            landscape_data_output = landscape_data_output)
      NULL
    }, error = function(e) e)
    if (!inherits(result, "error")) {
      attempts[[length(attempts) + 1L]] <- data.frame(
        minobs = current_minobs, status = "completed", message = NA_character_
      )
      return(list(
        success = TRUE,
        selected_minobs = current_minobs,
        attempted_minobs = minobs_candidates[seq_len(i)],
        fallback_used = i > 1L,
        attempts = bind_status(attempts),
        message = NA_character_
      ))
    }
    message <- conditionMessage(result)
    attempts[[length(attempts) + 1L]] <- data.frame(
      minobs = current_minobs, status = "failed", message = message
    )
    if (i == length(minobs_candidates)) {
      if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
      return(list(
        success = FALSE,
        selected_minobs = NA_integer_,
        attempted_minobs = minobs_candidates[seq_len(i)],
        fallback_used = i > 1L,
        attempts = bind_status(attempts),
        message = message
      ))
    }
  }
  stop("unreachable ALFAK fallback state.", call. = FALSE)
}

prepare_observed_input <- function(observation_path, observation_steps = alfak_observation_steps) {
  observed <- readRDS(observation_path)
  if (!is.list(observed) || is.null(observed$x) || is.null(observed$dt) || is.null(observed$passage_times)) {
    stop("Invalid ABM observation object: ", observation_path, call. = FALSE)
  }
  x <- as.matrix(observed$x)
  if (is.null(rownames(x)) || any(!nzchar(rownames(x)))) {
    stop("ABM observation matrix must have karyotype row names: ", observation_path, call. = FALSE)
  }
  if (is.null(colnames(x)) || any(!nzchar(colnames(x)))) {
    stop("ABM observation matrix must have timepoint column names: ", observation_path, call. = FALSE)
  }
  if (length(observed$passage_times) != ncol(x)) {
    stop("ABM observation passage_times length does not match observation columns: ", observation_path, call. = FALSE)
  }
  step_labels <- as.character(observation_steps)
  keep_cols <- match(step_labels, colnames(x))
  if (anyNA(keep_cols)) {
    stop("ABM observation is missing required timepoint column(s): ",
         paste(step_labels[is.na(keep_cols)], collapse = ", "),
         " in ", observation_path, call. = FALSE)
  }
  x <- x[, keep_cols, drop = FALSE]
  passage_times <- observed$passage_times[keep_cols]
  # ALFAK input support is the union of karyotypes observed at the selected steps.
  observed_rows <- rowSums(x) > 0
  if (!any(observed_rows)) stop("No observed karyotypes for ", observation_path, call. = FALSE)
  x <- x[observed_rows, , drop = FALSE]
  zero_depth <- colSums(x) <= 0
  if (any(zero_depth)) {
    stop("ABM observation has zero-depth timepoint(s): ",
         paste(colnames(x)[zero_depth], collapse = ", "),
         " in ", observation_path, call. = FALSE)
  }
  list(yi = list(x = x, dt = observed$dt), passage_times = passage_times,
       observed_karyotypes = rownames(x),
       observation_steps = as.numeric(colnames(x)),
       karyotype_selection = alfak_karyotype_selection,
       fit_mode = alfak_fit_mode)
}

read_abm_missegregation_rate <- function(observation_path) {
  metadata_path <- file.path(dirname(observation_path), "run_metadata.rds")
  if (!file.exists(metadata_path)) stop("ABM run metadata is missing for ", observation_path, call. = FALSE)
  metadata <- readRDS(metadata_path)
  candidates <- list(
    if (!is.null(metadata$p_mis_phase2)) metadata$p_mis_phase2 else NULL,
    if (is.list(metadata$task) && !is.null(metadata$task$p_mis)) metadata$task$p_mis else NULL,
    if (is.list(metadata$provenance) && is.list(metadata$provenance$task) && !is.null(metadata$provenance$task$p_mis_phase2)) metadata$provenance$task$p_mis_phase2 else NULL,
    if (is.list(metadata$provenance) && is.list(metadata$provenance$task) && !is.null(metadata$provenance$task$p_mis)) metadata$provenance$task$p_mis else NULL
  )
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) stop("ABM missegregation rate is missing from ", metadata_path, call. = FALSE)
  rate <- as.numeric(candidates[[1L]])[[1L]]
  if (!is.finite(rate) || rate < 0 || rate > 1) {
    stop("Invalid ABM missegregation rate in ", metadata_path, call. = FALSE)
  }
  rate
}

run_alfak_abm <- function(...) {
  fn <- get("run_karyotype_abm", envir = asNamespace("alfakR"), inherits = FALSE)
  fn(...)
}

bind_status <- function(rows) {
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  do.call(rbind, lapply(rows, function(x) {
    x[setdiff(columns, names(x))] <- NA
    x[columns]
  }))
}

read_phase1_parameters <- function(project_dir) {
  path <- file.path(project_dir, "outputs", "bounded_grf", "p_mis_lhs.csv")
  parameters <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("p_index", "p_mis") %in% names(parameters))) stop("Invalid p_mis table: ", path, call. = FALSE)
  parameters[order(parameters$p_index), , drop = FALSE]
}

p_mis_dir <- function(p_index, p_mis) {
  sprintf("p_mis_%02d_%.8f", as.integer(p_index), as.numeric(p_mis))
}

expected_phase2_grid <- function(parameters, landscape_id, n_replicates = 10L) {
  p_dirs <- p_mis_dir(parameters$p_index, parameters$p_mis)
  grid <- expand.grid(
    p1_dir = p_dirs,
    replicate_id = sprintf("replicate_%02d", seq_len(n_replicates)),
    p2_dir = p_dirs,
    stringsAsFactors = FALSE
  )
  grid$landscape_id <- landscape_id
  grid <- grid[, c("landscape_id", "p1_dir", "replicate_id", "p2_dir")]
  grid$relative <- file.path(grid$landscape_id, grid$p1_dir, grid$replicate_id, grid$p2_dir)
  grid
}

stop_if_missing_files <- function(paths, label, max_examples = 10L) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    examples <- paste(utils::head(missing, max_examples), collapse = "\n  ")
    stop(sprintf("Missing %d %s file(s). First missing paths:\n  %s", length(missing), label, examples), call. = FALSE)
  }
  invisible(paths)
}

validate_phase2_grid_rows <- function(data, expected, label) {
  keys <- c("landscape_id", "p1_dir", "replicate_id", "p2_dir")
  missing_columns <- setdiff(keys, names(data))
  if (length(missing_columns)) stop(label, " is missing key columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  observed_key <- do.call(paste, c(data[keys], sep = "\r"))
  expected_key <- do.call(paste, c(expected[keys], sep = "\r"))
  duplicate_key <- unique(observed_key[duplicated(observed_key)])
  missing_key <- setdiff(expected_key, observed_key)
  extra_key <- setdiff(observed_key, expected_key)
  if (length(duplicate_key) || length(missing_key) || length(extra_key)) {
    stop(sprintf(
      "%s does not match the expected phase-2 grid: %d duplicates, %d missing, %d unexpected rows.",
      label, length(duplicate_key), length(missing_key), length(extra_key)
    ), call. = FALSE)
  }
  invisible(TRUE)
}
