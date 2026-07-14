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

prepare_terminal_input <- function(observation_path) {
  observed <- readRDS(observation_path)
  final_path <- file.path(dirname(observation_path), "final_karyotypes.csv")
  final <- read.csv(final_path, stringsAsFactors = FALSE)
  if (!all(c("karyotype", "count") %in% names(final))) stop("Invalid final population for ", observation_path, call. = FALSE)
  terminal <- unique(final$karyotype[final$count > 0])
  if (!length(terminal)) stop("No nonzero endpoint karyotypes for ", observation_path, call. = FALSE)
  missing <- setdiff(terminal, rownames(observed$x))
  if (length(missing)) stop("Endpoint karyotypes are absent from the recorded trajectory: ", paste(missing, collapse = ", "), call. = FALSE)
  list(yi = list(x = observed$x[terminal, , drop = FALSE], dt = observed$dt), passage_times = observed$passage_times, terminal = terminal)
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
