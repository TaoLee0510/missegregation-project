#!/usr/bin/env Rscript

# Generate a chromosome-specifically bounded GRF ABM. Initial karyotype
# candidates are sampled from the bounded GRF, then converted to cell counts
# from an x0 frequency vector using the same integer-allocation logic as
# alfakR's ABM wrapper. Every in-bounds descendant receives fitness directly
# from the GRF. Bounds are empirical PDX support, not the small FQ + NN
# inference shell.

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]]) else getwd()
workers <- if (length(args) >= 2L) as.integer(args[[2]]) else max(1L, parallel::detectCores(logical = TRUE) - 1L)
landscape_index <- if (length(args) >= 3L) as.integer(args[[3]]) else NA_integer_
p_index_filter <- if (length(args) >= 5L) as.integer(args[[4]]) else NA_integer_
replicate_filter <- if (length(args) >= 5L) as.integer(args[[5]]) else NA_integer_
profile_enabled <- tolower(Sys.getenv("PROFILE_PHASE1_ABM", unset = "false")) %in% c("1", "true", "yes")
profile_clock <- function() unname(proc.time()[["elapsed"]])
profile_events <- list()
profile_run_id <- Sys.getenv("PHASE1_ABM_PROFILE_ID", unset = "")
if (profile_enabled && !nzchar(profile_run_id)) {
  profile_run_id <- paste(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    if (is.na(landscape_index)) "all_landscapes" else sprintf("landscape_%02d", landscape_index),
    if (is.na(p_index_filter)) "all_p" else sprintf("p_%02d", p_index_filter),
    if (is.na(replicate_filter)) "all_replicates" else sprintf("replicate_%02d", replicate_filter),
    sep = "_"
  )
}
profile_root <- Sys.getenv(
  "PHASE1_ABM_PROFILE_DIR",
  unset = file.path(project_dir, "outputs", "diagnostics", "phase1_abm_profile")
)
profile_output_dir <- if (profile_enabled) file.path(profile_root, profile_run_id) else NA_character_
if (profile_enabled) dir.create(profile_output_dir, recursive = TRUE, showWarnings = FALSE)

profile_step <- function(stage, expr, task = NULL, message = NA_character_, output_path = NA_character_) {
  if (!profile_enabled) return(force(expr))
  started_at <- Sys.time()
  started_elapsed <- profile_clock()
  status <- "completed"
  error_message <- message
  value <- tryCatch(
    force(expr),
    error = function(e) {
      status <<- "failed"
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  ended_at <- Sys.time()
  task <- if (is.null(task)) list() else task
  profile_events[[length(profile_events) + 1L]] <<- data.frame(
    stage = stage,
    status = status,
    elapsed_seconds = profile_clock() - started_elapsed,
    started_at = format(started_at, "%Y-%m-%d %H:%M:%OS3 %Z"),
    ended_at = format(ended_at, "%Y-%m-%d %H:%M:%OS3 %Z"),
    landscape_id = if (!is.null(task$landscape_id)) task$landscape_id else NA_character_,
    p_index = if (!is.null(task$p_index)) task$p_index else NA_integer_,
    p_mis = if (!is.null(task$p_mis)) task$p_mis else NA_real_,
    replicate_id = if (!is.null(task$replicate_id)) task$replicate_id else NA_integer_,
    n_steps = NA_integer_,
    output_path = output_path,
    message = error_message,
    stringsAsFactors = FALSE
  )
  if (identical(status, "failed")) stop(error_message, call. = FALSE)
  value
}

write_profile_outputs <- function(status = NULL, status_path = NA_character_) {
  if (!profile_enabled) return(invisible(NULL))
  rows <- if (length(profile_events)) do.call(rbind, profile_events) else data.frame()
  if (nrow(rows) && exists("n_steps", inherits = TRUE)) rows$n_steps[is.na(rows$n_steps)] <- get("n_steps", inherits = TRUE)
  timing_path <- file.path(profile_output_dir, "timing.csv")
  utils::write.csv(rows, timing_path, row.names = FALSE)

  top_stages <- if (nrow(rows)) {
    totals <- stats::aggregate(elapsed_seconds ~ stage, data = rows, FUN = sum)
    totals <- totals[order(-totals$elapsed_seconds), , drop = FALSE]
    utils::head(totals, 20L)
  } else {
    data.frame(stage = character(), elapsed_seconds = numeric())
  }
  summary_path <- file.path(profile_output_dir, "summary.txt")
  writeLines(c(
    "phase1 bounded-GRF ABM profile",
    paste("project_dir:", project_dir),
    paste("results_dir:", if (exists("results_dir", inherits = TRUE)) get("results_dir", inherits = TRUE) else NA_character_),
    paste("profile_dir:", profile_output_dir),
    paste("status_path:", status_path),
    "",
    "top stages by elapsed seconds:",
    capture.output(print(top_stages, row.names = FALSE)),
    "",
    "run status:",
    if (is.null(status)) "not available" else capture.output(print(status, row.names = FALSE))
  ), summary_path)
  message("Profile timing written to: ", timing_path)
  invisible(timing_path)
}
if (!is.finite(workers) || workers < 1L) stop("`workers` must be a positive integer.", call. = FALSE)
if (profile_enabled && workers > 1L) {
  warning("PROFILE_PHASE1_ABM is enabled; forcing workers = 1 so timing events are ordered.", call. = FALSE)
  workers <- 1L
}
if (!is.na(landscape_index) && (!is.finite(landscape_index) || landscape_index < 1L)) {
  stop("`landscape_index` must be a positive integer.", call. = FALSE)
}
if (!is.na(p_index_filter) && (!is.finite(p_index_filter) || p_index_filter < 1L)) {
  stop("`p_index` must be a positive integer.", call. = FALSE)
}
invisible(profile_step("load_project_helpers", source(file.path(project_dir, "R", "project_helpers.R"))))
invisible(profile_step("load_project_alfak", load_project_alfak(project_dir)))
if (!is.na(replicate_filter) && (!is.finite(replicate_filter) || replicate_filter < 1L)) {
  stop("`replicate_id` must be a positive integer.", call. = FALSE)
}

invisible(profile_step("load_reference_inputs", {
  landscape_dir <- file.path(project_dir, "data", "landscapes")
  reference <- readRDS(file.path(project_dir, "data", "reference_ploidy", "reference_ploidy.rds"))
  reference_digest <- object_digest(reference)
  bound_profile <- reference$bound_profile
  bound_profile_digest <- object_digest(bound_profile)
  lower_copy_numbers <- as.integer(bound_profile$lower_copy_number)
  upper_copy_numbers <- as.integer(bound_profile$upper_copy_number)
  reference_mean_ploidy <- reference$reference_mean_ploidy
  arm_loci <- readRDS(file.path(project_dir, "data", "salehi_reference", "raw", "arm_loci.Rds"))
  length_table <- aggregate(end ~ chrom, data = subset(arm_loci, chrom %in% as.character(1:22)), FUN = max)
  chromosome_lengths <- length_table$end[match(as.character(seq_len(22L)), length_table$chrom)]
}))
results_dir <- Sys.getenv("BOUNDED_GRF_RESULTS_DIR", unset = file.path(project_dir, "outputs", "bounded_grf"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

n_cells <- 10000L
# With per-ABM task args, the sixth CLI argument is an optional short-run
# override.  The older four-argument form still treats the fourth argument as
# the smoke-test step count.
n_steps <- if (length(args) >= 6L) {
  as.integer(args[[6]])
} else if (length(args) == 4L) {
  as.integer(args[[4]])
} else {
  2000L
}
if (!is.finite(n_steps) || n_steps < 0L) stop("`n_steps` must be a non-negative integer.", call. = FALSE)
dt <- 0.1
record_interval <- 50L # Preserve 5-day summaries; integration uses every 0.1-day step.
n_replicates <- 10L
n_p_mis <- 20L
n_frequent_karyotypes <- 50L
diploid_tag <- paste(rep.int(2L, 22L), collapse = ".")
model_version <- "bounded_grf_v7_alfak_x0_initialization_culling_cap_no_diploid_state"
landscape_digest_cache <- new.env(parent = emptyenv())

landscape_digest <- function(landscape_id) {
  if (!exists(landscape_id, envir = landscape_digest_cache, inherits = FALSE)) {
    assign(landscape_id, file_digest(file.path(landscape_dir, paste0(landscape_id, ".rds"))),
           envir = landscape_digest_cache)
  }
  get(landscape_id, envir = landscape_digest_cache, inherits = FALSE)
}

fitness_grf <- function(karyotypes, centroids, lambda) {
  karyotypes <- as.matrix(karyotypes)
  vapply(seq_len(nrow(karyotypes)), function(i) {
    d <- sqrt(rowSums((centroids - rep(karyotypes[i, ], each = nrow(centroids)))^2))
    sum(sin(d / lambda)) / (pi * sqrt(nrow(centroids)))
  }, numeric(1))
}

weighted_ploidy <- function(karyotypes) {
  karyotypes <- as.matrix(karyotypes)
  drop(karyotypes %*% chromosome_lengths / sum(chromosome_lengths))
}

make_reference_ploidy_seed <- function() {
  state <- pmax(lower_copy_numbers, pmin(upper_copy_numbers, rep.int(floor(reference_mean_ploidy), length(lower_copy_numbers))))
  while (weighted_ploidy(matrix(state, nrow = 1L)) < reference_mean_ploidy) {
    candidates <- which(state < upper_copy_numbers)
    if (!length(candidates)) break
    chromosome <- sample(candidates, 1L)
    state[[chromosome]] <- state[[chromosome]] + 1L
  }
  state
}

# Sample distinct, plausible initial FQ states from the bounded GRF.  This is
# a Metropolis chain on one-missegregation moves; the ploidy factor only
# defines the starting x0 ensemble and does not modify the GRF during the ABM.
sample_initial_fq <- function(landscape, n_fq, seed, burnin = 1000L, thin = 20L, ploidy_sd = 0.35, ploidy_band = 0.35) {
  set.seed(seed)
  state <- make_reference_ploidy_seed()
  score <- function(x) {
    if (abs(weighted_ploidy(matrix(x, nrow = 1L)) - reference_mean_ploidy) > ploidy_band) return(0)
    if (identical(paste(x, collapse = "."), diploid_tag)) return(0)
    f <- fitness_grf(matrix(x, nrow = 1L), landscape$centroids, landscape$lambda)[[1L]]
    max(f, .Machine$double.eps) * exp(-0.5 * ((weighted_ploidy(matrix(x, nrow = 1L)) - reference_mean_ploidy) / ploidy_sd)^2)
  }
  current_score <- score(state)
  attempts <- 0L
  accepted <- character()
  accepted_scores <- numeric()
  max_attempts <- burnin + thin * n_fq * 100L
  while (length(accepted) < n_fq && attempts < max_attempts) {
    attempts <- attempts + 1L
    candidate <- state
    chromosome <- sample.int(length(candidate), 1L)
    candidate[[chromosome]] <- candidate[[chromosome]] + sample(c(-1L, 1L), 1L)
    if (candidate[[chromosome]] >= lower_copy_numbers[[chromosome]] && candidate[[chromosome]] <= upper_copy_numbers[[chromosome]]) {
      candidate_score <- score(candidate)
      if (candidate_score > 0 && stats::runif(1) < min(1, candidate_score / current_score)) {
        state <- candidate
        current_score <- candidate_score
      }
    }
    if (attempts > burnin && attempts %% thin == 0L && current_score > 0) {
      tag <- paste(state, collapse = ".")
      if (!(tag %in% accepted)) {
        accepted <- c(accepted, tag)
        accepted_scores <- c(accepted_scores, current_score)
      }
    }
  }
  if (length(accepted) < n_fq) stop("Could not sample the requested number of distinct in-bounds FQ states; increase bounds or sampler budget.", call. = FALSE)
  data.frame(
    karyotype = accepted[seq_len(n_fq)],
    x0_weight = accepted_scores[seq_len(n_fq)],
    stringsAsFactors = FALSE
  )
}

largest_remainder_allocate_alfak <- function(prob, total_size) {
  if (!length(prob)) stop("`prob` must not be empty.", call. = FALSE)
  if (any(!is.finite(prob)) || any(prob < 0)) stop("`prob` must contain finite non-negative values.", call. = FALSE)
  if (sum(prob) <= 0) stop("`prob` must sum to a positive value.", call. = FALSE)
  if (!is.finite(total_size) || total_size < 0 || total_size != floor(total_size)) stop("`total_size` must be a non-negative integer.", call. = FALSE)
  original_names <- names(prob)
  prob <- prob / sum(prob)
  raw <- prob * total_size
  counts <- floor(raw)
  remaining <- total_size - sum(counts)
  if (!is.finite(remaining) || remaining < 0 || remaining != floor(remaining) || remaining > length(prob)) {
    stop("Internal error: largest remainder allocation produced an invalid remainder.", call. = FALSE)
  }
  if (remaining > 0) {
    fractional <- raw - counts
    order_idx <- order(-fractional, seq_along(fractional))
    counts[order_idx[seq_len(remaining)]] <- counts[order_idx[seq_len(remaining)]] + 1
  }
  if (length(counts) != length(prob) || any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts)) || sum(counts) != total_size) {
    stop("Internal error: largest remainder allocation produced invalid integer-valued counts.", call. = FALSE)
  }
  if (!is.null(original_names)) names(counts) <- original_names
  counts
}

prepare_initial_population_from_x0 <- function(x0, total_size) {
  if (is.null(names(x0)) || any(!nzchar(names(x0)))) stop("`x0` must be a named numeric vector.", call. = FALSE)
  if (anyDuplicated(names(x0))) stop("`x0` must not contain duplicate karyotypes.", call. = FALSE)
  initial_counts <- largest_remainder_allocate_alfak(x0, total_size)
  positive <- initial_counts > 0
  if (!any(positive)) stop("Initial population for ABM is zero after filtering zero counts.", call. = FALSE)
  positive_counts <- as.integer(initial_counts[positive])
  names(positive_counts) <- names(initial_counts)[positive]
  allocated_counts <- as.integer(initial_counts)
  names(allocated_counts) <- names(initial_counts)
  list(
    tags = names(initial_counts)[positive],
    counts = positive_counts,
    x0 = x0 / sum(x0),
    allocated_counts = allocated_counts
  )
}

estimate_peak_reference <- function(landscape, seed, n_draws = 100000L) {
  set.seed(seed)
  draws <- sapply(seq_len(landscape$n_chromosomes), function(i) {
    sample.int(upper_copy_numbers[[i]] - lower_copy_numbers[[i]] + 1L, n_draws, replace = TRUE) + lower_copy_numbers[[i]] - 1L
  })
  max(fitness_grf(draws, landscape$centroids, landscape$lambda))
}

state_summary <- function(counts, landscape, peak_threshold) {
  if (!length(counts)) {
    return(data.frame(population = 0, diversity = 0, peak_cells = 0, peak_percent = NA_real_))
  }
  k <- do.call(rbind, strsplit(names(counts), ".", fixed = TRUE))
  storage.mode(k) <- "numeric"
  fit <- fitness_grf(k, landscape$centroids, landscape$lambda)
  population <- sum(counts)
  peak_cells <- sum(counts[fit >= peak_threshold])
  data.frame(population = population, diversity = sum(counts > 0), peak_cells = peak_cells,
             peak_percent = 100 * peak_cells / population)
}

write_initial_population <- function(landscape, results_landscape_dir) {
  profile_task <- list(landscape_id = landscape$landscape_id)
  initial_path <- file.path(results_landscape_dir, "initial_population.csv")
  candidate_path <- file.path(results_landscape_dir, "initial_x0_candidates.csv")
  meta_path <- file.path(results_landscape_dir, "initialization.rds")
  landscape_number <- as.integer(sub(".*_", "", landscape$landscape_id))
  seed <- 710000L + landscape_number
  expected_provenance <- list(
    model_version = model_version,
    landscape_id = landscape$landscape_id,
    landscape_digest = landscape_digest(landscape$landscape_id),
    reference_digest = reference_digest,
    bound_profile_digest = bound_profile_digest,
    seed = seed,
    n_cells = n_cells,
    n_frequent_karyotypes = n_frequent_karyotypes,
    initial_population_method = "alfakR_prepare_abm_initial_population_largest_remainder",
    x0_weight_method = "bounded-GRF Metropolis score: max(GRF fitness, eps) times Gaussian ploidy weight",
    allocation_method = "alfakR_largest_remainder_allocate",
    target_ploidy = reference_mean_ploidy,
    ploidy_band = 0.35
  )
  if (file.exists(initial_path) && file.exists(candidate_path) && file.exists(meta_path)) {
    existing <- readRDS(meta_path)
    if (provenance_matches(existing, expected_provenance)) return(invisible(NULL))
  }
  candidates <- profile_step(
    "initial_population/sample_initial_fq",
    sample_initial_fq(landscape, n_frequent_karyotypes, seed = seed),
    task = profile_task
  )
  x0 <- stats::setNames(candidates$x0_weight, candidates$karyotype)
  allocated <- profile_step(
    "initial_population/allocate_counts",
    prepare_initial_population_from_x0(x0, n_cells),
    task = profile_task
  )
  counts <- allocated$counts
  fq_tags <- allocated$tags
  names(counts) <- fq_tags
  initial <- list(
    counts = unname(counts),
    tags = fq_tags,
    x0 = allocated$x0,
    x0_candidate_counts = allocated$allocated_counts,
    sampler = list(method = "bounded-GRF Metropolis x0 candidate sample", n_fq = n_frequent_karyotypes,
                   x0_weight_method = "bounded-GRF score times Gaussian ploidy weight",
                   allocation_method = "alfakR largest_remainder_allocate with zero-count filtering",
                   target_ploidy = reference_mean_ploidy, ploidy_band = 0.35, model_version = model_version),
    provenance = expected_provenance
  )
  candidate_k <- do.call(rbind, strsplit(candidates$karyotype, ".", fixed = TRUE))
  storage.mode(candidate_k) <- "numeric"
  candidate_fitness <- profile_step(
    "initial_population/score_candidates",
    fitness_grf(candidate_k, landscape$centroids, landscape$lambda),
    task = profile_task
  )
  candidate_counts <- allocated$allocated_counts[candidates$karyotype]
  profile_step(
    "initial_population/write_candidate_table",
    write_csv_atomic(data.frame(karyotype = candidates$karyotype, x0_weight = candidates$x0_weight,
                                x0_frequency = allocated$x0[candidates$karyotype],
                                allocated_count = as.integer(candidate_counts),
                                fitness = candidate_fitness,
                                weighted_ploidy = weighted_ploidy(candidate_k)),
                     candidate_path, row.names = FALSE),
    task = profile_task,
    output_path = candidate_path
  )
  k <- do.call(rbind, strsplit(fq_tags, ".", fixed = TRUE))
  storage.mode(k) <- "numeric"
  fitness <- profile_step(
    "initial_population/score_allocated_population",
    fitness_grf(k, landscape$centroids, landscape$lambda),
    task = profile_task
  )
  profile_step(
    "initial_population/write_initial_table",
    write_csv_atomic(data.frame(karyotype = initial$tags, count = initial$counts, fitness = fitness,
                                x0_frequency = allocated$x0[fq_tags],
                                weighted_ploidy = weighted_ploidy(k)), initial_path, row.names = FALSE),
    task = profile_task,
    output_path = initial_path
  )
  profile_step(
    "initial_population/write_metadata",
    write_rds_atomic(initial, meta_path),
    task = profile_task,
    output_path = meta_path
  )
}

prepare_landscape <- function(landscape_id) {
  profile_task <- list(landscape_id = landscape_id)
  landscape <- profile_step(
    "prepare_landscape/read_landscape",
    readRDS(file.path(landscape_dir, paste0(landscape_id, ".rds"))),
    task = profile_task
  )
  landscape_result_dir <- file.path(results_dir, landscape_id)
  dir.create(landscape_result_dir, recursive = TRUE, showWarnings = FALSE)
  lock_dir <- file.path(landscape_result_dir, ".prepare_landscape.lock")
  if (!dir.create(lock_dir, showWarnings = FALSE)) {
    for (attempt in seq_len(1800L)) {
      if (!dir.exists(lock_dir)) return(prepare_landscape(landscape_id))
      Sys.sleep(2)
    }
    stop("Timed out waiting for landscape preparation lock: ", lock_dir, call. = FALSE)
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  profile_step(
    "prepare_landscape/write_initial_population",
    write_initial_population(landscape, landscape_result_dir),
    task = profile_task
  )
  peak_path <- file.path(landscape_result_dir, "peak_reference.rds")
  peak_seed <- 810000L + as.integer(sub(".*_", "", landscape_id))
  peak_n_draws <- 100000L
  peak_provenance <- list(
    model_version = "peak_reference_v1",
    landscape_id = landscape_id,
    landscape_digest = landscape_digest(landscape_id),
    reference_digest = reference_digest,
    bound_profile_digest = bound_profile_digest,
    seed = peak_seed,
    n_uniform_draws = peak_n_draws,
    threshold_fraction = 0.95
  )
  if (file.exists(peak_path)) {
    existing <- readRDS(peak_path)
    if (provenance_matches(existing, peak_provenance)) return(invisible(NULL))
  }
  peak_max <- profile_step(
    "peak_reference/estimate_domain_max",
    estimate_peak_reference(landscape, seed = peak_seed, n_draws = peak_n_draws),
    task = profile_task
  )
  profile_step(
    "peak_reference/write_metadata",
    write_rds_atomic(list(estimated_domain_max_fitness = peak_max, peak_threshold = 0.95 * peak_max,
                          domain = "empirical chromosome-specific bounded GRF", n_uniform_draws = peak_n_draws,
                          provenance = peak_provenance), peak_path),
    task = profile_task,
    output_path = peak_path
  )
}

manifest <- profile_step("load_landscape_manifest", read.csv(file.path(landscape_dir, "manifest.csv"), stringsAsFactors = FALSE))
if (!is.na(landscape_index)) {
  if (landscape_index > nrow(manifest)) stop("`landscape_index` exceeds the number of landscapes.", call. = FALSE)
  manifest <- manifest[landscape_index, , drop = FALSE]
}
set.seed(20260713)
# Endpoint-inclusive one-dimensional Latin hypercube, with the literature-range
# endpoints explicitly retained and one random point in each interior stratum.
lhs_u <- c(0, (seq_len(n_p_mis - 2L) + stats::runif(n_p_mis - 2L)) / n_p_mis, 1)
p_mis <- 0.00025 + lhs_u * (0.01 - 0.00025)
param_table <- data.frame(p_index = seq_along(p_mis), p_mis = p_mis)
param_path <- file.path(results_dir, "p_mis_lhs.csv")
invisible(profile_step("prepare_p_mis_table", {
  if (!file.exists(param_path)) {
    param_tmp <- tempfile("p_mis_lhs_", tmpdir = results_dir)
    utils::write.csv(param_table, param_tmp, row.names = FALSE)
    if (!file.rename(param_tmp, param_path) && !file.exists(param_path)) {
      stop("Could not create p_mis_lhs.csv.", call. = FALSE)
    }
  } else {
    existing_param <- read.csv(param_path, stringsAsFactors = FALSE)
    invalid_param <- !all(names(param_table) %in% names(existing_param)) ||
      nrow(existing_param) != nrow(param_table) ||
      !identical(as.integer(existing_param$p_index), as.integer(param_table$p_index)) ||
      !isTRUE(all.equal(existing_param$p_mis, param_table$p_mis, tolerance = 1e-12, check.attributes = FALSE))
    if (invalid_param) {
      stop("Existing p_mis_lhs.csv does not match the current deterministic p_mis table.", call. = FALSE)
    }
  }
}, output_path = param_path))
p_mis_lhs_digest <- file_digest(param_path)
if (!is.na(p_index_filter) && p_index_filter > nrow(param_table)) stop("`p_index` exceeds the p_mis table.", call. = FALSE)
if (!is.na(replicate_filter) && replicate_filter > n_replicates) stop("`replicate_id` exceeds the replicate count.", call. = FALSE)
p_indices <- if (is.na(p_index_filter)) seq_along(p_mis) else p_index_filter
replicate_ids <- if (is.na(replicate_filter)) seq_len(n_replicates) else replicate_filter

tasks <- profile_step("build_task_grid", {
  tasks <- list()
  for (landscape_id in manifest$landscape_id) {
    for (p_index in p_indices) {
      for (replicate_id in replicate_ids) {
        tasks[[length(tasks) + 1L]] <- list(landscape_id = landscape_id, p_index = p_index, p_mis = p_mis[[p_index]], replicate_id = replicate_id)
      }
    }
  }
  tasks
})
for (landscape_id in unique(vapply(tasks, `[[`, character(1), "landscape_id"))) {
  invisible(profile_step("prepare_landscape/total", prepare_landscape(landscape_id), task = list(landscape_id = landscape_id)))
}

run_task <- function(task) {
  landscape <- profile_step(
    "task/read_landscape",
    readRDS(file.path(landscape_dir, paste0(task$landscape_id, ".rds"))),
    task = task
  )
  landscape_result_dir <- file.path(results_dir, task$landscape_id)
  dir.create(landscape_result_dir, recursive = TRUE, showWarnings = FALSE)
  initial <- profile_step(
    "task/read_initial_population",
    readRDS(file.path(landscape_result_dir, "initialization.rds")),
    task = task
  )
  initial_counts <- initial$counts
  names(initial_counts) <- initial$tags
  p_dir <- file.path(landscape_result_dir, sprintf("p_mis_%02d_%.8f", task$p_index, task$p_mis))
  replicate_dir <- file.path(p_dir, sprintf("replicate_%02d", task$replicate_id))
  dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)
  final_path <- file.path(replicate_dir, "final_karyotypes.csv")
  summary_path <- file.path(replicate_dir, "trajectory_summary.csv")
  observations_path <- file.path(replicate_dir, "abm_observations.rds")
  metadata_path <- file.path(replicate_dir, "run_metadata.rds")
  initial_path <- file.path(landscape_result_dir, "initialization.rds")
  expected_provenance <- profile_step(
    "task/build_expected_provenance",
    list(
      model_version = model_version,
      task = task,
      landscape_digest = landscape_digest(task$landscape_id),
      reference_digest = reference_digest,
      bound_profile_digest = bound_profile_digest,
      p_mis_lhs_digest = p_mis_lhs_digest,
      initialization_digest = file_digest(initial_path),
      n_cells = n_cells,
      n_steps = n_steps,
      dt = dt,
      record_interval = record_interval,
      culling_survival_fraction = 0.999,
      excluded_karyotypes = diploid_tag
    ),
    task = task
  )
  if (file.exists(final_path) && file.exists(summary_path) && file.exists(observations_path) && file.exists(metadata_path)) {
    previous <- readRDS(metadata_path)
    if (provenance_matches(previous, expected_provenance)) return(data.frame(status = "skipped", task))
  }

  peak_reference <- profile_step(
    "task/read_peak_reference",
    readRDS(file.path(landscape_result_dir, "peak_reference.rds")),
    task = task
  )
  peak_max <- peak_reference$estimated_domain_max_fitness
  peak_threshold <- peak_reference$peak_threshold
  seed <- 910000L + as.integer(sub(".*_", "", task$landscape_id)) * 1000L + task$p_index * 10L + task$replicate_id
  raw <- profile_step(
    "task/run_alfak_abm",
    run_alfak_abm(
      initial_population_r = as.list(initial_counts),
      fitness_map_r = stats::setNames(list(), character(0)),
      p_missegregation = task$p_mis,
      dt = dt,
      n_steps = n_steps,
      max_population_size = n_cells,
      culling_survival_fraction = 0.999,
      record_interval = record_interval,
      seed = seed,
      grf_centroids = landscape$centroids,
      grf_lambda = landscape$lambda,
      lower_copy_numbers = lower_copy_numbers,
      upper_copy_numbers = upper_copy_numbers,
      excluded_karyotypes = diploid_tag
    ),
    task = task
  )
  boundary_rejections <- unname(attr(raw, "boundary_rejections")[[1L]])
  observations <- profile_step("task/raw_to_observations", raw_to_observations(raw, dt), task = task)
  profile_step("task/write_observations", saveRDS(observations, observations_path), task = task, output_path = observations_path)
  summaries <- profile_step(
    "task/build_trajectory_summary",
    {
      summaries <- lapply(names(raw), function(step) {
        counts <- raw[[step]]
        s <- state_summary(counts, landscape, peak_threshold)
        cbind(time_days = as.numeric(step) * dt, step = as.integer(step), s)
      })
      do.call(rbind, summaries)
    },
    task = task
  )
  final_counts <- raw[[tail(names(raw), 1L)]]
  k <- do.call(rbind, strsplit(names(final_counts), ".", fixed = TRUE))
  storage.mode(k) <- "numeric"
  final_fitness <- profile_step("task/score_final_population", fitness_grf(k, landscape$centroids, landscape$lambda), task = task)
  profile_step(
    "task/write_final_karyotypes",
    utils::write.csv(data.frame(karyotype = names(final_counts), count = as.numeric(final_counts), fitness = final_fitness,
                                is_peak = final_fitness >= peak_threshold), final_path, row.names = FALSE),
    task = task,
    output_path = final_path
  )
  profile_step("task/write_trajectory_summary", utils::write.csv(summaries, summary_path, row.names = FALSE), task = task, output_path = summary_path)
  profile_step(
    "task/write_run_metadata",
    saveRDS(list(task = task, seed = seed, peak_reference_max = peak_max, peak_threshold = peak_threshold,
                 integration_dt_days = dt, n_steps = n_steps, record_interval_steps = record_interval,
                 carrying_capacity = n_cells, culling_survival_fraction = 0.999,
                 fitness_mode = "full GRF evaluated on every in-bounds karyotype",
                 excluded_karyotypes = diploid_tag,
                 bounds = bound_profile, reference_mean_ploidy = reference_mean_ploidy,
                 boundary_rejections = boundary_rejections, model_version = model_version,
                 provenance = expected_provenance),
            metadata_path),
    task = task,
    output_path = metadata_path
  )
  data.frame(status = "completed", task,
             final_population = tail(summaries$population, 1L),
             final_diversity = tail(summaries$diversity, 1L),
             final_peak_percent = tail(summaries$peak_percent, 1L))
}

run_task_safely <- function(task) {
  tryCatch(run_task(task), error = function(e) data.frame(
    status = "failed", landscape_id = task$landscape_id, p_index = task$p_index,
    p_mis = task$p_mis, replicate_id = task$replicate_id, message = conditionMessage(e)
  ))
}
status <- if (.Platform$OS.type == "windows" || workers <= 1L) {
  lapply(tasks, run_task_safely)
} else {
  parallel::mclapply(tasks, run_task_safely,
                     mc.cores = workers, mc.preschedule = FALSE)
}
status <- bind_status(status)
status_path <- if (is.na(landscape_index) && is.na(p_index_filter) && is.na(replicate_filter)) {
  file.path(results_dir, "run_status.csv")
} else {
  status_suffix <- paste(
    if (is.na(landscape_index)) "all_landscapes" else sprintf("landscape_%02d", landscape_index),
    if (is.na(p_index_filter)) "all_p" else sprintf("p_%02d", p_index_filter),
    if (is.na(replicate_filter)) "all_replicates" else sprintf("replicate_%02d", replicate_filter),
    sep = "_"
  )
  file.path(results_dir, paste0("run_status_", status_suffix, ".csv"))
}
invisible(profile_step("write_run_status", utils::write.csv(status, status_path, row.names = FALSE), output_path = status_path))
write_profile_outputs(status, status_path)
if (any(status$status == "failed")) quit(status = 1L)
