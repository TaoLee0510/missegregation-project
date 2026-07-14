#!/usr/bin/env Rscript

# Generate a chromosome-specifically bounded GRF ABM. The initial population
# is dispersed across 50 in-bounds FQ states; every in-bounds descendant
# receives fitness directly from the GRF. Bounds are empirical PDX support,
# not the small FQ + NN inference shell.

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]]) else getwd()
workers <- if (length(args) >= 2L) as.integer(args[[2]]) else max(1L, parallel::detectCores(logical = TRUE) - 1L)
landscape_index <- if (length(args) >= 3L) as.integer(args[[3]]) else NA_integer_
if (!is.finite(workers) || workers < 1L) stop("`workers` must be a positive integer.", call. = FALSE)
if (!is.na(landscape_index) && (!is.finite(landscape_index) || landscape_index < 1L)) {
  stop("`landscape_index` must be a positive integer.", call. = FALSE)
}
library(alfakR)
source(file.path(project_dir, "R", "project_helpers.R"))

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
results_dir <- file.path(project_dir, "outputs", "bounded_grf")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

n_cells <- 10000L
# Fourth CLI argument is an optional short-run override for smoke tests.
n_steps <- if (length(args) >= 4L) as.integer(args[[4]]) else 2000L
dt <- 0.1
record_interval <- 50L # Preserve 5-day summaries; integration uses every 0.1-day step.
n_replicates <- 10L
n_p_mis <- 20L
n_frequent_karyotypes <- 50L
diploid_tag <- paste(rep.int(2L, 22L), collapse = ".")
model_version <- "bounded_grf_v5_no_diploid_state"
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

# Sample 50 distinct, plausible initial FQ states from the bounded GRF.  This
# is a Metropolis chain on one-missegregation moves; the ploidy factor only
# defines the starting ensemble and does not modify the GRF during the ABM.
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
    if (attempts > burnin && attempts %% thin == 0L) accepted <- unique(c(accepted, paste(state, collapse = ".")))
  }
  if (length(accepted) < n_fq) stop("Could not sample 50 distinct in-bounds FQ states; increase bounds or sampler budget.", call. = FALSE)
  accepted[seq_len(n_fq)]
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
  initial_path <- file.path(results_landscape_dir, "initial_population.csv")
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
    target_ploidy = reference_mean_ploidy,
    ploidy_band = 0.35
  )
  if (file.exists(initial_path) && file.exists(meta_path)) {
    existing <- readRDS(meta_path)
    if (provenance_matches(existing, expected_provenance)) return(invisible(NULL))
  }
  fq_tags <- sample_initial_fq(landscape, n_frequent_karyotypes, seed = seed)
  counts <- rep.int(n_cells %/% n_frequent_karyotypes, n_frequent_karyotypes)
  counts[[1L]] <- counts[[1L]] + n_cells - sum(counts)
  initial <- list(
    counts = counts,
    tags = fq_tags,
    sampler = list(method = "bounded-GRF Metropolis sample dispersed across 50 FQ states", n_fq = n_frequent_karyotypes,
                   target_ploidy = reference_mean_ploidy, ploidy_band = 0.35, model_version = model_version),
    provenance = expected_provenance
  )
  k <- do.call(rbind, strsplit(fq_tags, ".", fixed = TRUE))
  storage.mode(k) <- "numeric"
  fitness <- fitness_grf(k, landscape$centroids, landscape$lambda)
  utils::write.csv(data.frame(karyotype = initial$tags, count = initial$counts, fitness = fitness,
                              weighted_ploidy = weighted_ploidy(k)), initial_path, row.names = FALSE)
  saveRDS(initial, meta_path)
}

prepare_landscape <- function(landscape_id) {
  landscape <- readRDS(file.path(landscape_dir, paste0(landscape_id, ".rds")))
  landscape_result_dir <- file.path(results_dir, landscape_id)
  dir.create(landscape_result_dir, recursive = TRUE, showWarnings = FALSE)
  write_initial_population(landscape, landscape_result_dir)
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
  peak_max <- estimate_peak_reference(landscape, seed = peak_seed, n_draws = peak_n_draws)
  saveRDS(list(estimated_domain_max_fitness = peak_max, peak_threshold = 0.95 * peak_max,
               domain = "empirical chromosome-specific bounded GRF", n_uniform_draws = peak_n_draws,
               provenance = peak_provenance), peak_path)
}

manifest <- read.csv(file.path(landscape_dir, "manifest.csv"), stringsAsFactors = FALSE)
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
p_mis_lhs_digest <- file_digest(param_path)

tasks <- list()
for (landscape_id in manifest$landscape_id) {
  for (p_index in seq_along(p_mis)) {
    for (replicate_id in seq_len(n_replicates)) {
      tasks[[length(tasks) + 1L]] <- list(landscape_id = landscape_id, p_index = p_index, p_mis = p_mis[[p_index]], replicate_id = replicate_id)
    }
  }
}
for (landscape_id in unique(vapply(tasks, `[[`, character(1), "landscape_id"))) prepare_landscape(landscape_id)

run_task <- function(task) {
  landscape <- readRDS(file.path(landscape_dir, paste0(task$landscape_id, ".rds")))
  landscape_result_dir <- file.path(results_dir, task$landscape_id)
  dir.create(landscape_result_dir, recursive = TRUE, showWarnings = FALSE)
  initial <- readRDS(file.path(landscape_result_dir, "initialization.rds"))
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
  expected_provenance <- list(
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
  )
  if (file.exists(final_path) && file.exists(summary_path) && file.exists(observations_path) && file.exists(metadata_path)) {
    previous <- readRDS(metadata_path)
    if (provenance_matches(previous, expected_provenance)) return(data.frame(status = "skipped", task))
  }

  peak_reference <- readRDS(file.path(landscape_result_dir, "peak_reference.rds"))
  peak_max <- peak_reference$estimated_domain_max_fitness
  peak_threshold <- peak_reference$peak_threshold
  seed <- 910000L + as.integer(sub(".*_", "", task$landscape_id)) * 1000L + task$p_index * 10L + task$replicate_id
  raw <- alfakR::run_karyotype_abm(
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
  )
  boundary_rejections <- unname(attr(raw, "boundary_rejections")[[1L]])
  saveRDS(raw_to_observations(raw, dt), observations_path)
  summaries <- lapply(names(raw), function(step) {
    counts <- raw[[step]]
    s <- state_summary(counts, landscape, peak_threshold)
    cbind(time_days = as.numeric(step) * dt, step = as.integer(step), s)
  })
  summaries <- do.call(rbind, summaries)
  final_counts <- raw[[tail(names(raw), 1L)]]
  k <- do.call(rbind, strsplit(names(final_counts), ".", fixed = TRUE))
  storage.mode(k) <- "numeric"
  final_fitness <- fitness_grf(k, landscape$centroids, landscape$lambda)
  utils::write.csv(data.frame(karyotype = names(final_counts), count = as.numeric(final_counts), fitness = final_fitness,
                              is_peak = final_fitness >= peak_threshold), final_path, row.names = FALSE)
  utils::write.csv(summaries, summary_path, row.names = FALSE)
  saveRDS(list(task = task, seed = seed, peak_reference_max = peak_max, peak_threshold = peak_threshold,
               integration_dt_days = dt, n_steps = n_steps, record_interval_steps = record_interval,
               carrying_capacity = n_cells, culling_survival_fraction = 0.999,
               fitness_mode = "full GRF evaluated on every in-bounds karyotype",
               excluded_karyotypes = diploid_tag,
               bounds = bound_profile, reference_mean_ploidy = reference_mean_ploidy,
               boundary_rejections = boundary_rejections, model_version = model_version,
               provenance = expected_provenance),
          metadata_path)
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
status_path <- if (is.na(landscape_index)) {
  file.path(results_dir, "run_status.csv")
} else {
  file.path(results_dir, sprintf("run_status_landscape_%02d.csv", landscape_index))
}
utils::write.csv(status, status_path, row.names = FALSE)
if (any(status$status == "failed")) quit(status = 1L)
