#!/usr/bin/env Rscript

# Generate the GRF landscapes used for the CIN study. Parameterization follows
# the ALFA-K source accompanying s41467-025-67750-0:
# 22 chromosome dimensions, 10 centroids, centroid coordinates in -10:20,
# and wavelength range 0.2--1.6.

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]]) else getwd()
out_dir <- file.path(project_dir, "data", "landscapes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

arg_value <- function(args, name) {
  equals_prefix <- paste0("--", name, "=")
  equals_match <- args[startsWith(args, equals_prefix)]
  if (length(equals_match)) return(sub(equals_prefix, "", equals_match[[1L]], fixed = TRUE))
  flag_index <- match(paste0("--", name), args)
  if (!is.na(flag_index) && flag_index < length(args)) return(args[[flag_index + 1L]])
  NULL
}

parse_positive_integer <- function(value, name) {
  out <- suppressWarnings(as.integer(value))
  if (!is.finite(out) || out < 1L || out != as.numeric(value)) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  out
}

set.seed(20260712)
n_landscapes_arg <- arg_value(args, "n-landscapes")
if (is.null(n_landscapes_arg) && length(args) >= 2L && !startsWith(args[[2L]], "--")) {
  n_landscapes_arg <- args[[2L]]
}
n_landscapes <- if (is.null(n_landscapes_arg)) 200L else parse_positive_integer(n_landscapes_arg, "n_landscapes")
n_chromosomes <- 22L
n_centroids <- 10L
lambda_min <- 0.2
lambda_max <- 1.6
centroid_coordinate_min <- -10L
centroid_coordinate_max <- 20L
base_seed <- 20260712L

# One-dimensional Latin hypercube: exactly one lambda from each stratum.
lambda <- ((seq_len(n_landscapes) - stats::runif(n_landscapes)) / n_landscapes)
lambda <- sample(lambda_min + lambda * (lambda_max - lambda_min))

manifest <- vector("list", n_landscapes)
generation_parameters <- vector("list", n_landscapes)
centroid_column_names <- as.vector(t(outer(
  sprintf("centroid_%02d", seq_len(n_centroids)),
  sprintf("chr_%02d", seq_len(n_chromosomes)),
  paste,
  sep = "_"
)))
for (i in seq_len(n_landscapes)) {
  landscape_seed <- base_seed + i
  set.seed(landscape_seed)
  centroids <- matrix(
    sample(centroid_coordinate_min:centroid_coordinate_max, n_centroids * n_chromosomes, replace = TRUE),
    nrow = n_centroids,
    ncol = n_chromosomes
  )
  landscape_id <- sprintf("landscape_%02d", i)
  landscape_file <- paste0(landscape_id, ".rds")
  centroid_file <- paste0(landscape_id, "_centroids.csv")
  landscape <- list(
    landscape_id = landscape_id,
    centroids = centroids,
    lambda = lambda[[i]],
    n_centroids = n_centroids,
    n_chromosomes = n_chromosomes,
    centroid_coordinate_range = c(centroid_coordinate_min, centroid_coordinate_max),
    generation_seed = landscape_seed,
    n_landscapes = n_landscapes,
    lambda_range = c(lambda_min, lambda_max),
    lambda_sampling = "one-dimensional Latin hypercube with one value per requested landscape",
    centroid_sampling = "independent integer uniform draws over centroid_coordinate_range",
    source = "ALFA-K supplementary source: scripts/S01_run_abm_sims.R and R/utils_karyo.R"
  )
  landscape_path <- file.path(out_dir, landscape_file)
  centroid_path <- file.path(out_dir, centroid_file)
  saveRDS(landscape, landscape_path)
  utils::write.csv(
    data.frame(centroid_id = seq_len(n_centroids), centroids, check.names = FALSE),
    centroid_path,
    row.names = FALSE
  )
  landscape_digest <- unname(tools::md5sum(landscape_path))
  centroid_digest <- unname(tools::md5sum(centroid_path))
  manifest[[i]] <- data.frame(
    landscape_id = landscape_id,
    lambda = lambda[[i]],
    n_centroids = n_centroids,
    n_chromosomes = n_chromosomes,
    centroid_min = min(centroids),
    centroid_max = max(centroids),
    generation_seed = landscape_seed,
    landscape_file = landscape_file,
    centroid_file = centroid_file,
    landscape_digest = landscape_digest,
    centroid_digest = centroid_digest,
    stringsAsFactors = FALSE
  )
  parameter_row <- data.frame(
    landscape_id = landscape_id,
    landscape_index = i,
    n_landscapes = n_landscapes,
    base_seed = base_seed,
    generation_seed = landscape_seed,
    lambda = lambda[[i]],
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    lambda_sampling = "one-dimensional Latin hypercube",
    n_centroids = n_centroids,
    n_chromosomes = n_chromosomes,
    centroid_coordinate_min = centroid_coordinate_min,
    centroid_coordinate_max = centroid_coordinate_max,
    centroid_sampling = "integer uniform with replacement",
    landscape_file = landscape_file,
    centroid_file = centroid_file,
    landscape_digest = landscape_digest,
    centroid_digest = centroid_digest,
    stringsAsFactors = FALSE
  )
  parameter_row[centroid_column_names] <- as.list(as.integer(as.vector(t(centroids))))
  generation_parameters[[i]] <- parameter_row
}

manifest <- do.call(rbind, manifest)
generation_parameters <- do.call(rbind, generation_parameters)
utils::write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)
saveRDS(manifest, file.path(out_dir, "manifest.rds"))
utils::write.csv(generation_parameters, file.path(out_dir, "grf_generation_parameters.csv"), row.names = FALSE)
saveRDS(generation_parameters, file.path(out_dir, "grf_generation_parameters.rds"))
writeLines(c(
  "# Synthetic GRF landscapes",
  "",
  sprintf("%d landscapes generated by Latin-hypercube sampling of lambda on [0.2, 1.6].", n_landscapes),
  "All other GRF settings are fixed to the authors' synthetic-validation design:",
  "22 chromosome dimensions, 10 centroids, and centroids sampled from integers -10 through 20.",
  "grf_generation_parameters.csv records the exact lambda, seed, fixed settings, file digests, and centroid coordinates for every landscape.",
  "Each RDS contains the complete reproducible landscape definition."
), file.path(out_dir, "README.md"))
