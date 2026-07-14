#!/usr/bin/env Rscript

# Run once on an HPC login/build node after loading the desired R module.
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]]) else getwd()
library_dir <- Sys.getenv("R_LIBS_USER")
if (!nzchar(library_dir)) stop("Set R_LIBS_USER to a writable persistent library directory.", call. = FALSE)

dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))
repos <- c(CRAN = Sys.getenv("R_REPOSITORY", "https://cloud.r-project.org"))
required <- c("Rcpp", "deSolve", "fields", "quadprog", "RSpectra", "tidyr", "igraph", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = library_dir, repos = repos, dependencies = TRUE)
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
  stop("One or more required R packages could not be installed.", call. = FALSE)
}

status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--preclean", paste0("--library=", library_dir),
    file.path(project_dir, "packages", "alfakR"))
)
if (status != 0L || !requireNamespace("alfakR", quietly = TRUE)) {
  stop("Local alfakR installation failed.", call. = FALSE)
}
message("HPC R environment is ready in: ", library_dir)
