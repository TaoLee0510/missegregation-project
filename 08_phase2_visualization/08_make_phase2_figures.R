#!/usr/bin/env Rscript

# Visualize original-GRF recovery, phase-1/phase-2 fidelity, and support change.
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]]) else getwd()
topology_dir <- file.path(project_dir, "outputs", "phase2_topology")
out_dir <- file.path(project_dir, "outputs", "phase2_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)

read_shards <- function(pattern) {
  paths <- sort(Sys.glob(file.path(topology_dir, pattern)))
  if (length(paths) != 200L) stop("Expected 200 topology shards matching ", pattern, "; found ", length(paths), call. = FALSE)
  do.call(rbind, lapply(paths, read.csv, stringsAsFactors = FALSE))
}
rate <- function(x) factor(sub("^p_mis_[0-9]+_", "", x), levels = sort(unique(sub("^p_mis_[0-9]+_", "", x))))
metrics <- read_shards("phase2_topology_metrics_landscape_*.csv")
terminal <- read_shards("phase2_terminal_fitness_landscape_*.csv")
metrics$p_mis_phase1 <- rate(metrics$p1_dir); metrics$p_mis_phase2 <- rate(metrics$p2_dir)
terminal$p_mis_phase1 <- rate(terminal$p1_dir); terminal$p_mis_phase2 <- rate(terminal$p2_dir)

p_terminal <- ggplot2::ggplot(terminal, ggplot2::aes(x = p_mis_phase2, y = mean_original_grf_fitness, fill = p_mis_phase2)) +
  ggplot2::geom_boxplot(outlier.size = 0.3) + ggplot2::facet_wrap(~p_mis_phase1, ncol = 4) +
  ggplot2::labs(x = expression(p[mis][2]), y = "Terminal cell-weighted original GRF fitness") + ggplot2::theme_classic() + ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(out_dir, "phase2_terminal_original_grf_fitness_by_rate_pair.png"), p_terminal, width = 14, height = 12, dpi = 300)

p_fidelity <- ggplot2::ggplot(metrics, ggplot2::aes(x = p_mis_phase1, y = p_mis_phase2, fill = phase1_phase2_shared_spearman)) +
  ggplot2::stat_summary(fun = mean, geom = "tile") + ggplot2::scale_fill_viridis_c(na.value = "grey90") +
  ggplot2::labs(x = expression(p[mis][1]), y = expression(p[mis][2]), fill = "Phase-1/2 shared\nSpearman") + ggplot2::theme_classic()
ggplot2::ggsave(file.path(out_dir, "phase2_inference_fidelity_heatmap.png"), p_fidelity, width = 8, height = 7, dpi = 300)

p_overlap <- ggplot2::ggplot(metrics, ggplot2::aes(x = p_mis_phase1, y = p_mis_phase2, fill = jaccard)) +
  ggplot2::stat_summary(fun = mean, geom = "tile") + ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") +
  ggplot2::labs(x = expression(p[mis][1]), y = expression(p[mis][2]), fill = "Support\nJaccard") + ggplot2::theme_classic()
ggplot2::ggsave(file.path(out_dir, "phase2_support_overlap_heatmap.png"), p_overlap, width = 8, height = 7, dpi = 300)

phase1_truth <- metrics[!duplicated(metrics[c("landscape_id", "p1_dir", "replicate_id")]), ]
p_truth <- ggplot2::ggplot(phase1_truth, ggplot2::aes(x = p_mis_phase1, y = phase1_original_grf_spearman)) +
  ggplot2::geom_boxplot(outlier.size = 0.3) + ggplot2::labs(x = expression(p[mis][1]), y = "Phase-1 inference vs original GRF (Spearman)") + ggplot2::theme_classic()
ggplot2::ggsave(file.path(out_dir, "phase1_original_grf_recovery_by_rate.png"), p_truth, width = 10, height = 6, dpi = 300)

peak_long <- rbind(data.frame(metrics[, c("p_mis_phase1", "p_mis_phase2")], metric = "Phase 2 inferred", maxima = metrics$phase2_inferred_local_maxima),
  data.frame(metrics[, c("p_mis_phase1", "p_mis_phase2")], metric = "Original GRF on phase-2 support", maxima = metrics$original_grf_on_phase2_support_local_maxima))
p_peaks <- ggplot2::ggplot(peak_long, ggplot2::aes(x = p_mis_phase2, y = maxima, fill = metric)) +
  ggplot2::geom_boxplot(outlier.size = 0.3) + ggplot2::facet_wrap(~p_mis_phase1, ncol = 4) +
  ggplot2::labs(x = expression(p[mis][2]), y = "Local maxima", fill = "Landscape") + ggplot2::theme_classic()
ggplot2::ggsave(file.path(out_dir, "phase2_local_maxima_by_rate_pair.png"), p_peaks, width = 14, height = 12, dpi = 300)

utils::write.csv(data.frame(figure = c("phase2_terminal_original_grf_fitness_by_rate_pair.png", "phase2_inference_fidelity_heatmap.png", "phase2_support_overlap_heatmap.png", "phase1_original_grf_recovery_by_rate.png", "phase2_local_maxima_by_rate_pair.png")), file.path(out_dir, "figure_manifest.csv"), row.names = FALSE)
