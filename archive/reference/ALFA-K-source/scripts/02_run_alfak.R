#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# ALFA‑K batch runner
#
# This script wraps a *single* ALFA‑K execution so that it can be launched
# as an sbatch array element.  It:
#   1.Takes two command‑line arguments –`minobs` and the filename (`fi`)
#   2.Runs `alfak()` and both ODE and ABM forward predictions
#   3.Writes all outputs underdata/processed/salehi/alfak_outputs/<id>/minobs_<minobs>
#
# ---- SLURM template ------------------------------------------------------
## see companion sh script
# -------------------------------------------------------------------------
if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_env.R")
ensure_packages("alfakR")

# ---- Parse arguments -----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript run_alfak.R <minobs> <input_file>")

minobs <- as.numeric(args[1])
if (is.na(minobs) || minobs <= 0) stop("minobs must be a positive integer")

fi_full <- args[2]
if (!file.exists(fi_full)) stop("Input file not found: ", fi_full)

fi <- basename(fi_full)                        # keep only filename
id <- sub("\\.Rds$", "", fi)                  # sample identifier


in_dir   <- "data/processed/salehi/alfak_inputs/"
meta     <- read.csv("data/raw/salehi/metadata.csv")
lineages <- readRDS("data/processed/salehi/lineages.Rds")

out_dir  <- file.path("data/processed/salehi/alfak_outputs",
                      id, paste0("minobs_", minobs))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Running ALFA‑K for ", fi, "(minobs=", minobs, ")")

# ---- ALFA‑K parameters ---------------------------------------------------
pred_iters <- 100; nboot <- 100
n0 <- 1e5; nb <- 1e7; pm <- 5e-5          # NB pm = 0.00005
yi <- readRDS(file.path(in_dir, fi))

ix <- strsplit(fi, "_")[[1]]
dataset <- ix[1]; passage <- ix[2]
if (grepl("^SA535", fi)) { dataset <- paste(ix[1:3], collapse = "_"); passage <- ix[4] }

pass_times <- min(as.numeric(colnames(yi$x))):max(as.numeric(colnames(yi$x)))

# ---- Run ALFA‑K ----------------------------------------------------------
message("alfak()")
tryCatch(
  alfak(yi, out_dir, pass_times, minobs = minobs,
        nboot = nboot, n0 = n0, nb = nb, pm = pm),
  error = function(e) { message("alfak() failed: ", e$message); quit(save = "no", status = 0) }
)

# ---- Predict next passage if lineage info available ----------------------
lin <- lineages[[id]]
if (length(lin$dec1) == 0) quit(save = "no", status = 0)  # nothing to predict

parent_p  <- as.numeric(sub("X", "", meta$timepoint[meta$uid %in% tail(lin$ids, 1)]))
daughter_p<- as.numeric(sub("X", "", meta$timepoint[meta$uid %in% lin$dec1]))
days_per_p<- ifelse(meta$PDX_id[meta$uid %in% tail(lin$ids, 1)] %in% c("SA906", "SA039"), 5, 15)

sim_time  <- (max(daughter_p) - parent_p) * days_per_p
sim_times <- seq(0, sim_time, by = 1)

ls_file <- file.path(out_dir, "landscape.Rds")
if (!file.exists(ls_file)) quit(save = "no", status = 0)
lscape <- readRDS(ls_file)

x0 <- rep(0, nrow(lscape)); names(x0) <- lscape$k
input_last <- yi$x[rownames(yi$x) %in% lscape$k, ncol(yi$x)]
x0[names(input_last)] <- input_last; x0 <- x0 / sum(x0)

# ---- ODE prediction ------------------------------------------------------
message("ODE prediction")
ode_res <- predict_evo(lscape, p = pm, times = sim_times,
                       x0 = x0, prediction_type = "ODE")
saveRDS(ode_res, file.path(out_dir, "ode_preds.Rds"))

# ---- ABM predictions -----------------------------------------------------
abm_dir <- file.path(out_dir, "abm_preds"); dir.create(abm_dir, showWarnings = FALSE)
for (rep in 1:5) {
  abm_res <- predict_evo(lscape, p = pm, times = sim_times, x0 = x0,
                         prediction_type = "ABM",
                         abm_pop_size = 1e5, abm_delta_t = 0.1,
                         abm_max_pop = nb, abm_record_interval = 10,abm_culling_survival = 0.01)
  saveRDS(abm_res, file.path(abm_dir, paste0("rep_", rep, ".Rds")))
}

message("Done")


