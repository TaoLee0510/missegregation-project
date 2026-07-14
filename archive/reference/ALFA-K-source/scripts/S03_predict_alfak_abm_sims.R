if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_env.R")
ensure_packages(c("alfakR","stringr"))
Nreps <- 5 ## number of ABM replicates
rawDir <- "data/raw/ABM"
procDir <- "data/processed/ABM"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: Rscript name_of_this_script.R <input_file>")

fi_full <- args[1]
if (!file.exists(fi_full)) stop("Input file not found: ", fi_full)

fi0 <- basename(fi_full)# keep only filename

fi_stripped <- gsub(".Rds","",fi0)

fi <- file.path(rawDir,fi0)
base_outdir <- file.path(procDir,fi_stripped)
fit_paths <- list.files(base_outdir)

##determine initial conditions and eval times (same for all fits per ABM run)
xi <- readRDS(fi)
all_times <- as.numeric(colnames(xi$abm_output$x))
start_time <- max(all_times[all_times<120])
report_times <- head(all_times[all_times>start_time],10)-start_time
x0 <- xi$abm_output$x[,colnames(xi$abm_output$x)==start_time]
names(x0) <- rownames(xi$abm_output$x)

x0_base <- x0



## run predictions for each fit
for(path in fit_paths){
  for(rep in 1:Nreps){
    alfak_fit_path <- file.path(base_outdir,path,"landscape.Rds")
    if(!file.exists(alfak_fit_path)) next;
    abm_out_dir <- file.path(base_outdir,path,"abm_preds")
    dir.create(abm_out_dir)  
    
    abm_save_path <- file.path(abm_out_dir,
                               paste0("rep_",stringr::str_pad(rep,2,pad="0"),".Rds"))
    
    xj <- readRDS(alfak_fit_path)
    
    x0 <- x0_base
    x0 <- x0[names(x0)%in%xj$k]
    x0 <-x0/sum(x0)
    xtmp <- 0*xj$mean
    names(xtmp) <- xj$k
    xtmp[names(x0)]<-x0
    x0 <- xtmp
    
    ## run prediction (pop size and delta_t parms etc. to match initial sims)
    abm_delta_t <- 0.1
    y <- predict_evo(lscape=xj,p=5e-05,times=c(0,report_times),x0=x0,prediction_type = "ABM",
                     abm_pop_size = 2e5,abm_max_pop = 2e7,abm_delta_t = abm_delta_t,
                     abm_culling_survival = 0.01,
                     abm_record_interval = 10,
                     abm_seed = rep)
    
    ## match to evaluation timepoints
    ix <- sapply(report_times,function(i) which.min(abs(i-y$time)))
    y <- y[c(1,ix),]
    
    saveRDS(y,abm_save_path)
  }
}













