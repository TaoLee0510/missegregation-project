if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_env.R")
ensure_packages(c("transport"))
source("R/utils_karyo.R")

R2R <- function(obs, pred) {
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
}

summarize_results <- function(dj,fi,
                              rawDir="data/raw/ABM",
                              fitDir = "data/processed/ABM")
{
  ## setup paths to fitted landscape, predictions, raw data etc.
  fi <- basename(fi)
  abm_sim_id <- gsub(".Rds","",fi)
  fi <- file.path(rawDir,fi)
  fj <- file.path(fitDir,abm_sim_id,dj,"landscape.Rds")
  fx <- file.path(fitDir,abm_sim_id,dj,"xval.Rds")
  preds_path <- file.path(fitDir,abm_sim_id,dj,"abm_preds/")
  preds_paths <- paste0(preds_path,list.files(preds_path))
  save_path <- file.path(fitDir,abm_sim_id,dj,"abm_preds_summary.Rds")
  
  ## extract metadata from filenames
  abm_ids <- strsplit(abm_sim_id,split="_") |> unlist()
  wavelength <- gsub("p",".",abm_ids[2]) |> as.numeric()
  abmrep <- abm_ids[6]
  parm_ids <- strsplit(dj,split="_") |> unlist()
  
  minobs <- parm_ids[2]
  ntp <- parm_ids[4]
  
  info <- data.frame(w=wavelength,abmrep,minobs,ntp)
  
  xi <- readRDS(fi)
  lscape <- readRDS(fj)
  tru_lscape <- xi$true_landscape
  Rxv <- readRDS(fx)
  
  f <- sapply(lscape$k,function(k){
    kvec <- strsplit(k,split="[.]") |> unlist() |> as.numeric()
    getf(kvec,wavelength,tru_lscape)
  })
  
  summary <- data.frame(r=cor(lscape$mean,f),
                        rfq=cor(lscape$mean[lscape$fq],f[lscape$fq]),
                        rnn=cor(lscape$mean[lscape$nn],f[lscape$nn]),
                        rd2=cor(lscape$mean[!(lscape$fq|lscape$nn)],f[!(lscape$fq|lscape$nn)]),
                        rho=cor(lscape$mean,f,method = "spearman"),
                        rhofq=cor(lscape$mean[lscape$fq],f[lscape$fq],method = "spearman"),
                        rhonn=cor(lscape$mean[lscape$nn],f[lscape$nn],method = "spearman"),
                        rhod2=cor(lscape$mean[!(lscape$fq|lscape$nn)],f[!(lscape$fq|lscape$nn)],method = "spearman"),
                        R=R2R(f,lscape$mean),
                        Rfq=R2R(f[lscape$fq],lscape$mean[lscape$fq]),
                        Rnn=R2R(f[lscape$nn],lscape$mean[lscape$nn]),
                        Rd2=R2R(f[!(lscape$fq|lscape$nn)],lscape$mean[!(lscape$fq|lscape$nn)]),
                        Rxv=Rxv,
                        nfq=sum(lscape$fq))
  
  preds <- lapply(preds_paths,readRDS)
  
  preds <- as.data.frame(Reduce(`+`, lapply(preds, function(df) as.matrix(df[-1]))) / length(preds))
  
  yi <- xi$abm_output$x
  
  tt <- as.numeric(colnames(yi))
  t0 <- max(tt[tt<120])
  measure_tps <- head(tt[tt>=t0],11)
  
  yi <- yi[,colnames(yi)%in%measure_tps]
  
  y <- preds*0
  rx <- rownames(yi)[rownames(yi)%in%lscape$k]
  
  for(i in 1:length(measure_tps)){
    y[i,rx] <- yi[rx,i]/sum(yi[rx,i])
  }
  
  y0 <- y[1,]
  
  k <- do.call(rbind,lapply(lscape$k,function(k){strsplit(k,split="[.]") |> unlist() |> as.numeric()}))
  
  k0 <- colSums(as.numeric(y0)*k)
  
  mag <- function(v) sqrt(sum(v^2))
  getangle <- function(a,b) 180*acos(sum(a*b)/(mag(a)*mag(b)))/pi
  df_a <- do.call(rbind,lapply(1:nrow(y),function(i){
    ky <- colSums(as.numeric(y[i,])*k)
    kp <- colSums(as.numeric(preds[i,])*k)
    
    data.frame(angle=getangle(ky-k0,kp-k0))
    
  }))
  rownames(df_a) <- measure_tps
  
  df_de <- do.call(rbind,lapply(1:nrow(y),function(i){
    ky <- colSums(as.numeric(y[i,])*k)
    kp <- colSums(as.numeric(preds[i,])*k)
    
    data.frame(de_p=sqrt(sum((kp-ky)^2)),
               de_0=sqrt(sum((k0-ky)^2)))
    
  }))
  rownames(df_de) <- colnames(yi$data$x)[measure_tps]
  
  s0 <- wpp(k,mass=y0)
  df_dw <- do.call(rbind,lapply(1:nrow(y),function(i){
    
    tryCatch({
      sy <- wpp(k,mass=y[i,])
      sp <- wpp(k,mass=preds[i,])
      data.frame(dw_p=wasserstein(sy,sp),
                 dw_0=wasserstein(sy,s0))
    },error=function(e){
      data.frame(dw_p=NaN,
                 dw_0=NaN)
    })
    
    
  }))
  rownames(df_dw) <- measure_tps
  
  df_sc <- do.call(rbind,lapply(1:nrow(y),function(i){
    sc_0 <- sum(y0*y[i,])/prod(sqrt(sum(y0^2)) , sqrt(sum(y[i,]^2)))
    sc_p <- sum(preds[i,]*y[i,])/prod(sqrt(sum(preds[i,]^2)) , sqrt(sum(y[i,]^2)))
    data.frame(sc_p,sc_0)
  }))
  rownames(df_sc) <- measure_tps
  df_so <- do.call(rbind,lapply(1:nrow(y),function(i){
    so_0 <- sum(pmin(y0, y[i,]))/min(sum(y0), sum(y[i,]))
    so_p <- sum(pmin(preds[i,], y[i,]))/min(sum(preds[i,]), sum(y[i,]))
    data.frame(so_p,so_0)
  }))
  rownames(df_so) <- measure_tps
  
  
  result <- list(info=info,summary=summary,df_so=df_so,df_sc=df_sc,df_de=df_de,df_dw=df_dw,df_a=df_a)
  saveRDS(result,save_path)
  
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: Rscript name_of_this_script.R <input_file>")

fi_full <- args[1]
if (!file.exists(fi_full)) stop("Input file not found: ", fi_full)


fi <- basename(fi_full)# keep only filename
fitDir = "data/processed/ABM"
fi_strp <- gsub(".Rds","",fi)
targetDir <- file.path(fitDir,fi_strp)

pardirs <- list.files(targetDir)

lapply(pardirs,function(dj){
  tryCatch(summarize_results(dj,fi),error=function(e) print(paste("error for",fi,dj,":\n",e)))
})
