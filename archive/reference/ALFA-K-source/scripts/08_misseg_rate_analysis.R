if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_env.R")
ensure_packages("alfakR")

procDir <- "data/processed/salehi/alfak_outputs_proc/"
files <- list.files(procDir, full.names = TRUE)
x_list <- lapply(files, readRDS)

x <- do.call(rbind, x_list)

# Filter data and select observation with maximum xv per file and per clean_id
x <- x[!is.na(x$xv) & x$xv > 0 & x$test_treat != "NaN", ]
x <- do.call(rbind, lapply(split(x, x$fi), function(df) df[df$xv == max(df$xv), ]))
rownames(x) <- NULL
x$clean_id <- sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", x$fi)
x <- do.call(rbind, lapply(split(x, x$clean_id), function(df) head(df[df$xv == max(df$xv), ],1)))
x <- x[!duplicated(x$fi),]


df <- lapply(1:nrow(x),function(i){
  inpath <- paste0("data/processed/salehi/alfak_outputs/",x$fi[i],"/minobs_",
                   x$min_obs[i],"/landscape.Rds")
  
  
  if(!file.exists(inpath)) return(NULL)
  lscape <- readRDS(inpath)
  
  
  p <- 0.00002*2^(1:8)
  res <- do.call(cbind,pbapply::pblapply(p,function(p0){
    ## Nmax is the maximum allowable number of missegregations.
    alfakR::find_steady_state(lscape = lscape,p = p0,Nmax=1)

  }))
  colnames(res) <- p
  res <- res[order(rowSums(res),decreasing=T),]
  print(c(i,x$fi[i]))
  print(head(round(res,digits = 2)))
  
  dfi <- data.frame(fi=x$fi[i],nmax=length(unique(apply(res,2,which.max))))
  print(dfi)
  return(res)
})

names(df) <- x$fi

for(i in 1:length(df)){
  attr(df[[i]],"min_obs") <- x$min_obs[i]
  attr(df[[i]],"fi") <- x$fi[i]
}

saveRDS(df,"data/processed/salehi/steadyStatePredictions.Rds")

