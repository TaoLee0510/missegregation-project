# ------------------------------------------------------------------------------
# ALFA‑K DATA‑PREPARATION PIPELINE
#
# Purpose
#   Convert the longitudinal single‑cell copy‑number dataset of Salehi et.al.
#   (2021) into the chromosome‑level karyotype and lineage files required by
#   the ALFA‑K fitness‑inference framework.
#
# Workflow
#   1. Validate project root and load helper utilities.
#   2. Aggregate per‑bin absolute copy‑number calls to:
#        • whole‑chromosome integer states (used by ALFA‑K)
#        • chromosome‑arm states (generated but not used downstream)
#   3. Reconstruct ancestry paths for every single‑cell library (one uid per
#      sample) and retain all lineage segments ≥ 2 generations deep.
#   4. For every lineage, build a karyotype‑by‑passage count matrix and save it
#      with the assumed passage interval for ALFA‑K.
#
# Key Assumptions
#   • Passage interval = 15days for all models except SA039 & SA906 (i.e. in vitro lines), which are
#     assigned 5days (exact intervals unavailable).
#   • Sex chromosomes and rows with malformed genomic‑bin labels are ignored;
#     sample “SA004” is discarded due to corrupted input files.
#
# Inputs (relative to project root)
#   data/raw/salehi/
#     ├─ arm_loci.Rds            genomic bin lookup table
#     ├─ metadata.csv            sample metadata with parent links & time points
#     └─ raw_post_jump/          <sample>/named_mat.csv.gz per single‑cell library
#
# Outputs (written to data/processed/salehi/)
#   chrom_level_cn.Rds           list of per‑sample 22‑chromosome copy numbers
#   arm_level_cn.Rds             list of per‑sample arm‑level copy numbers
#   lineages.Rds                 named list of reconstructed lineage objects
#   alfak_inputs/
#     └─ <lineage_id>.Rds        list with:
#                                  $x  = karyotype counts by passage
#                                  $dt = assumed passage interval
#
#   The arm‑level file is produced for completeness but is not consumed by
#   downstream ALFA‑K analyses.
# ------------------------------------------------------------------------------


if(!basename(getwd())=="ALFA-K") stop("Ensure working directory set to ALFA-K root")
source("R/utils_env.R")
ensure_packages(c("stringr","R.utils","data.table"))

# Helper to wrap lineage IDs into a list.
# Args:
#   ids: Vector of identifiers for a lineage segment.
# Returns: List containing 'ids'.
lineage_wrapper <- function(ids){
  list(ids=ids)
}

# Generates ancestral lineages for a UID using metadata 'm'.
# Args:
#   uid: Starting unique identifier.
#   m: Metadata dataframe (must contain 'uid' and 'parent' columns).
# Returns: List of lineage segments.
lineage_generator <- function(uid, m){
  ids <- c()
  lineages <- list()
  while(!is.na(uid)){
    ids <- c(uid,ids)
    # Find parent UID from metadata 'm'.
    uid_match_indices <- which(m$uid == uid)
    if (length(uid_match_indices) > 0) {
      uid <- m$parent[uid_match_indices[1]] # Take the first match if multiple
    } else {
      uid <- NA # Stop if uid not found
    }
    if(length(ids)>1) lineages <- c(lineages,lineage_wrapper(ids))
  }
  return(lineages)
}

# Checks for direct (dec1) and second-gen (dec2) descendants.
# Args:
#   ids: Vector of UIDs representing a lineage path.
#   m: Metadata dataframe.
# Returns: List with original 'ids', 'dec1' (direct descendants), 'dec2' (second-gen descendants).
descendent_checker <- function(ids, m){
  last_id <- tail(ids,1)
  dec1 <- m$uid[which(m$parent == last_id & !is.na(m$parent))]
  dec2 <- m$uid[which(m$parent %in% dec1 & !is.na(m$parent))] # Ensure parent is not NA for dec2 as well
  list(ids=ids,dec1=dec1,dec2=dec2)
}

# Creates a standardized name for a lineage.
# Args:
#   lineage: A lineage object (output from descendent_checker).
#   m: Metadata dataframe.
# Returns: String name for the lineage.
lineage_namer <- function(lineage, m){
  suffix <- paste0("_l_",length(lineage$ids),"_d1_",length(lineage$dec1),
                   "_d2_",length(lineage$dec2))
  id <- tail(lineage$ids,1)
  id_match_indices <- which(m$uid == id)
  sid <- "" # Default sid
  if (length(id_match_indices) > 0) {
    idx <- id_match_indices[1] # Take the first match
    sid_datasetname <- gsub(" ", "", m$datasetname[idx], fixed = TRUE)
    sid_timepoint <- m$timepoint[idx]
    sid <- paste(sid_datasetname, sid_timepoint, sep = "_")
    sid <- gsub("/", "", sid, fixed = TRUE)
  }
  paste0(sid,suffix)
}

# Generates, checks descendants, and names all lineages from metadata 'm'.
# Args:
#   m: Metadata dataframe.
# Returns: Named list of lineage objects.
process_lineages <- function(m_arg){
  # Generate all lineages.
  lineages_generated <- do.call(c, lapply(m_arg$uid, function(u) lineage_generator(u, m_arg)))
  # Add descendant information to each lineage.
  # lineage_generator returns list(ids=...). We need to pass l_item$ids.
  lineages_checked <- lapply(lineages_generated, function(l_item) descendent_checker(l_item, m_arg))
  # Name lineages.
  lnames <- sapply(lineages_checked, function(l_chk_item) lineage_namer(l_chk_item, m_arg))
  names(lineages_checked) <- lnames
  return(lineages_checked)
}

# Calculates the statistical mode (most frequent value) of a vector.
# Args:
#   x: A vector.
# Returns: The mode of 'x'.
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Aggregates copy number (CN) matrix at chromosome or arm level using mode.
# Args:
#   x: CN matrix (rows: genomic regions, cols: samples). Row names like "chrom_chrompos_...".
#   arms: Dataframe with chromosome arm info ('chrom', 'arm', 'start').
#   level: "chrom" or "arm" for aggregation level.
# Returns: Aggregated CN matrix (rows: samples, cols: regions).
call_cn <- function(x,arms,level="chrom"){
  df <- do.call(rbind,lapply(rownames(x), function(xi) unlist(strsplit(xi,split="_"))))
  df <- data.frame(df[,1:2])
  colnames(df) <- c("chrom","chrompos")
  df$arm <- sapply(1:nrow(df), function(i){
    chrom_val <- df$chrom[i]
    chrompos_val <- df$chrompos[i]
    # Ensure arms lookup is safe for chrom_val
    q_arm_info <- arms[arms$chrom == chrom_val & arms$arm == "q", "start"]
    is_q <- FALSE # Default if no q_arm_info
    if(length(q_arm_info) > 0) {
      is_q <- q_arm_info[1] < chrompos_val
    }
    is_p <- !is_q
    c("q","p")[c(is_q,is_p)][1] # Return first element, ensure single value
  })
  
  x <- cbind(df[,c("chrom","arm")],x)
  x <- x[!df$chrom%in%c("X","Y"),]
  xcolnames <- colnames(x)
  if(level=="chrom")  x <- split(x,f=x$chrom)
  if(level=="arm") x <- split(x,f=interaction(x$chrom,x$arm))
  lx <- sapply(x,nrow)
  x <- x[lx>0]
  
  x <- data.frame(do.call(rbind,lapply(x,  function(xi){
    if(level=="chrom")  id_val <- xi[1,1]
    if(level=="arm") id_val <- paste(xi[1,1:2],collapse="")
    y <- c(id_val,as.numeric(apply(xi[,3:ncol(xi)],2,Mode)))
    y <- as.character(y)
  })))
  colnames(x) <- c("chrom",xcolnames[-c(1,2)])
  chrom_col_data <- x$chrom # Changed variable name from 'chrom' to 'chrom_col_data'
  x <- t(x[,!colnames(x)=="chrom"])
  colnames(x) <- stringr::str_pad(as.character(chrom_col_data),width=4) # Use 'chrom_col_data' here
  x <- x[,order(colnames(x))]
  return(x)
}

# Extracts and processes CN profiles from raw data.
# Saves chromosome-level and arm-level CN profiles.
# Args:
#   raw_data_root: Root directory for raw data.
#   processed_data_root: Root directory for processed data.
#   arm_loci_filename: Filename for arm loci RDS file (within raw_data_root).
#   raw_cn_dirname: Sub-directory name for raw CN data (within raw_data_root).
# Output: Saves "chrom_level_cn.Rds" and "arm_level_cn.Rds".
extract_cn_profiles <- function(
    raw_data_root = "data/raw/salehi",
    processed_data_root = "data/processed/salehi",
    arm_loci_filename = "arm_loci.Rds",
    raw_cn_dirname = "raw_post_jump"
){
  path2armloci <- file.path(raw_data_root, arm_loci_filename)
  path2raw_base <- file.path(raw_data_root, raw_cn_dirname)
  # Ensure processed_data_root exists for saving output files
  if (!dir.exists(processed_data_root)) {
    dir.create(processed_data_root, recursive = TRUE, showWarnings = FALSE)
  }
  
  arms <- readRDS(path2armloci)
  
  ff_names <- list.files(path2raw_base)
  ff_names <- ff_names[!ff_names=="SA004"] # Exclude problematic file "SA004".
  ff_full_paths <- file.path(path2raw_base, ff_names)
  
  # Chromosome-level CN profiles
  x_chrom <- do.call(rbind,lapply(ff_full_paths, function(fi_path){
    print(fi_path)
    named_mat_path <- file.path(fi_path, "named_mat.csv.gz")
    if(!file.exists(named_mat_path)) return(NULL)
    # Use tryCatch for robust file reading
    current_x <- tryCatch({
      data.frame(data.table::fread(named_mat_path))
    }, error = function(e) {
      warning(paste("Failed to read or process:", named_mat_path, "-", e$message))
      return(NULL)
    })
    if(is.null(current_x) || nrow(current_x) == 0) return(NULL)
    
    binformat <- unlist(strsplit(as.character(current_x$bin_name[1]),split="_"))
    if(length(binformat)<3) return(NULL)
    rownames(current_x) <- current_x$bin_name
    current_x <- current_x[,!colnames(current_x)=="bin_name"]
    call_cn(current_x,arms)
  }))
  
  if (!is.null(x_chrom) && nrow(x_chrom) > 0) {
    cx_chrom <- colnames(x_chrom)
    x_chrom <- t(apply(x_chrom,1,as.numeric))
    colnames(x_chrom) <- cx_chrom
    ids_chrom <- rownames(x_chrom)
    ids_chrom <- sapply(ids_chrom, function(idi) unlist(strsplit(idi,split="[.]"))[2])
    x_chrom_df <- data.frame(x_chrom,check.names=F)
    x_chrom_list <- split(x_chrom_df,f=ids_chrom)
    saveRDS(x_chrom_list,file=file.path(processed_data_root,"chrom_level_cn.Rds"))
  } else {
    warning("No chromosome-level CN data processed or available to save.")
  }
  
  # Arm-level CN profiles
  x_arm <- do.call(rbind,lapply(ff_full_paths, function(fi_path){
    print(fi_path)
    named_mat_path <- file.path(fi_path, "named_mat.csv.gz")
    if(!file.exists(named_mat_path)) return(NULL)
    current_x <- tryCatch({
      data.frame(data.table::fread(named_mat_path))
    }, error = function(e) {
      warning(paste("Failed to read or process:", named_mat_path, "-", e$message))
      return(NULL)
    })
    if(is.null(current_x) || nrow(current_x) == 0) return(NULL)
    
    binformat <- unlist(strsplit(as.character(current_x$bin_name[1]),split="_"))
    if(length(binformat)<3) return(NULL)
    rownames(current_x) <- current_x$bin_name
    current_x <- current_x[,!colnames(current_x)=="bin_name"]
    call_cn(current_x,arms,level = "arm")
  }))
  
  if (!is.null(x_arm) && nrow(x_arm) > 0) {
    cx_arm <- colnames(x_arm)
    x_arm <- t(apply(x_arm,1,as.numeric))
    colnames(x_arm) <- cx_arm
    ids_arm <- rownames(x_arm)
    ids_arm <- sapply(ids_arm, function(idi) unlist(strsplit(idi,split="[.]"))[2])
    x_arm_df <- data.frame(x_arm,check.names=F)
    x_arm_list <- split(x_arm_df,f=ids_arm)
    saveRDS(x_arm_list,file=file.path(processed_data_root,"arm_level_cn.Rds"))
  } else {
    warning("No arm-level CN data processed or available to save.")
  }
}

# Extracts lineages, associates with CN profiles, formats for ALFA-K.
# Args:
#   raw_data_root: Root directory for raw data.
#   processed_data_root: Root directory for processed data.
#   metadata_filename: Filename for metadata CSV (within raw_data_root).
#   lineages_rds_name: Filename for saving processed lineages (within processed_data_root).
#   alfak_subdir_name: Subdirectory for ALFA-K inputs (within processed_data_root).
#   chrom_cn_rds_name: Filename for chromosome-level CN RDS (within processed_data_root).
# Output: Saves lineages RDS and ALFA-K input files.
extract_lineages <- function(
    raw_data_root = "data/raw/salehi",
    processed_data_root = "data/processed/salehi",
    metadata_filename = "metadata.csv",
    lineages_rds_name = "lineages.Rds",
    alfak_subdir_name = "alfak_inputs",
    chrom_cn_rds_name = "chrom_level_cn.Rds"
){
  m_path <- file.path(raw_data_root, metadata_filename)
  # Ensure processed_data_root exists for saving output files
  if (!dir.exists(processed_data_root)) {
    dir.create(processed_data_root, recursive = TRUE, showWarnings = FALSE)
  }
  
  m_local <- read.csv(m_path) # Load metadata locally.
  lineages_data <- process_lineages(m_local) # Pass local 'm_local'.
  saveRDS(lineages_data, file.path(processed_data_root, lineages_rds_name))
  
  outdir_path <- file.path(processed_data_root, alfak_subdir_name)
  dir.create(outdir_path, showWarnings = FALSE, recursive = TRUE)
  
  cnmat_path <- file.path(processed_data_root, chrom_cn_rds_name)
  if (!file.exists(cnmat_path)) {
    warning(paste("CN matrix file not found, skipping ALFA-K input generation:", cnmat_path))
    return()
  }
  cnmat <- readRDS(cnmat_path)
  
  for(i in 1:length(lineages_data)){
    lineage_id_name <- names(lineages_data)[i]
    print(i)
    current_lineage <- lineages_data[[i]]
    
    # Use local 'm_local' for metadata operations within this function.
    msel <- m_local[m_local$uid %in% current_lineage$ids,]
    if (nrow(msel) == 0) next # Skip if no metadata for this lineage
    
    pdx_id <- msel$PDX_id[1]
    dt <- 15
    if(pdx_id %in% c("SA039","SA906")) dt <- 5
    
    msel_expanded <- do.call(rbind,lapply(1:nrow(msel), function(row_idx){
      lid <- unlist(strsplit(as.character(msel$library_ids[row_idx]),split=";"))
      data.frame(timepoint=msel$timepoint[row_idx],library_id=lid,row.names = NULL, stringsAsFactors = FALSE)
    }))
    
    if(nrow(msel_expanded) == 0 || mean(msel_expanded$library_id %in% names(cnmat))<1) next
    
    msel_filtered <- msel_expanded[msel_expanded$library_id %in% names(cnmat),]
    if (nrow(msel_filtered) == 0) next
    
    x_cn_subset <- cnmat[msel_filtered$library_id]
    
    # Check if x_cn_subset contains non-NULL and non-empty data frames before proceeding
    valid_cn_data_indices <- sapply(x_cn_subset, function(df) !is.null(df) && nrow(df) > 0)
    if (!any(valid_cn_data_indices)) next
    
    x_cn_subset <- x_cn_subset[valid_cn_data_indices]
    msel_filtered <- msel_filtered[valid_cn_data_indices, ] # Keep msel_filtered consistent
    
    if (nrow(msel_filtered) == 0) next # Check again after filtering based on valid_cn_data
    
    # Ensure number of rows in each df of x_cn_subset matches before rep
    # This assumes all dataframes in x_cn_subset for a given timepoint group should have same number of rows (e.g. bins/chroms)
    # For simplicity, we'll take the number of rows from the first valid df, assuming consistency
    # A more robust check might be needed if row counts vary unexpectedly.
    
    ids_map <- unlist(lapply(1:nrow(msel_filtered), function(j) {
      if (!is.null(x_cn_subset[[j]]) && nrow(x_cn_subset[[j]]) > 0) {
        rep(msel_filtered$timepoint[j], nrow(x_cn_subset[[j]]))
      } else {
        NULL
      }
    }))
    
    x_combined_cn <- do.call(rbind,x_cn_subset)
    if (is.null(x_combined_cn) || nrow(x_combined_cn) == 0) next
    
    k_profiles <- apply(x_combined_cn,1,paste,collapse=".")
    x0_unique_profiles <- unique(k_profiles)
    
    x_timepoints_split <- split(k_profiles,f=ids_map)
    
    x_final_matrix <- do.call(cbind,lapply(x_timepoints_split, function(xi_timepoint){
      sapply(x0_unique_profiles, function(xj_profile) sum(xi_timepoint==xj_profile))
    }))
    if (is.null(x_final_matrix) || ncol(x_final_matrix) == 0) next
    
    x_final_matrix <- x_final_matrix[order(rowSums(x_final_matrix),decreasing=T),,drop=FALSE] # keep as matrix
    colnames(x_final_matrix) <- gsub("X","",colnames(x_final_matrix))
    # Ensure column names are numeric before ordering
    numeric_colnames <- suppressWarnings(as.numeric(colnames(x_final_matrix)))
    if (!any(is.na(numeric_colnames))) {
      x_final_matrix <- x_final_matrix[,order(numeric_colnames),drop=FALSE]
    }
    
    x_output <- list(x=x_final_matrix,pop.fitness=NULL,dt=dt)
    saveRDS(x_output,file.path(outdir_path,paste0(lineage_id_name,".Rds")))
  }
}

# --- Main Script Execution ---
# Set appropriate working directory.
# setwd("~/projects/ALFA-K/") # Example: Uncomment and set your path

# Load metadata once for operations outside function-specific loads.
# This 'm_main' is passed to functions that require metadata.
# It's assumed "data/raw/salehi/metadata.csv" is relative to the working directory.
# Make sure the file paths for raw data used by extract_cn_profiles and extract_lineages
# are correct relative to your working directory, or provide absolute paths.

# Default file structure:
# getwd()
#   |- data/
#      |- salehi/
#         |- raw/
#         |  |- arm_loci.Rds
#         |  |- metadata.csv
#         |  |- raw_post_jump/ (contains subdirs with named_mat.csv)
#         |- processed/
#            |- (outputs will be saved here)

main_metadata_path <- "data/raw/salehi/metadata.csv" # Adjusted path assuming metadata is with other raw files.

if (file.exists(main_metadata_path)) {
  m_main <- read.csv(main_metadata_path)
  
  # Process lineages using the main metadata.
  lineages_processed <- process_lineages(m_main)
  
  extract_cn_profiles(
    raw_data_root = "data/raw/salehi",
    processed_data_root = "data/processed/salehi",
    arm_loci_filename = "arm_loci.Rds",
    raw_cn_dirname = "raw_post_jump"
  )
  
  # Extract lineage-specific data for ALFA-K using streamlined path arguments.
  extract_lineages(
    raw_data_root = "data/raw/salehi",
    processed_data_root = "data/processed/salehi",
    metadata_filename = "metadata.csv", # This metadata is loaded inside extract_lineages
    lineages_rds_name = "lineages.Rds",
    alfak_subdir_name = "alfak_inputs",
    chrom_cn_rds_name = "chrom_level_cn.Rds"
  )
} else {
  warning(paste("Main metadata file not found:", main_metadata_path, "- Script execution might be incomplete."))
}