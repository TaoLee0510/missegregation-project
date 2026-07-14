###############################################################################
# # Utility functions for “angle‑vs‑random‑orientation” hypothesis tests.
#
# Reference problem
# -----------------
# We have several *clusters* (experimental lineages, sister‑pair IDs, or
# unrelated‑pair IDs).  Each cluster contributes one or more angles θ∈[0°,180°].
# Question: is the **median‑of‑cluster‑medians** significantly smaller than the
# 90° median expected if directions were random in an N‑dimensional space?
#
# Functions
# ---------
# r_sphere_angle(n, N = 22)
#   • Draw `n` random angles (degrees) from the analytic distribution of angles
#     between two random unit vectors on the (N–1)‑sphere.
#
# sphere_null_test(angles, cluster, N = 22, B = 1e4, side = "left")
#   • angles  : numeric vector of angles (0–180).
#   • cluster : factor/character vector; same length as `angles`; defines which
#               rows belong to the same independent cluster.
#   • N       : embedding dimension (22 for karyotype vectors).
#   • B       : Monte‑Carlo sample size (≥ 10 000 recommended).
#   • side    : "left" → H₁: median < 90° (default); "two" for two‑sided test.
#
#   Returns a named list with
#     $stat   – observed median‑of‑medians
#     $p      – Monte‑Carlo p‑value
#     $null   – vector of length B (null distribution, invisible unless needed)
#
# Example usage
# -------------
# pred_res  <- sphere_null_test(z_euc$angle, z_euc$type)           # predictions
# sis_res   <- sphere_null_test(y$value[y$metric=="angle"],
#                               y$lineage[y$metric=="angle"])      # sisters
# unrel_res <- sphere_null_test(z_pairs$angle,
#                               paste(z_pairs$PDX1, z_pairs$PDX2)) # unrelated
###############################################################################

r_sphere_angle <- function(n, N = 22){
  r <- (N - 1) / 2                     # Beta shape parameter
  acos(2 * rbeta(n, r, r) - 1) * 180 / pi
}

sphere_null_test <- function(angles, cluster, N = 22, B = 1e4, side = "left",return_null=F){
  stopifnot(length(angles) == length(cluster))
  cl_sizes <- table(cluster)                          # n_g per cluster
  ## 1. observed statistic
  obs_stat <- median(tapply(angles, cluster, median))
  ## 2. Monte‑Carlo null
  null_stat <- replicate(B, {
    fake_meds <- vapply(cl_sizes,
                        function(n) median(r_sphere_angle(n, N)),
                        numeric(1))
    median(fake_meds)
  })
  ## 3. p‑value
  if(side == "left"){
    p_val <- mean(null_stat <= obs_stat)              # H₁: smaller
  } else if(side == "two"){
    p_val <- mean(abs(null_stat - 90) >= abs(obs_stat - 90))
  } else stop("side must be 'left' or 'two'")
  if(!return_null) null_stat <- NULL
  list(stat = obs_stat, p = p_val, null = null_stat)
}
###############################################################################
