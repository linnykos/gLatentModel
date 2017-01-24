gLatentModel <- function(dat, K, verbose = F, seed = NA, num_subsample = NA){
  n <- nrow(dat); d <- ncol(dat)

  if(!is.na(seed)) set.seed(seed)
  cov_noise <- .estimate_cov_noise(dat, num_subsample)
  cov_dat <- stats::cov(dat)
  cov_denoise <- cov_dat - diag(cov_noise)
  if(verbose) print("Finished step 1: estimating gamma")

  if(!is.na(seed)) set.seed(seed)
  pure_idx <- pure_nodes(cov_denoise, K)
  idx_rest <- c(1:d)[-pure_idx]; new_idx <- c(pure_idx, idx_rest); idx_original <- order(new_idx)
  if(verbose) print("Finished step 2: finding pure nodes")

  cov_latent <- .estimate_cov_latent(cov_denoise, pure_idx)
  assignment_mat <- .estimate_assignment(cov_latent, cov_denoise, pure_idx)
  assignment_mat <- assignment_mat[idx_original,]

  if(!is.na(seed)) set.seed(seed)
  partition_mat <- .estimate_partition_mat(assignment_mat, cov_latent)
  cluster_vec <- .estimate_cluster(partition_mat, K)
  if(verbose) print("Finished step 3: estimating clusters")

  assignment_mat <- .reestimate_assignment(cluster_vec)
  cov_denoise <- .reestimate_cov_denoise(dat, cluster_vec, cov_denoise)
  cov_latent <- .reestimate_cov_latent(cov_dat, diag(cov_denoise), assignment_mat)
  if(verbose) print("Finished step 4: reestimating")

  structure(list(cov_latent = cov_latent, cluster = cluster_vec), class = "gLatentModel")
}
