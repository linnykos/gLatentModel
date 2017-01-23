gLatentModel <- function(dat, K, verbose = F, seed = NA, num_subsample = NA){
  n <- nrow(dat); d <- ncol(dat)

  if(!is.na(seed)) set.seed(seed)
  cov_noise <- .estimate_cov_noise(dat, num_subsample)
  cov_dat <- stats::cov(dat)
  cov_denoise <- cov_dat - diag(cov_noise)
  if(verbose) print("Finished step 1: estimating gamma")

  if(!is.na(seed)) set.seed(seed)
  pure_idx <- pure_nodes(cov_denoise, K)
  idx_rest <- c(1:d)[-pure_idx]; new_idx <- c(pure_idx, idx_rest)
  idx_original <- order(new_idx)
  if(verbose) print("Finished step 2: finding pure nodes")

  cov_latent <- .estimate_cov_latent(cov_denoise, pure_idx)
  assignment_mat <- .estimate_assignment(cov_latent, cov_denoise, pure_idx)
  assignment_mat <- assignment_mat[idx_original,]
  partition_list <- .partition_cluster(assignment_mat)
  if(verbose) print("Finished step 3: estimating A")

  cov_denoise <- .reestimate_gamma(dat, partition_list, cov_denoise)
  if(verbose) print("Finished step 4: reestimating gamma")

  cov_latent <- .reestimate_theta(cov_dat, diag(cov_denoise), partition_list)

  cluster <- rep(0, d)
  for(i in 1:K){
    cluster[partition_list[[i]]] <- i
  }

  structure(list(cov_latent = cov_latent, assignment_mat = assignment_mat,
    partition_list = partition_list, cluster = cluster, cov_denoise = cov_denoise),
    class = "gLatentModel")
}
