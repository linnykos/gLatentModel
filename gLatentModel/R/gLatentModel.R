#' Estimate G Latent Model
#'
#' @param dat n by d matrix where there are n samples and d variables
#' @param K number of clusters for the d variables
#' @param verbose boolean to print out progress statements
#' @param seed seed to fix the random subsampling and kmeans within estimation
#' @param num_subsample number of subsamples when initially estimate variance of noise
#'
#' @return a gLatentModel object
#' @export
gLatentModel <- function(dat, K, verbose = F, seed = NA, num_subsample = NA,
                         debugging = F){
  n <- nrow(dat); d <- ncol(dat)

  if(!is.na(seed)) set.seed(seed)
  cov_noise <- .estimate_cov_noise(dat, num_subsample)
  cov_dat <- stats::cov(dat)
  cov_denoise <- cov_dat - diag(cov_noise)
  if(verbose) print("Finished step 1: estimating noise variance")

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

  assignment_mat_2 <- .reestimate_assignment(cluster_vec)
  cov_noise_2 <- .reestimate_cov_denoise(dat, cluster_vec, cov_noise)
  cov_latent_2 <- .reestimate_cov_latent(cov_dat, diag(cov_noise_2), assignment_mat_2)
  if(verbose) print("Finished step 4: reestimating")

  if(!debugging) {
    structure(list(cov_latent = cov_latent_2, cluster = cluster_vec), class = "gLatentModel")
  } else {
    structure(list(cov_latent = cov_latent_2, cluster = cluster_vec,
                   assignment_mat_org = assignment_mat, cov_latent_org = cov_latent,
                   cov_noise = cov_noise_2, cov_noise_org = cov_noise,
                   pure_idx = pure_idx), class = "gLatentModel")
  }
}
