gLatentModel <- function(dat, K, verbose = F){
  n <- nrow(dat); d <- ncol(dat)

  gamma_vec <- estimate_gamma(dat)
  cov_mat <- cov(dat)
  c_mat <- cov_mat - diag(gamma_vec)
  if(verbose) print("Finished step 1: estimating gamma")

  pure_idx <- pure_nodes(c_mat, K)
  idx_rest <- c(1:d)[-pure_idx]; new_idx <- c(pure_idx, idx_rest)
  idx_original <- order(new_idx)
  if(verbose) print("Finished step 2: finding pure nodes")

  theta_mat <- estimate_theta(c_mat, pure_idx)
  a_mat <- estimate_a(theta_mat, c_mat, pure_idx)
  a_mat <- a_mat[idx_original,]
  partition_list <- partition_cluster(a_mat)
  if(verbose) print("Finished step 3: estimating A")

  gamma_vec <- reestimate_gamma(dat, partition_list)
  if(verbose) print("Finished step 4: reestimating gamma")

  theta_mat <- reestimate_theta(cov_mat, diag(gamma_vec), partition_list)

  structure(list(theta = theta_mat, a = a_mat,
    partition_list = partition_list, gamma = gamma_vec),
    class = "gLatentModel")
}
