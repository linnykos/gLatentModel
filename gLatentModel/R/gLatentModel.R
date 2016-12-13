gLatentModel <- function(dat, K){
  n <- nrow(dat); d <- ncol(dat)

  gamma_vec <- estimate_gamma(dat)
  cov_mat <- cov(dat)
  c_mat <- cov_mat - diag(gamma_vec)

  pure_idx <- pure_nodes(c_mat, K)
  idx_rest <- c(1:d)[-pure_idx]; new_idx <- c(pure_idx, idx_rest)
  dat <- dat[new_idx, ]; cov_mat <- cov_mat[new_idx, new_idx]

  theta_mat <- estimate_theta(c_mat, pure_idx)
  a_mat <- estimate_a(theta_mat, c_mat, pure_idx)
  partition_list <- partition_cluster(a_mat)

  gamma_vec <- reestimate_gamma(dat, partition_list)

  theta <- reestimate_theta(cov_mat, gamma_mat, a_mat)
}
