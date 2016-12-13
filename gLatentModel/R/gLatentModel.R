gLatentModel <- function(dat, K){
  n <- nrow(dat); d <- ncol(dat)

  gamma_mat <- estimate_gamma(dat)
  cov_mat <- cov(dat)
  c_mat <- cov_mat - gamma_mat

  pure_idx <- pure_nodes(c_mat, K)
  idx_rest <- c(1:d)[-pure_idx]; new_idx <- c(pure_idx, idx_rest)
  dat <- dat[new_idx, ]

  theta_mat <- estimate_theta(c_mat, pure_idx)
  a_mat <- estimate_a(theta_mat, c_mat, pure_idx)
  group_list <- group_cluster(a_mat)

  gamma_mat <- reestimate_gamma(dat, group_list)
}
