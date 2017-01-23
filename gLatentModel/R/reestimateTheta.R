.reestimate_theta <- function(cov_mat, gamma_mat, partition_list){
  stopifnot(ncol(gamma_mat) == nrow(gamma_mat))
  stopifnot(nrow(cov_mat) == nrow(gamma_mat))

  a_ind <- matrix(0, nrow(cov_mat), length(partition_list))
  for(i in 1:length(partition_list)){
    a_ind[partition_list[[i]],i] <- 1
  }

  a_inv <- solve(t(a_ind) %*% a_ind)
  a_inv %*% t(a_ind) %*% (cov_mat - gamma_mat) %*% a_ind %*% a_inv
}
