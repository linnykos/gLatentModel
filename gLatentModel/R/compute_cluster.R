.estimate_partition_mat <- function(assignment_mat, cov_latent){
  assignment_mat %*% cov_latent %*% t(assignment_mat)
}

.estimate_cluster <- function(partition_mat, K){
  stopifnot(nrow(partition_mat) == ncol(partition_mat))

  res <- LICORS::kmeanspp(partition_mat, K)
  res$cluster
}
