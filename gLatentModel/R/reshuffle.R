.reshuffle <- function(obj, cluster_vec){
  stopifnot(length(cluster_vec) == nrow(obj$assignment_mat))

  K <- max(cluster_vec)
  tab <- table(obj$cluster, cluster_vec)
  tab <- tab + 0.01*stats::rnorm(K^2)

  idx_pairing <- matrix(0, 2, K)
  unmatched_rows <- 1:K; unmatched_cols <- 1:K

  for(i in 1:K){
    idx_pairing[,i] <- which(tab == max(tab[unmatched_rows, unmatched_cols]),
                             arr.ind = T)[1,]
    unmatched_rows <- unmatched_rows[unmatched_rows != idx_pairing[1,i]]
    unmatched_cols <- unmatched_cols[unmatched_cols != idx_pairing[2,i]]
  }

  idx_pairing <- idx_pairing[,order(idx_pairing[1,])]
  obj$cov_latent <- obj$cov_latent[idx_pairing[2,] ,idx_pairing[2,]]
  obj$cluster <- plyr::mapvalues(obj$cluster, 1:K, idx_pairing[2,])

  obj
}
