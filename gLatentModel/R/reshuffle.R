.reshuffle <- function(gmodel, a_mat){
  stopifnot(all(dim(a_mat) == dim(gmodel$assignment_mat)))

  a_mat_est <- gmodel$assignment_mat; K <- ncol(a_mat)

  #each column is estimate, each row is truth
  dif_mat <- apply(a_mat_est, 2, function(x){
    apply(a_mat, 2, function(y){.l2norm(x-y)})
  })

  idx_pairing <- matrix(0, 2, K)
  unmatched_rows <- 1:K; unmatched_cols <- 1:K

  for(i in 1:K){
    idx_pairing[,i] <- which(dif_mat == min(dif_mat[unmatched_rows, unmatched_cols]),
      arr.ind = T)[1,]
    unmatched_rows <- unmatched_rows[unmatched_rows != idx_pairing[1,i]]
    unmatched_cols <- unmatched_cols[unmatched_cols != idx_pairing[2,i]]
  }

  idx_pairing <- idx_pairing[,order(idx_pairing[1,])]
  gmodel$cov_latent <- gmodel$cov_latent[idx_pairing[2,] ,idx_pairing[2,]]
  gmodel$assignment_mat <- gmodel$assignment_mat[,idx_pairing[2,]]
  gmodel$partition_list <- gmodel$partition_list[idx_pairing[2,]]
  gmodel$cluster <- plyr::mapvalues(gmodel$cluster, 1:4, idx_pairing[2,])

  gmodel
}


.reshuffle_cluster <- function(obj, true_cluster){
  stopifnot(class(obj) == "gLatentModel", length(true_cluster) == nrow(obj$assignment_mat))

  K <- max(true_cluster)
  tab <- table(obj$cluster, true_cluster)
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
  obj$cluster <- plyr::mapvalues(obj$cluster, 1:4, idx_pairing[2,])

  obj
}
