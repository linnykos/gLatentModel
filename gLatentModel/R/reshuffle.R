reshuffle <- function(gmodel, a_mat){
  stopifnot(all(dim(a_mat) == dim(gmodel$a)))

  a_mat_est <- gmodel$a; K <- ncol(a_mat)

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
  gmodel$theta <- gmodel$theta[idx_pairing[2,] ,idx_pairing[2,]]
  gmodel$a <- gmodel$a[,idx_pairing[2,]]
  gmodel$partition_list <- gmodel$partition_list[idx_pairing[2,]]

  gmodel
}
