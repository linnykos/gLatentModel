pure_nodes <- function(c_mat, K){
  res <- stats::kmeans(c_mat, K)

  idx <- sapply(1:K, function(x){
    loc <- which(res$cluster == x)
    loc[which.max(diag(c_mat)[loc])]
  })

  sort(idx)
}
