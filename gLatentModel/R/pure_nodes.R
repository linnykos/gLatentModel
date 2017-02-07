pure_nodes <- function(cov_denoise, K){
  res <- stats::kmeans(scale(cov_denoise), K)

  idx <- sapply(1:K, function(x){
    loc <- which(res$cluster == x)
    loc[which.max(diag(cov_denoise)[loc])]
  })

  sort(idx)
}

