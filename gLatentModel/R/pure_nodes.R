pure_nodes <- function(cov_denoise, K){
  res <- stats::kmeans(scale(cov_denoise), K)

  #return the node closest to the cluster center
  idx <- sapply(1:K, function(x){
    node <- which(res$cluster == x)
    l2vec <- apply(cov_denoise[,node,drop = F], 2, function(y){
      .l2norm(y - res$centers[x,,drop = F])
    })

    node[which.min(l2vec)]
  })

  sort(idx)
}

