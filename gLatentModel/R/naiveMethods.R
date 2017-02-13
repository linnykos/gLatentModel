naive_clustering_hclust <- function(dat, K){
  res <- stats::hclust(dist(t(dat)))
  cluster <- cutree(res, K)

  theta_mat <- stats::cov(.average_data(dat, cluster))

  list(cov_latent = theta_mat, cluster = cluster)
}

# motivation from https://projecteuclid.org/euclid.aos/1418135620
naive_clustering_sbm <- function(dat, K){
  d <- ncol(dat)

  cor_est <- stats::cor(dat)

  #compute confidence interval
  conf_mat <- matrix(psych::r.con(cor_est,n), ncol=2)
  idx1 <- which(apply(conf_mat, 1, function(x){
    if(any(is.nan(x))) return(T)
    if(x[1]<0 & x[2]>0) return(F) else (T)
  }))

  g_est <- matrix(0, d, d)
  g_est[idx1] <- 1

  eig <- eigen(adj_mat)
  new_dat <- eig$vectors[,1:K]

  res <- stats::kmeans(new_dat, K)

  theta_mat <- stats::cov(.average_data(dat, res$cluster))

  list(cov_latent = theta_mat, cluster = res$cluster)
}
