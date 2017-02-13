.average_data <- function(dat, idx){
  K <- max(idx)

  sapply(1:K, function(x){
    apply(dat[,which(idx == x),drop = F], 1, mean)
  })
}

.cluster_distance <- function(vec1, vec2, samples = 500){
  stopifnot(length(vec1) == length(vec2))

  len <- length(vec1)

  if(len*(len-1)/2 > samples){
    cluster <- sapply(1:samples, function(x){sample(1:len, 2)})
  } else {
    cluster <- combn(len, 2)
  }

  res <- apply(cluster, 2, function(x){
    bool1 <- (vec1[x[1]] == vec1[x[2]])
    bool2 <- (vec2[x[1]] == vec2[x[2]])
    if(bool1 == bool2) return(0) else return(1)
  })

  sum(res)/ncol(cluster)
}

.jaccard_distance <- function(vec1, vec2){
  stopifnot(all(sort(unique(vec1)) == sort(unique(vec2))))

  res <- sapply(unique(vec1), function(x) {
    idx1 <- which(vec1 == x); idx2 <- which(vec2 == x)
    1-length(intersect(idx1, idx2))/length(unique(c(idx1, idx2)))
  })

  median(res)
}

.graph_edge_distance <- function(cov_est, cov_true, n, tol = 1e-4){
  stopifnot(nrow(cov_est) == ncol(cov_est), nrow(cov_true) == ncol(cov_true),
            nrow(cov_est) == nrow(cov_true))

  d <- nrow(cov_est)
  cor_est <- stats::cov2cor(cov_est)

  #compute confidence interval
  conf_mat <- matrix(psych::r.con(cor_est,n),ncol=2)
  idx1 <- which(apply(conf_mat, 1, function(x){
    if(any(is.nan(x))) return(T)
    if(x[1]<0 & x[2]>0) return(F) else (T)
  }))

  g_est <- matrix(0, d, d)
  g_est[idx1] <- 1

  #compute true graph
  idx2 <- which(abs(cov_true) > tol)
  g_true <- matrix(0, d, d)
  g_true[idx2] <- 1

  sum(abs(g_est - g_true))/d^2
}

.max_norm_mat <- function(mat1, mat2){
  stopifnot(all(dim(mat1) == dim(mat2)))
  max(abs(mat1 - mat2))
}

.forbenius_norm_mat <- function(mat1, mat2){
  stopifnot(all(dim(mat1) == dim(mat2)))
  sum((mat1 - mat2)^2)/prod(dim(mat1))
}

.spectral_norm_mat <- function(mat1, mat2){
  stopifnot(all(dim(mat1) == dim(mat2)))
  max(svd(mat1 - mat2)$d)
}

.L1_norm_mat <- function(mat1, mat2){
  stopifnot(all(dim(mat1) == dim(mat2)))
  sum(abs(mat1 - mat2))/prod(dim(mat1))
}
