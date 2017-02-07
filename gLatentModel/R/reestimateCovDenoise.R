.reestimate_cov_denoise <- function(dat, cluster_vec, gamma_vec = NA){
  if(any(table(cluster_vec) == 1)){
    warning("Singular partitions detected, returning original gamma_vec")
    if(any(is.na(gamma_vec))) stop() else return(gamma_vec)
  }

  n <- nrow(dat); d <- ncol(dat); K <- max(cluster_vec)

  var_vec <- sapply(1:K, function(x){
    idx <- which(cluster_vec == x)
    var_vec <- apply(dat[,idx, drop = F], 2, stats::var)
    sum(var_vec)/(n * length(idx))
  })

  var_vec[cluster_vec]
}
