.reestimate_assignment <- function(cluster_vec){
  K <- max(cluster_vec)
  sapply(1:K, function(x){
    as.numeric(cluster_vec == x)
  })
}
