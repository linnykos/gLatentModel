reestimate_gamma <- function(dat, partition_list, gamma_vec = NA){
  if(any(sapply(partition_list, length) == 1)){
    warning("Singular partitions detected, returning original gamma_vec")
    if(any(is.na(gamma_vec))) stop() else return(gamma_vec)
  }

  n <- nrow(dat); d <- ncol(dat); K <- length(partition_list)

  var_vec <- sapply(partition_list, function(x){
    var_vec <- apply(dat[,x], 2, stats::var)
    sum(var_vec)/(n * length(x))
  })

  idx <- rep(0, d)
  for(i in 1:K){
    idx[partition_list[[i]]] <- i
  }

  var_vec[idx]
}
