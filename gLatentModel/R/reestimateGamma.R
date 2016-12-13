reestimate_gamma <- function(dat, partition_list){
  if(any(sapply(partition_list, length) == 1)){
    print(partition_list)
    stop("Singular partitions detected, causing errors")
  }

  n <- nrow(dat); d <- ncol(dat); K <- length(partition_list)

  var_vec <- sapply(partition_list, function(x){
    var_vec <- apply(dat[x,], 2, var)
    sum(var_vec)/(n * length(x))
  })

  idx <- rep(0, d)
  for(i in 1:K){
    idx[partition_list[[i]]] <- i
  }

  var_vec[idx]
}
