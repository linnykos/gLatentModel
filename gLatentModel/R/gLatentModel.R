#' Estimate G Latent Model
#'
#' @param dat n by d matrix where there are n samples and d variables
#' @param K number of clusters for the d variables
#' @param ... additional parameters for binary search or cord
#'
#' @return a gLatentModel object
#' @export
gLatentModel <- function(dat, K, ...){
  n <- nrow(dat); d <- ncol(dat); iter <- 1

  clust <- .cord_binary_search(dat, K, ...)
  clust <- .adjustment_cluster(clust, K, dat)

  dat_avg <- .average_data(dat, clust)
  cov_latent <- stats::cov(dat_avg)

  structure(list(cov_latent = cov_latent, cluster = clust, method = "cord"),
            class = "gLatentModel")
}

.cord_binary_search <- function(dat, K, max_iter = 5, low_coef = 0.5, high_coef = 3.5, ...){
  iter <- 1

  while(T){
    mid_coef <- mean(c(low_coef, high_coef))
    res <- cord::cord(dat, tau = mid_coef*sqrt(log(ncol(dat))/nrow(dat)), ...)$cluster

    iter <- iter + 1

    if(length(unique(res)) > K & iter <= max_iter){
      low_coef <- mid_coef
    } else if(length(unique(res)) < K & iter <= max_iter){
      high_coef <- mid_coef
    } else {
      break()
    }
  }

  res
}

.average_data <- function(dat, idx){
  K <- max(idx)

  sapply(1:K, function(x){
    apply(dat[,which(idx == x),drop = F], 1, mean)
  })
}
