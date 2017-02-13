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

.adjustment_cluster <- function(clust, K, dat = NA){
  while(T){
    if(length(unique(clust)) > K){
      clust <- .adjustment_cluster_remove(clust, dat = NA)
    } else if(length(unique(clust)) < K){
      clust <- .adjustment_cluster_add(clust, dat = NA)
    } else{
      break()
    }
  }

  clust
}

.adjustment_cluster_add <- function(clust, dat = NA){
  if(all(!is.na(dat))){
    val <- sapply(1:max(clust), function(x){
      idx <- which(clust == x)
      res <- stats::kmeans(dat[,idx,dropped = F], 2, nstart = 10)
      res$betweenss/res$totss
    })

    x <- which.max(val)
    res <- stats::kmeans(dat[,which(clust == x),dropped = F], 2, nstart = 10)
    res <- res$cluster
    idx_split <- which(cluster == x)[which(res == 2)]

  } else {
    tab <- table(clust)
    idx <- which.max(tab)
    val <- as.numeric(names(tab)[idx])

    idx <- which(clust == val)
    idx_split <- sample(idx, floor(length(idx)/2))
  }

  clust[idx_split] <- length(unique(clust))+1
  clust
}

.adjustment_cluster_remove <- function(clust, dat = NA){
  if(all(!is.na(dat))){
    combin_mat <- utils::combn(max(clust), 2)
    val <- apply(combin_mat, 2, function(x){
      idx <- which(clust %in% x)
      res <- stats::kmeans(dat[,idx], 2, nstart = 10)
      res$betweenss/res$totss
    })

    x <- which.min(val)
    clust[which(clust %in% combin_mat[,x])] <- combin_mat[1,x]

  } else {
    tab <- table(clust)
    idx <- order(tab, decreasing = F)[1:2]
    val <- as.numeric(names(tab)[idx])

    idx_merge <- which(clust %in% val)
    clust[idx_merge] <- val[1]
  }

  #reindex
  uniq <- unique(clust)
  clust2 <- clust
  for(i in 1:length(uniq)){
    clust2[clust == uniq[i]] <- i
  }

  clust2
}

.average_data <- function(dat, idx){
  K <- max(idx)

  sapply(1:K, function(x){
    apply(dat[,which(idx == x),drop = F], 1, mean)
  })
}
