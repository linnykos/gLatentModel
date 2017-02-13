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
      res <- stats::kmeans(dat[,idx,drop = F], 2, nstart = 10)
      res$betweenss/res$totss
    })

    x <- which.max(val)
    res <- stats::kmeans(t(dat[,which(clust == x),drop = F]), 2, nstart = 10)$cluster
    idx_split <- which(clust == x)[which(res == 2)]

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
      res <- stats::kmeans(t(dat[,idx]), 2, nstart = 10)
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
