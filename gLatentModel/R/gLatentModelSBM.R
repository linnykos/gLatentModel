#' Estimate G Latent Model
#'
#' @param dat n by d matrix where there are n samples and d variables
#' @param K number of clusters for the d variables
#' @param debugging boolean on whether or not to pass additional parameters out
#' @param ... additional parameters for binary search or cord
#'
#' @return a gLatentModel object
#' @export
gLatentModel_sbm <- function(dat, K, debugging = F, ...){
  n <- nrow(dat); d <- ncol(dat)

  cov_mat <- stats::cov(dat)
  latent_cov_vec <- .estimate_cov_noise(dat, ...)

  cov_denoise_mat <- cov_mat - diag(latent_cov_vec)

  eigenvec <- eigen(cov_denoise_mat)$vectors[,1:K]

  #sort eigenvec to ensure reproducibility
  idx <- order(eigenvec[,1])
  eigenvec <- eigenvec[idx,]
  eigenvec <- eigenvec[,order(eigenvec[1,])]

  res <- stats::kmeans(eigenvec, K, nstart = 10)
  clust <- res$cluster
  clust <- clust[order(idx)] #reverse the ordering

  dat_avg <- .average_data(dat, clust)
  cov_latent <- stats::cov(dat_avg)

  if(!debugging){
    structure(list(cov_latent = cov_latent, cluster = clust,
                   method = "sbm"), class = "gLatentModel")
  } else {
    structure(list(cov_latent = cov_latent, cluster = clust,
                   eigenvec = eigenvec, latent_cov_vec = latent_cov_vec,
                   method = "sbm"), class = "gLatentModel")
  }
}
