source("simulation_header.R")

trials <- 100
paramMat <- cbind(6, ceiling(exp(seq(log(5), log(500), length.out = 30))), 1)
colnames(paramMat) <- c("K", "times_dim", "n", "noise")

rule_closure <- function(){
  function(vec){
    K <- vec[1]; n <- vec[2]
    L <- huge::huge.generator(n = 5, d = K, graph = "scale-free", verbose = F)
    latent_dat <- MASS::mvrnorm(n, rep(0,K), L$omega)

    block_idx <- c(0,sapply(2:(K+1), function(x){sum(1:x)}))
    a_mat <- sapply(2:(K+1), function(i){
      vec <- rep(0, nrow(cov_mat))
      vec[(block_idx[i-1]+1):block_idx[i]] <- 1
      vec
    })
    true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

    dat <- latent_dat%*%t(a_mat)
    dat <- dat + vec[3]*rnorm(prod(dim(dat)))

    res <- gLatentModel(dat, K, seed = 10, num_subsample = 50)
    res <- gLatentModel:::.reshuffle(res, true_cluster)

    c(.max_norm_mat(res$cov_latent, L$sigma), .forbenius_norm_mat(res$cov_latent, L$sigma),
      .spectral_norm_mat(res$cov_latent, L$sigma), .L1_norm_mat(res$cov_latent, L$sigma),
      .cluster_distance(res$cluster, true_cluster),
      .jaccard_distance(res$cluster, true_cluster))
  }
}

rule_closure_naive <- function(method = naive_clustering_hclust){
  function(vec){
    K <- vec[1]; n <- vec[3]; times <- vec[2]
    L <- huge::huge.generator(n = n, d = K, graph = "hub", g = 3, verbose = F)
    latent_dat <- L$data

    a_mat <- diag(K)
    for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
    true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

    dat <- latent_dat%*%t(a_mat)
    dat <- dat + vec[4]*rnorm(prod(dim(dat)))

    res <- method(dat, K)
    res <- gLatentModel:::.reshuffle(res, true_cluster)

    c(.max_norm_mat(res$theta, L$sigma), .forbenius_norm_mat(res$theta, L$sigma),
      .spectral_norm_mat(res$theta, L$sigma), .L1_norm_mat(res$theta, L$sigma),
      .cluster_distance(res$cluster, true_cluster),
      .jaccard_distance(res$cluster, true_cluster))
  }
}
