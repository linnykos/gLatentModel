source("simulation_header.R")

trials <- 100
paramMat <- cbind(4, 6, ceiling(exp(seq(log(5), log(500), length.out = 30))), 1)
colnames(paramMat) <- c("K", "times_dim", "n", "noise")

rule_closure <- function(){
  function(vec){
    K <- vec[1]; n <- vec[3]; times <- vec[2]
    L <- huge::huge.generator(n = n, d = K, graph = "hub", g = 3, verbose = F)
    latent_dat <- L$data

    a_mat <- diag(K)
    for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
    true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

    dat <- latent_dat%*%t(a_mat)
    dat <- dat + vec[4]*rnorm(prod(dim(dat)))

    res <- gLatentModel(dat, K, seed = 10, num_subsample = 50)
    res <- gLatentModel:::.reshuffle(res, true_cluster)
 
    c(.max_norm_mat(res$cov_latent, L$sigma), .forbenius_norm_mat(res$cov_latent, L$sigma),
      .spectral_norm_mat(res$cov_latent, L$sigma), .L1_norm_mat(res$cov_latent, L$sigma),
      .cluster_distance(res$cluster, true_cluster))
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
      .cluster_distance(res$cluster, true_cluster))
  }
}

#################################################

rule_glatent <- rule_closure()
rule_hclust <- rule_closure_naive(method = naive_clustering_hclust)
rule_sbm <- rule_closure_naive(method = naive_clustering_sbm)
criterion <- function(x, vec){x}

glatent_res <- simulationGenerator(rule_glatent, paramMat, criterion,
  trials, 20)
save.image("results.RData", safe = F)

hclust_res <- simulationGenerator(rule_hclust, paramMat, criterion,
  trials, 20)
save.image("results.RData", safe = F)

sbm_res <- simulationGenerator(rule_sbm, paramMat, criterion,
  trials, 20)
save.image("results.RData", safe = F)
