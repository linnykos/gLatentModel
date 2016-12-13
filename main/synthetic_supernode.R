source("source_header.R")

trials <- 10
paramMat <- cbind(4, 6, 100, seq(0.05, 10, length.out = 5))
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

    res <- gLatentModel(dat, K, seed = 10, num_subsample = 10)
    res <- reshuffle(res, a_mat)

    c(.max_norm_mat(res$theta, L$sigma), .forbenius_norm_mat(res$theta, L$sigma),
      .spectral_norm_mat(res$theta, L$sigma), .L1_norm_mat(res$theta, L$sigma),
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
    res <- reshuffle_cluster(res, true_cluster)

    c(.max_norm_mat(res$theta, L$sigma), .forbenius_norm_mat(res$theta, L$sigma),
      .spectral_norm_mat(res$theta, L$sigma), .L1_norm_mat(res$theta, L$sigma),
      .cluster_distance(res$cluster, true_cluster))
  }
}

rule_closure_oracle <- function(){
  function(vec){
    K <- vec[1]; n <- vec[3]; times <- vec[2]
    L <- huge::huge.generator(n = n, d = K, graph = "hub", g = 3, verbose = F)
    latent_dat <- L$data

    a_mat <- diag(K)
    for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
    true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

    dat <- latent_dat%*%t(a_mat)
    dat <- dat + vec[4]*rnorm(prod(dim(dat)))

    c(.max_norm_mat(L$sigmahat, L$sigma), .forbenius_norm_mat(L$sigmahat, L$sigma),
      .spectral_norm_mat(L$sigmahat, L$sigma), .L1_norm_mat(L$sigmahat, L$sigma))
  }
}

#################################################

rule_glatent <- rule_closure()
rule_hclust <- rule_closure_naive(method = naive_clustering_hclust)
rule_sbm <- rule_closure_naive(method = naive_clustering_sbm)
rule_oracle <- rule_closure_oracle()
criterion <- function(x, vec){x}

glatent_res <- simulationGenerator(rule_glatent, paramMat, criterion,
  trials, NA)
hclust_res <- simulationGenerator(rule_hclust, paramMat, criterion,
  trials, NA)
sbm_res <- simulationGenerator(rule_sbm, paramMat, criterion,
  trials, NA)
oracle_res <- simulationGenerator(rule_oracle, paramMat[1,,drop=F], criterion,
  trials, NA)

save.image("results.RData", safe = F)
