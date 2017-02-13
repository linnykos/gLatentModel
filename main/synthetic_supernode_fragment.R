source("../main/simulation_header.R")

trials <- 100
paramMat <- cbind(8, 3, ceiling(exp(seq(log(5), log(500), length.out = 30))), 1)
colnames(paramMat) <- c("K", "times_dim", "n", "noise")

rule_closure <- function(method = gLatentModel, ...){
  function(vec){
    K <- vec[1]; n <- vec[3]; times <- vec[2]
    L <- huge::huge.generator(n, d = K, graph = "scale-free", verbose = F)
    latent_dat <- MASS::mvrnorm(n, rep(0,K), L$omega)

    a_mat <- diag(K)
    for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
    true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})
    a_mat <- a_mat[,K:1]

    dat <- latent_dat%*%t(a_mat)
    dat <- dat + vec[4]*rnorm(prod(dim(dat)))

    res <- method(dat, K, ...)
    res <- gLatentModel:::.reshuffle(res, true_cluster)

    c(.max_norm_mat(res$cov_latent, L$omega), .forbenius_norm_mat(res$cov_latent, L$omega),
      .spectral_norm_mat(res$cov_latent, L$omega), .L1_norm_mat(res$cov_latent, L$omega),
      .graph_edge_distance(res$cov_latent, L$omega, n),
      .cluster_distance(res$cluster, true_cluster),
      .jaccard_distance(res$cluster, true_cluster),
      sum(abs(dat))) #the last entry for the signature
  }
}

#################################################

rule_glatent <- rule_closure(method = gLatentModel)
rule_glatent_sbm <- rule_closure(method = gLatentModel_sbm)
rule_hclust <- rule_closure(method = naive_clustering_hclust)
rule_sbm <- rule_closure(method = naive_clustering_sbm)
criterion <- function(x, vec){x}

glatent_res <- simulationGenerator(rule_glatent, paramMat, criterion,
                                   trials, 20)
glatent_sbm_res <- simulationGenerator(rule_glatent_sbm, paramMat, criterion,
                                       trials, 20)
hclust_res <- simulationGenerator(rule_hclust, paramMat, criterion,
                                  trials, 20)
sbm_res <- simulationGenerator(rule_sbm, paramMat, criterion,
                               trials, 20)

save.image("../results/results_fragment.RData", safe = F)
