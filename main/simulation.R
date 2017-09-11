source("../main/simulation_header.R")

trials <- 100
d <- 40
n_seq <- c(10, 40, 100, 250)
strength_seq <- c(0, 0.3, 0.6, 0.9)
param_mat <- as.matrix(expand.grid(n_seq, strength_seq))

rule_closure <- function(d, bootstrap_trials = 200){
  function(vec){
    cov_mat <- generate_matrix(d = d, strength = vec[2])
    cov_mat <- clean_matrix(cov_mat)

    dat <- MASS::mvrnorm(vec[1], mu = rep(0,d), Sigma = cov_mat)

    ## use cck
    combn_mat <- combn(d, 2)
    g_list <- lapply(1:ncol(combn_mat), function(x){
      gLatentModel::row_difference_closure(combn_mat[1,x], combn_mat[2,x], d)})
    cck_idx <- gLatentModel::stepdown(dat, g_list, cores = 4, alpha = 0.05, trials = bootstrap_trials)
    cck_res <- gLatentModel::connected_components(d, combn_mat[,cck_idx])

    cck_idx2 <- gLatentModel::stepdown(dat, g_list, cores = 4, alpha = 0.2, trials = bootstrap_trials)
    cck_res2 <- gLatentModel::connected_components(d, combn_mat[,cck_idx2])

    mcord_res <- cord::cord(dat)$cluster
    mcord_res <- lapply(unique(mcord_res), function(x){which(mcord_res == x)})

    dd <- as.dist((1 - cor(dat))/2)
    obj <- hclust(dd, method = "single")
    hier_res <- cutree(obj, 4)
    hier_res <- lapply(unique(hier_res), function(x){which(hier_res == x)})

    list(cck_res = cck_res, mcord_res = mcord_res, hier_res = hier_res)
  }
}

rule <- rule_closure(d)
criterion <- function(x, vec){x}

res <- simulationGenerator(rule, param_mat, criterion, trials, cores = NA)

save.image("../main/results.RData")

