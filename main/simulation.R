rm(list=ls())
devtools::install_github("linnylin92/gLatentModel", ref = "cck",
                         subdir = "gLatentModel", force = T)

source("../main/simulation_helper.R")
d <- 40; n <- 100

cov_mat <- generate_matrix(n = n, d = d, strength = 0)
cov_mat <- clean_matrix(cov_mat)

# cov_mat_simplified <- matrix(sapply(c(0:3)*d4+1, function(x){
#   cov_mat[c(0:3)*d4+2,x]
# }), ncol = 4, nrow = 4)

set.seed(10)
dat <- MASS::mvrnorm(n, mu = rep(0,d), Sigma = cov_mat)

## use cck
combn_mat <- combn(d, 2)
g_list <- lapply(1:ncol(combn_mat), function(x){
  gLatentModel::row_difference_closure(combn_mat[1,x], combn_mat[2,x], d)})
cck_idx <- gLatentModel::stepdown(dat, g_list, cores = 4, alpha = 0.01)
cck_res <- gLatentModel::connected_components(d, combn_mat[,cck_idx])

## use mcord
mcord_res <- cord::cord(dat)

## use hierarchical clust
dd <- as.dist((1 - cor(dat))/2)
obj <- hclust(dd, method = "single")
hier_res <- cutree(obj, 4)

save.image("../main/results.RData")

