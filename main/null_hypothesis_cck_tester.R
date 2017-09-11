source("../main/simulation_header.R")

trials <- 1000
vec <- rep(NA, trials)

cov_mat <- matrix(0, 16, 16)
cov_mat[,1:4] <- rep(c(0, 0.25, 0.5, 0.75), each = 4)
cov_mat[,5:8] <- rep(c(0.25, 0.25, 0.5, 0.75), each = 4)
cov_mat[,9:12] <- rep(c(0.5, 0.5, 0.5, 0.75), each = 4)
cov_mat[,13:16] <- 0.75

cov_mat <- clean_matrix(cov_mat)

for(i in 1:trials){
  if(i %% floor(trials/10) == 0) cat('*')

  dat <- MASS::mvrnorm(1000, mu = rep(0,16), Sigma = cov_mat)
  g <- gLatentModel::row_difference_closure(1,2,16)
  vec[i] <- gLatentModel::cck(dat, g = g, trials = 500, cores = 20)$pval
}

save.image("null_hypothesis.RData")
