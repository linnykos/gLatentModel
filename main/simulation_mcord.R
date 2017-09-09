source("../main/simulation_helper.R")
strength_vec <- c(0, 0.3, 0.6, 0.9)
trials <- 10
clust_list <- lapply(1:length(strength_vec), function(x){rep(NA, trials)})
error_list <- lapply(1:length(strength_vec), function(x){rep(NA, trials)})

n <- 500; d <- 100
true_cluster <- rep(1:4, each = d/4)

for(i in 1:length(strength_vec)){
  cov_mat <- generate_matrix(n = n, d = d, strength = strength_vec[i])
  cov_mat <- clean_matrix(cov_mat)

  for(j in 1:trials){
    set.seed(10*j)
    dat <- MASS::mvrnorm(n, mu = rep(0,d), Sigma = cov_mat)
    res <- cord::cord(dat)
    clust_list[[i]][j] <- length(unique(res$cluster))
    error_list[[i]][j] <- error_function(true_cluster, res$cluster)

    print('*')
  }
  print(i)
}

sapply(clust_list, mean)
sapply(error_list, mean)


