source("../main/simulation_header.R")

trials <- 500; d <- 16; n <- 500
res <- vector("list", trials)

cov_mat <- matrix(.2, d, d)
diag(cov_mat) <- 1
cov_mat <- clean_matrix(cov_mat)

for(i in 1:trials){
  set.seed(10*i)
  if(i %% floor(trials/10) == 0) cat('*')

  dat <- MASS::mvrnorm(n, mu = rep(0,16), Sigma = cov_mat)

  combn_mat <- combn(d, 2)
  g_list <- lapply(1:ncol(combn_mat), function(x){
    row_difference_closure(combn_mat[1,x], combn_mat[2,x], d)})

  res[[i]] <- stepdown(dat, g_list, cores = 20)
}

percentage <- length(unlist(res))/(d*(d-1)/2*trials)
percentage
save.image("null_stepdown.RData")
