source("../main/simulation_header.R")

trials <- 1000
vec <- rep(NA, trials)

for(i in 1:trials){
  set.seed(10*i)
  if(i %% floor(trials/10) == 0) cat('*')

  dat <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = diag(4))
  g <- gLatentModel::row_difference_closure(1,4,4)
  vec[i] <- gLatentModel::cck(dat, g = g, trials = 200, cores = 20)$pval
}

save.image("null_hypothesis_small.RData")

# plot(sort(vec), seq(0, 1, length.out = length(vec)), asp = T)
# lines(c(0,1),c(0,1), lwd = 2, col = "red")
