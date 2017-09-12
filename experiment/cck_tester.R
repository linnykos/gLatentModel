vec <- rep(NA, 500)

for(i in 1:500){
  if(i %% floor(500/10) == 0) cat('*')

  dat <- MASS::mvrnorm(500, mu = rep(0,4), Sigma = diag(4))
  g <- row_difference_closure(1,4,4)
  res <- cck(dat, g = g)
  vec[i] <- res$pval
}
