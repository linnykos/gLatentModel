vec <- rep(NA, 500)

for(i in 1:500){
  if(i %% floor(500/10) == 0) cat('*')

  dat <- MASS::mvrnorm(500, mu = rep(0,4), Sigma = diag(4))
  g <- row_difference_closure(1,4,4)
  res <- cck(dat, g = g)
  vec[i] <- res$pval
}

############

cov_mat <- diag(4)
cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5
mean_vec <- rep(0,4) #c(1:4)

trials <- 1000
vec <- rep(NA, trials)
for(i in 1:trials){
  if(i %% floor(trials/10) == 0) cat('*')
  set.seed(i)

  dat <- MASS::mvrnorm(500, mu = c(1:4), Sigma = cov_mat)
  g <- row_difference_closure(1,4,4)
  dat <- scale(dat, scale = F)
  sigma_vec <- apply(dat, 2, stats::sd)
  psi <- cor_vec(dat, sigma_vec = sigma_vec)
  t0 <- max(abs(g(psi)))
  vec[i] <- t0
}
hist(vec)

######################

trials <- 1000; n <- 500
cov_mat <- diag(4)
cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5
mean_vec <- rep(0,4) #c(1:4)

set.seed(1)
dat <- MASS::mvrnorm(n, mu = mean_vec, Sigma = cov_mat)
g <- row_difference_closure(1,4,4)
dat <- scale(dat, scale = F)
sigma_vec <- apply(dat, 2, stats::sd)
psi <- cor_vec(dat, sigma_vec = sigma_vec)
theta <- g(psi)
t0 <- max(abs(theta))

doMC::registerDoMC(cores = 2)

func <- function(i){
  idx <- sample(1:n, n, replace = T)
  dat_tmp <- dat[idx,]
  dat_tmp <- scale(dat_tmp, scale = F)
  sigma_tmp <- apply(dat_tmp, 2, stats::sd)
  psi_tmp <- cor_vec(dat_tmp, sigma_vec = sigma_tmp)
  g(psi_tmp)
}

i <- 1 #debugging purposes
theta_boot <- foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                                        func(i))
t_boot <- sapply(theta_boot, function(x){max(abs(x - theta))})
hist(t_boot)

####################

minval <- min(c(t_boot, vec)); maxval <- max(c(t_boot, vec))
plot(sort(t_boot), sort(vec), asp = T, xlim = c(minval, maxval), ylim = c(minval, maxval))
lines(c(minval, maxval), c(minval, maxval), lwd = 2, col = "red")
