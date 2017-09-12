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
mean_vec <- c(1:4)

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
mean_vec <- c(1:4)

set.seed(1)
dat <- MASS::mvrnorm(n, mu = mean_vec, Sigma = cov_mat)
g <- row_difference_closure(1,4,4)
dat <- scale(dat, scale = F)
sigma_vec <- apply(dat, 2, stats::sd)
psi <- cor_vec(dat, sigma_vec = sigma_vec)
t0 <- max(abs(g(psi)))

doMC::registerDoMC(cores = 2)

func <- function(i){
  e <- stats::rnorm(n)
  g(cor_vec(dat, sigma_vec = sigma_vec, noise_vec = e),
    average_vec = psi*sum(e)/n)
}

i <- 1 #debugging purposes
t_boot <- as.numeric(foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                                        func(i)))
hist(t_boot)

####################

minval <- min(c(t_boot, vec)); maxval <- max(c(t_boot, vec))
plot(sort(t_boot), sort(vec), asp = T, xlim = c(minval, maxval), ylim = c(minval, maxval))
lines(c(minval, maxval), c(minval, maxval), lwd = 2, col = "red")
