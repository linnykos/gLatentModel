library(huge)

set.seed(10)
K <- 4; n <- 100; times <- 6
L <- huge.generator(n = n, d = K, graph = "hub", g = 3)
latent_dat <- L$data

a_mat <- diag(K)
for(i in 1:(times-1)){
  a_mat <- rbind(a_mat, diag(K))
}

dat <- latent_dat%*%t(a_mat)
dat <- dat + rnorm(prod(dim(dat)))
idx <- sample(1:ncol(dat))
dat <- dat[,idx]; a_mat <- a_mat[idx,]

#image(cor(dat))

res <- gLatentModel(dat, K, verbose = T, seed = 10, num_subsample = 10)
res <- reshuffle(res, a_mat)

res$theta
L$sigma
cov(latent_dat)
