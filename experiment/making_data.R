library(huge)

set.seed(10)
K <- 4; n <- 1000
L <- huge.generator(n = n, d = K, graph = "hub", g = 3)
latent_dat <- L$data

a_mat <- rbind(diag(K), diag(K), diag(K), diag(K), diag(K), diag(K))
dat <- latent_dat%*%t(a_mat)
dat <- dat + rnorm(prod(dim(dat)))
#dat <- dat[sample(1:n), sample(1:(K*6))]

#image(cor(dat))

res <- gLatentModel(dat, K, verbose = T)
res <- reshuffle(res, a_mat)

res$theta ## theta might be in a different order...
L$sigma
cov(latent_dat)
