context("Test reshuffle")

# .reshuffle is correct

test_that(".reshuffle works", {
  set.seed(10)
  K <- 4; n <- 100; times <- 3
  L <- huge::huge.generator(n = n, d = K, graph = "hub", g = 3, verbose = F)
  latent_dat <- L$data

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K)

  cluster_vec2 <- plyr::mapvalues(res$cluster, from = 1:4, to = c(4,2,1,3))
  res2 <- res
  res2$cluster <- cluster_vec2; res2$cov_latent <- res2$cov_latent[c(4,2,1,3), c(4,2,1,3)]

  res3 <- .reshuffle(res2, res$cluster)

  expect_true(all(res$cluster == res3$cluster))
  expect_true(any(res$cluster != res2$cluster))

  expect_true(sum(abs(res$cov_latent - res3$cov_latent))/K^2 < 1e-4)
  expect_true(sum(abs(res$cov_latent - res2$cov_latent))/K^2 > 1e-4)
})

test_that(".reshuffle works with a different K", {
  set.seed(10)
  K <- 6; n <- 100; times <- 3
  L <- huge::huge.generator(n = n, d = K, graph = "hub", g = 3, verbose = F)
  latent_dat <- L$data

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + 0.1*rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K)

  cluster_vec2 <- plyr::mapvalues(res$cluster, from = 1:K, to = K:1)
  res2 <- res
  res2$cluster <- cluster_vec2
  res2$cov_latent <- res2$cov_latent[c(K:1), c(K:1)]

  res3 <- .reshuffle(res2, res$cluster)

  expect_true(all(res$cluster == res3$cluster))
  expect_true(any(res$cluster != res2$cluster))

  expect_true(sum(abs(res$cov_latent - res3$cov_latent))/K^2 < 1e-4)
  expect_true(sum(abs(res$cov_latent - res2$cov_latent))/K^2 > 1e-4)
})

test_that(".reshuffle makes the right permutation", {
  vec <- c(4, 6, 500, 1); set.seed(1)

  #run the simulation
  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- huge::huge.generator(n, d = K, graph = "scale-free", verbose = F)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L$omega)

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
  a_mat <- a_mat[,K:1]
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + vec[4]*rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K)
  res <- gLatentModel:::.reshuffle(res, true_cluster)

  #ensure that the given covariance arrangment is the best one possible
  perm <- combinat::permn(1:4)

  .forbenius_norm_mat <- function(mat1, mat2){
    stopifnot(all(dim(mat1) == dim(mat2)))
    sum((mat1 - mat2)^2)/prod(dim(mat1))
  }

  error_vec <- sapply(perm, function(x){
    tmp <- res$cov_latent[x,x]
    .forbenius_norm_mat(tmp, L$omega)
  })

  expect_true(which.min(error_vec) == 1)
})
