context("Test gLatentModel_sbm")

## gLatentModel_sbm is correct

test_that("gLatentModel_sbm returns properly", {
  set.seed(10)
  K <- 2
  latent_dat <- MASS::mvrnorm(1000, rep(0, K), matrix(c(2,-1,-1,2),2,2))
  a_mat <- rbind(diag(K), diag(K), diag(K), diag(K), diag(K), diag(K))
  dat <- latent_dat%*%t(a_mat)
  dat <- dat + 3*rnorm(prod(dim(dat)))

  res <- gLatentModel_sbm(dat, K)

  expect_true(class(res) == "gLatentModel")
  expect_true(all(dim(res$cov_latent) == c(2,2)))
  expect_true(is.numeric(res$cov_latent))
  expect_true(is.matrix(res$cov_latent))
})

test_that("gLatentModel_sbm is unaffected (after reshuffling) by the initial order of columns",{
  set.seed(10)
  K <- 4; n <- 100; times <- 3
  L <- huge::huge.generator(n = n, d = K, graph = "hub", g = 3, verbose = F)
  latent_dat <- L$data

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
  cluster_vec <- rep(c(1:K), times = times)

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + rnorm(prod(dim(dat)))

  idx <- sample(1:ncol(dat))
  dat2 <- dat[,idx]

  set.seed(10); res <- gLatentModel_sbm(dat, K, debugging = T)
  set.seed(10); res2 <- gLatentModel_sbm(dat2, K, debugging = T)

  expect_true(all(sort(res$latent_cov_vec) == sort(res2$latent_cov_vec)))
  expect_true(sum(abs(sort(as.numeric(res$cov_latent)) - sort(as.numeric(res2$cov_latent))))/K^2 < 1e-2)
})

###################################

test_that("gLatentModel_sbm test case with no noise",{
  vec <- c(4, 6, 100, 0)
  set.seed(2)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- huge::huge.generator(n, d = K, graph = "scale-free", verbose = F)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L$omega)

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
  a_mat <- a_mat[,K:1]
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})


  dat <- latent_dat%*%t(a_mat)
  dat <- dat + vec[4]*rnorm(prod(dim(dat)))

  res <- gLatentModel_sbm(dat, K)
  res <- gLatentModel:::.reshuffle(res, true_cluster)

  expect_true(all(res$cluster == true_cluster))
})

test_that("gLatentModel_sbm test case with high SNR for even larger n",{
  vec <- c(4, 6, 500, 0.01)
  set.seed(2)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- huge::huge.generator(n, d = K, graph = "scale-free", verbose = F)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L$omega)

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
  a_mat <- a_mat[,K:1]
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})


  dat <- latent_dat%*%t(a_mat)
  dat <- dat + vec[4]*rnorm(prod(dim(dat)))

  res <- gLatentModel_sbm(dat, K)
  res <- gLatentModel:::.reshuffle(res, true_cluster)

  expect_true(all(res$cluster == true_cluster))
})

test_that("gLatentModel_sbm doesn't deterioate with small noise",{
  vec <- c(4, 6, 50, 0.01)
  set.seed(2)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- diag(K)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L)

  a_mat <- sapply(1:K, function(x){
    vec <- rep(0, K*times)
    vec[((x-1)*times+1):(x*times)] <- 1
    vec
  })
  a_mat <- a_mat[,K:1]
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

  dat <- latent_dat%*%t(a_mat)

  set.seed(10); res <- gLatentModel_sbm(dat, K, debugging = T)
  res <- gLatentModel:::.reshuffle(res, true_cluster)

  ###

  dat <- dat + vec[4]*rnorm(prod(dim(dat)))

  set.seed(10); res2 <- gLatentModel_sbm(dat, K, debugging = T)
  res2 <- gLatentModel:::.reshuffle(res2, true_cluster)

  expect_true(all(res$cluster == res2$cluster))
})

test_that("gLatentModel_sbm is invariant to the order of columns",{
  vec <- c(4, 6, 100)
  set.seed(2)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- diag(K)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L)

  a_mat <- sapply(1:K, function(x){
    vec <- rep(0, K*times)
    vec[((x-1)*times+1):(x*times)] <- 1
    vec
  })
  a_mat <- a_mat[,K:1]
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)[1]})

  dat <- latent_dat%*%t(a_mat)

  res <- gLatentModel_sbm(dat, K, debugging = T)
  res <- gLatentModel:::.reshuffle(res, true_cluster)

  #####

  a_mat2 <- diag(K)
  for(i in 1:(times-1)){a_mat2 <- rbind(a_mat2, diag(K))}
  a_mat2 <- a_mat2[,K:1]
  true_cluster2 <- apply(a_mat2, 1, function(x){which(x == 1)[1]})

  dat2 <- latent_dat%*%t(a_mat2)

  res2 <- gLatentModel_sbm(dat2, K, debugging = T)
  res2 <- gLatentModel:::.reshuffle(res2, true_cluster2)

  expect_true(all(sort(res$latent_cov_vec) == sort(res2$latent_cov_vec)))
  expect_true(all(table(res$cluster) == table(res2$cluster)))
})

test_that("gLatentModel_sbm will output the right number of clusters", {
  vec <- c(4, 6, 5, 1)
  set.seed(1)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- huge::huge.generator(n, d = K, graph = "scale-free", verbose = F)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L$omega)

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
  a_mat <- a_mat[,K:1]
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})


  dat <- latent_dat%*%t(a_mat)
  dat <- dat + vec[4]*rnorm(prod(dim(dat)))

  res <- gLatentModel_sbm(dat, K)

  expect_true(length(unique(res$cluster)) == 4)

})
