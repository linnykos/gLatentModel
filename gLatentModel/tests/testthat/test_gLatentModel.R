context("Test gLatentModel")

## gLatentModel is correct

test_that("gLatentModel returns properly", {
  set.seed(10)
  K <- 2
  latent_dat <- MASS::mvrnorm(1000, rep(0, K), matrix(c(2,-1,-1,2),2,2))
  a_mat <- rbind(diag(K), diag(K), diag(K), diag(K), diag(K), diag(K))
  dat <- latent_dat%*%t(a_mat)
  dat <- dat + 3*rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K)

  expect_true(class(res) == "gLatentModel")
  expect_true(all(dim(res$cov_latent) == c(2,2)))
  expect_true(is.numeric(res$cov_latent))
  expect_true(is.matrix(res$cov_latent))
})

test_that("gLatentModel is unaffected (after reshuffling) by the initial order of columns",{
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

  res <- gLatentModel(dat, K)
  res2 <- gLatentModel(dat2, K)

  expect_true(sum(abs(sort(as.numeric(res$cov_latent)) - sort(as.numeric(res2$cov_latent))))/K^2 < 1e-2)
})

###################################

test_that("gLatentModel test case with high SNR",{
  vec <- c(4, 6, 100, 0.01)
  set.seed(2)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- huge::huge.generator(n, d = K, graph = "scale-free", verbose = F)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L$omega)

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + vec[4]*rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K)
  res <- gLatentModel:::.reshuffle(res, true_cluster)

  expect_true(all(res$cluster == true_cluster))
})

test_that("gLatentModel test case with high SNR for even larger n",{
  vec <- c(4, 6, 500, 0.01)
  set.seed(2)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- huge::huge.generator(n, d = K, graph = "scale-free", verbose = F)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L$omega)

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + vec[4]*rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K)
  res <- gLatentModel:::.reshuffle(res, true_cluster)

  expect_true(all(res$cluster == true_cluster))
})

test_that("gLatentModel doesn't deterioate with small noise",{
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
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

  dat <- latent_dat%*%t(a_mat)

  res <- gLatentModel(dat, K)
  res <- gLatentModel:::.reshuffle(res, true_cluster)

  ###

  dat <- dat + vec[4]*rnorm(prod(dim(dat)))

  res2 <- gLatentModel(dat, K)
  res2 <- gLatentModel:::.reshuffle(res2, true_cluster)

  expect_true(all(res$cluster == res2$cluster))
})

test_that("gLatentModel is invariant to the order of columns",{
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
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)[1]})

  dat <- latent_dat%*%t(a_mat)

  res <- gLatentModel(dat, K)
  res <- gLatentModel:::.reshuffle(res, true_cluster)

  #####

  a_mat2 <- diag(K)
  for(i in 1:(times-1)){a_mat2 <- rbind(a_mat2, diag(K))}
  true_cluster2 <- apply(a_mat, 1, function(x){which(x == 1)[1]})

  dat2 <- latent_dat%*%t(a_mat2)

  res2 <- gLatentModel(dat2, K)
  res2 <- gLatentModel:::.reshuffle(res2, true_cluster2)

  expect_true(all(table(res$cluster) == table(res2$cluster)))
})

test_that("gLatentModel will output the right number of clusters", {
  vec <- c(4, 6, 5, 1)
  set.seed(1)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- huge::huge.generator(n, d = K, graph = "scale-free", verbose = F)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L$omega)

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}
  true_cluster <- apply(a_mat, 1, function(x){which(x == 1)})

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + vec[4]*rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K)

  expect_true(length(unique(res$cluster)) == 4)

})

##########################

## .adjustment_cluster is correct

test_that(".adjustment_cluster correctly removes one cluster", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster(clust, 5)

  expect_true(all(sort(unique(res)) == c(1:5)))
})

test_that(".adjustment_cluster correctly adds one cluster", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster(clust, 7)

  expect_true(all(sort(unique(res)) == c(1:7)))
})

test_that(".adjustment_cluster correctly removes many clusters", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster(clust, 4)

  expect_true(all(sort(unique(res)) == c(1:4)))
})

test_that(".adjustment_cluster correctly adds many clusters", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster(clust, 9)

  expect_true(all(sort(unique(res)) == c(1:9)))
})

############################

## .adjustment_cluster_add is correct

test_that(".adjustment_cluster_add adds", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster_add(clust)

  expect_true(all(sort(unique(res)) == c(1:7)))
})
