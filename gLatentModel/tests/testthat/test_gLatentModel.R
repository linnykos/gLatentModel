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

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + rnorm(prod(dim(dat)))

  idx <- sample(1:ncol(dat))
  a_mat2 <- a_mat[idx,]
  dat2 <- dat[,idx]

  res <- gLatentModel(dat, K, seed = 10)
  res <- reshuffle(res, a_mat)

  res2 <- gLatentModel(dat2, K, seed = 10)
  res2 <- reshuffle(res2, a_mat2)

  expect_true(sum(abs(res$cov_latent - res2$cov_latent))/(K^2) < 1e-2)
})
