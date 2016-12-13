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
  expect_true(all(dim(res$theta) == c(2,2)))
  expect_true(is.numeric(res$theta))
  expect_true(is.matrix(res$theta))
})
