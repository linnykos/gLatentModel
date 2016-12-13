context("Test reshuffle of gModel")

## reshuffle is correct

test_that("reshuffle relabels the clusters correctly", {
  set.seed(20)
  K <- 4; n <- 100
  L <- huge::huge.generator(n = n, d = K, graph = "hub", g = 3)
  latent_dat <- L$data

  a_mat <- rbind(diag(K), diag(K), diag(K), diag(K), diag(K), diag(K))
  dat <- latent_dat%*%t(a_mat)
  dat <- dat + 0.01*stats::rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K)
  res2 <- reshuffle(res, a_mat)

  expect_true(all(sort(res2$partition_list[[1]]) == c(1,5,9,13,17,21)))
  expect_true(all(sort(res2$partition_list[[2]]) == c(2,6,10,14,18,22)))
  expect_true(all(sort(res2$partition_list[[3]]) == c(3,7,11,15,19,23)))
  expect_true(all(sort(res2$partition_list[[4]]) == c(4,8,12,16,20,24)))
})
