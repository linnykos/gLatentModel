context("Test reshuffle of gModel")

## .reshuffle is correct

test_that(".reshuffle relabels the clusters correctly", {
  set.seed(10)
  K <- 4; n <- 100; times <- 6
  L <- huge::huge.generator(n = n, d = K, graph = "hub", g = 3, verbose = F)
  latent_dat <- L$data

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K)
  res2 <- .reshuffle(res, a_mat)

  expect_true(all(sort(res2$partition_list[[1]]) == c(1,5,9,13,17,21)))
  expect_true(all(sort(res2$partition_list[[2]]) == c(2,6,10,14,18,22)))
  expect_true(all(sort(res2$partition_list[[3]]) == c(3,7,11,15,19,23)))
  expect_true(all(sort(res2$partition_list[[4]]) == c(4,8,12,16,20,24)))
})
