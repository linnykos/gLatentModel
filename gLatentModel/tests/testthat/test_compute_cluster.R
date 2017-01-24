context("compute_cluster.R")

###################################################

## .estimate_partition_mat works properly

test_that(".estimate_partition_mat returns properly", {
  set.seed(10)
  a_mat <- matrix(stats::rnorm(50), 10, 5)

  cov_latent <- matrix(stats::rnorm(25), 5, 5)
  cov_latent <- (cov_latent + t(cov_latent))/2
  eig <- eigen(cov_latent); eig$values[eig$values < 0] <- 0
  cov_latent <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

  res <- .estimate_partition_mat(a_mat, cov_latent)
  res_eig <- eigen(res)

  expect_true(all(dim(res) == c(10,10)))
  expect_true(all(res_eig$values >= -1e-4))
  expect_true(sum(abs(res - t(res)))/100 <= 1e-4)
})

##################################################

## .estimate_cluster works properly

test_that(".estimate_cluster returns properly", {
  set.seed(10)
  partition_mat <- matrix(0, 10, 10)
  partition_mat[1:5, 1:5] <- 5 + stats::rnorm(25)
  partition_mat[6:10, 6:10] <- 5 + stats::rnorm(25)

  res <- .estimate_cluster(partition_mat, 2)

  expect_true(length(res) == 10)
  expect_true(length(unique(res[1:5])) == 1)
  expect_true(length(unique(res[6:10])) == 1)
  expect_true(length(unique(res)) == 2)
})
