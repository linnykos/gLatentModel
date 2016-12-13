context("Test reestimate gamma")

## reestimate_gamma is correct

test_that("reestimate_gamma returns properly", {
  set.seed(10)
  a_mat <- rbind(diag(10), diag(10))
  partition_list <- partition_cluster(a_mat)
  dat <- MASS::mvrnorm(20, rep(0, 10), diag(10))

  res <- reestimate_gamma(dat, partition_list)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 20)
  expect_true(all(res > 0))
  expect_true(length(unique(res)) == 10)
})
