context("Test reestimate theta")

## .reestimate_theta is correct

test_that(".reestimate_theta returns properly", {
  set.seed(10)
  dat <- MASS::mvrnorm(20, rep(0, 6), diag(6))
  cov_mat <- stats::cov(dat)
  gamma_mat <- diag(stats::rnorm(6))
  a_mat <- rbind(diag(3), diag(3))
  partition_list = .partition_cluster(a_mat)

  res <- .reestimate_theta(cov_mat, gamma_mat, partition_list)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(sum(abs(res - t(res))) < 1e-4)
  expect_true(all(dim(res) == c(3,3)))
})
