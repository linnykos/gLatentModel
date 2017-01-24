context("Test reestimate cov denoise")

## .reestimate_gamma is correct

test_that(".reestimate_gamma returns properly", {
  set.seed(10)
  cluster_vec <- c(1:10, 1:10)
  dat <- MASS::mvrnorm(50, rep(0, 20), diag(20))

  res <- .reestimate_gamma(dat, cluster_vec)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 20)
  expect_true(all(res > 0))
  expect_true(length(unique(res)) == 10)
})
