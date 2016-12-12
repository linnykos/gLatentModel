context("Test estimating A")

## pure_node is correct

test_that("pure_node returns properly", {
  set.seed(10)
  K <- 5
  dat <- MASS::mvrnorm(100, rep(0, 10), diag(c(rep(5,K), rep(0.1, 10-K))))
  c_mat <- cov(dat)

  res <- pure_nodes(c_mat, K)

  expect_true(is.numeric(res))
  expect_true(length(res) == K)
  expect_true(all(res%%1 == 0))
  expect_true(all(sort(res) == c(1:5)))
})
