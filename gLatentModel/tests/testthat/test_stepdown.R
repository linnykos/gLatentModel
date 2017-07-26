context("Test stepdown")

## stepdown is correct

test_that("stepdown works", {
  set.seed(10)
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

  dat <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = cov_mat)

  g_list <- list(row_difference_closure(1,4,4), row_difference_closure(1,2,4))
  res <- stepdown(dat, g_list)

  expect_true(is.numeric(res))
  expect_true(all(res %% 1 == 0))
})
