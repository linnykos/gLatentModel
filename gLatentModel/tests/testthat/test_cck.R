context("Test cck")

## row_difference_closure is correct

test_that("row_difference_closure works", {
  d <- 10
  dat <- matrix(0, d, d)
  dat[,1] <- 5; dat[1,] <- 5
  dat[,2] <- 1; dat[2,] <- 1
  vec <- dat[lower.tri(dat)]

  g <- row_difference_closure(1,2,d)
  res <- g(vec)

  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(res == 4)
})

test_that("row_difference_closure gives the right value", {
  d <- 10
  dat <- matrix(0, d, d)
  dat[,5] <- 1:d; dat[1,] <- 1:d
  dat[,2] <- 2*(1:d); dat[2,] <- 2*(1:d)
  vec <- dat[lower.tri(dat)]

  g <- row_difference_closure(2,5,d)
  res <- g(vec)

  expect_true(res == 10)
})

########################

## cor_vec is correct

test_that("cor_vec works", {
  dat <- matrix(rnorm(40), 8, 5)
  res <- cor_vec(dat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 5*4/2)
})

##########################

## cck is correct

test_that("cck works", {
  set.seed(10)
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

  dat <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = cov_mat)
  g <- row_difference_closure(1,4,4)
  res <- cck(dat, g = g)

  expect_true(length(res) == 3)
  expect_true(is.list(res))
  expect_true(res$pval <= 1)
  expect_true(res$pval >= 0)
})

test_that("cck has sensible p-values", {
  set.seed(10)
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

  dat <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = cov_mat)
  set.seed(10)
  g1 <- row_difference_closure(1,4,4)
  res1 <- cck(dat, g = g1)

  set.seed(10)
  g2 <- row_difference_closure(1,2,4)
  res2 <- cck(dat, g = g2)

  expect_true(res2$pval < res1$pval)
  expect_true(res1$t0 <= res2$t0)
})
