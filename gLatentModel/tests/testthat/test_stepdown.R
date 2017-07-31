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

test_that("stepdown finds the right partition", {
  set.seed(10)
  d <- 6
  cov_mat <- matrix(0, d, d)
  cov_mat[1:(d/2), 1:(d/2)] <- 0.75
  cov_mat[(d/2+1):d, (d/2+1):d] <- 0.75
  diag(cov_mat) <- 1

  dat <- MASS::mvrnorm(100, mu = rep(0,d), Sigma = cov_mat)

  combn_mat <- combn(d, 2)
  g_list <- lapply(1:ncol(combn_mat), function(x){
    row_difference_closure(combn_mat[1,x], combn_mat[2,x], d)})

  res <- stepdown(dat, g_list)
  correct_idx <- which(apply(combn_mat, 2, function(x){
    bool1 <- x[1] <= d/2
    bool2 <- x[2] <= d/2
    if(bool1 == bool2) TRUE else FALSE
  }))

  expect_true(all(sort(res) == sort(correct_idx)))
})

############

## connected_components is correct

test_that("connected_components works", {
  set.seed(10)
  d <- 6
  cov_mat <- matrix(0, d, d)
  cov_mat[1:(d/2), 1:(d/2)] <- 0.75
  cov_mat[(d/2+1):d, (d/2+1):d] <- 0.75
  diag(cov_mat) <- 1

  dat <- MASS::mvrnorm(100, mu = rep(0,d), Sigma = cov_mat)

  combn_mat <- combn(d, 2)
  g_list <- lapply(1:ncol(combn_mat), function(x){
    row_difference_closure(combn_mat[1,x], combn_mat[2,x], d)})

  idx <- stepdown(dat, g_list)
  res <- connected_components(d, combn_mat[,idx])

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(res[[1]] == c(1,2,3)))
  expect_true(all(res[[2]] == c(4,5,6)))
})
