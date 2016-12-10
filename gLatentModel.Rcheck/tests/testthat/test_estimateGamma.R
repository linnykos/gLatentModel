context("Test estimate gamma")

## estimate_gamma is correct

test_that("estimate_gamma returns a vector", {
  mat <- matrix(1:25, 5, 5)
  res <- estimate_gamma(mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == 5)
})

test_that("estimate_gamma is roughly correct", {
  set.seed(10)
  mat <- matrix(rnorm(1000), 100)
  res <- estimate_gamma(mat)

  mat2 <- matrix(rnorm(100), 10)
  res2 <- estimate_gamma(mat2)

  expect_true(sum((res - 1)^2) <= sum((res2 - 1)^2))
})

################################

## .l2norm is correct

test_that(".l2norm is correct", {
  vec <- c(1:5)
  res <- .l2norm(vec)

  expect_true(res == sqrt(1+4+9+16+25))
})

#############################

## .v_innerproduct is correct

test_that(".v_innerproduct can return a number", {
  mat <- matrix(1:25, 5, 5)
  res <- .v_innerproduct(mat, 1, 2)

  expect_true(is.numeric(res))
})

test_that(".v_innerproduct returns 0 properly", {
  mat <- matrix(0, 5, 5)
  mat[,1:2] <- 1
  res <- .v_innerproduct(mat, 1, 2)

  expect_true(res == 0)
})

test_that(".v_innerproduct can return values in right relation", {
  set.seed(10)
  mat <- matrix(rnorm(100), 10, 10)
  mat[,1] <- c(-5:4)
  for(i in 2:4){
    mat[,i] <- mat[,1] + 0.01*rnorm(10)
  }

  res <- .v_innerproduct(mat, 1, 2)
  res2 <- .v_innerproduct(mat, 1, 5)

  expect_true(res <= res2)
})

#####################################

## .neighbor_gamma_estimation is correct

test_that(".neighbor_gamma_estimation returns an index", {
  mat <- matrix(1:25, 5, 5)
  res <- .neighbor_gamma_estimation(mat, 1)

  expect_true(is.numeric(res))
  expect_true(res %% 1 == 0)
  expect_true(res >= 2)
  expect_true(res <= 5)
})

test_that(".neighbor_gamma_estimation can handle a vector for a.vec", {
  mat <- matrix(1:25, 5, 5)
  res <- .neighbor_gamma_estimation(mat, c(1,4))

  expect_true(is.numeric(res))
  expect_true(!res %in% c(1,4))
})

test_that(".neighbor_gamma_estimation returns the right relative order", {
 set.seed(10)
  mat <- matrix(rnorm(100), 10, 10)
  mat[,2] <- mat[,1] + 0.01*rnorm(10)

  res <- .neighbor_gamma_estimation(mat, 1)

  expect_true(res == 2)
})
