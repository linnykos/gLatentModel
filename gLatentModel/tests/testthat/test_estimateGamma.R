context("Test estimate gamma")

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
