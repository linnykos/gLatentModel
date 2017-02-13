context("Test estimate gamma")

## .estimate_cov_noise is correct

test_that(".estimate_cov_noise returns a vector", {
  mat <- matrix(1:25, 5, 5)
  res <- .estimate_cov_noise(mat)

  expect_true(is.numeric(res))
  expect_true(length(res) == 5)
})

test_that(".estimate_cov_noise is roughly correct", {
  set.seed(10)
  mat <- matrix(rnorm(1000), 100)
  res <- .estimate_cov_noise(mat)

  mat2 <- matrix(rnorm(100), 10)
  res2 <- .estimate_cov_noise(mat2)

  expect_true(sum((res - 1)^2) <= sum((res2 - 1)^2))
})

test_that(".estimate_cov_noise returns 0's when the matrix has duplicate cols", {
  vec <- c(4, 6, 100)
  set.seed(2)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- diag(K)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L)

  a_mat <- sapply(1:K, function(x){
    vec <- rep(0, K*times)
    vec[((x-1)*times+1):(x*times)] <- 1
    vec
  })

  dat <- latent_dat%*%t(a_mat)
  res <- .estimate_cov_noise(dat)

  expect_true(all(res == 0))
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

## .neighbor_noise_estimation is correct

test_that(".neighbor_noise_estimation returns an index", {
  mat <- matrix(1:25, 5, 5)
  res <- .neighbor_noise_estimation(mat, 1)

  expect_true(is.numeric(res))
  expect_true(res %% 1 == 0)
  expect_true(res >= 2)
  expect_true(res <= 5)
})

test_that(".neighbor_noise_estimation can handle a vector for a.vec", {
  mat <- matrix(1:25, 5, 5)
  res <- .neighbor_noise_estimation(mat, c(1,4))

  expect_true(is.numeric(res))
  expect_true(!res %in% c(1,4))
})

test_that(".neighbor_noise_estimation returns the right relative order", {
  set.seed(10)
  mat <- matrix(rnorm(100), 10, 10)
  mat[,2] <- mat[,1] + 0.01*rnorm(10)

  res <- .neighbor_noise_estimation(mat, 1)

  expect_true(res == 2)
})

test_that(".neighbor_noise_estimation can handle with duplicate columns", {
  vec <- c(4, 6, 100)
  set.seed(2)

  K <- vec[1]; n <- vec[3]; times <- vec[2]
  L <- diag(K)
  latent_dat <- MASS::mvrnorm(n, rep(0,K), L)

  a_mat <- sapply(1:K, function(x){
    vec <- rep(0, K*times)
    vec[((x-1)*times+1):(x*times)] <- 1
    vec
  })

  dat <- latent_dat%*%t(a_mat)

  x <- 19
  ne1 <- .neighbor_noise_estimation(dat, x, NA)
  ne2 <- .neighbor_noise_estimation(dat, c(x, ne1), NA)

  expect_true(all(c(ne1, ne2) >= 18))
})
