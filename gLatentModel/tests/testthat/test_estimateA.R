context("Test estimating A")

## estimate_theta is correct

test_that("estimate_theta returns properly", {
  set.seed(10)
  K <- 5
  dat <- MASS::mvrnorm(100, rep(0, 10), diag(c(rep(5,K), rep(0.1, 10-K))))
  c_mat <- cov(dat)
  idx <- c(1,5,10)

  res <- estimate_theta(c_mat, idx)
  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == rep(length(idx), 2)))
  expect_true(sum(abs(res - t(res))) < 1e-4)
  expect_true(sum(abs(res - c_mat[idx, idx])) < 1e-4)
})

##########################################

## .optim_solver_constrainLS is correct

test_that(".optim_solver_constrainLS returns properly", {
  set.seed(10)
  K <- 5
  vec <- stats::rnorm(K)
  mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)

  res <- .optim_solver_constrainLS(vec, mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == K)
 })


test_that(".optim_solver_constrainLS gets an obvious answer correct", {
  set.seed(10)
  vec <- c(1:5)
  mat <- matrix(rnorm(25),5,5); mat <- mat + t(mat)
  mat[,1] <- vec; mat[1,] <- vec

  res <- .optim_solver_constrainLS(vec, mat)

  expect_true(sum((vec - mat %*% res)^2) - 1e-4 <= sum((vec - mat %*% c(1,rep(0,4)))^2))
})

test_that(".optim_solver_constrainLS enforces primal feasibility",{
  for(trial in 1:10){
    set.seed(10*trial)
    K <- 5
    vec <- stats::rnorm(K)
    mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)

    res <- .optim_solver_constrainLS(vec, mat)

    expect_true(all(res >= -1e-4))
    expect_true(abs(sum(res) - 1) <= 1e-4)
  }
})

test_that(".optim_solver_contrainLS gives the same value of gamma", {
  for(trial in 1:10){
    set.seed(10*trial)
    K <- 5
    vec <- stats::rnorm(K)
    mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)

    res <- .optim_solver_constrainLS(vec, mat, rho = 10)

    zz = -2*(vec-res%*%mat)%*%mat
    idx = which(res != 0)
    gamma = mean(-zz[idx])

    expect_true(all(abs(-zz[idx] - gamma) <= 1e-2))
  }
})

test_that(".optim_solver_constrainLS enforces dual feasibility", {
  for(trial in 1:10){
    set.seed(10*trial)
    K <- 5
    vec <- stats::rnorm(K)
    mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)

    res <- .optim_solver_constrainLS(vec, mat, rho = 10)

    zz = -2*t(mat)%*%(vec-mat%*%res)
    idx = which(res != 0)
    gamma = mean(zz[idx])
    lambda = rep(0, K)
    lambda[-idx] = zz[-idx] - gamma

    expect_true(all(lambda >= -1e-2))
  }
})

test_that(".optim_solver_constrainLS enforces complimentary slackness", {
  for(trial in 1:10){
    set.seed(10*trial)
    K <- 5
    vec <- stats::rnorm(K)
    mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)

    res <- .optim_solver_constrainLS(vec, mat)

    zz = -2*(vec-res%*%mat)%*%mat
    idx = which(res != 0)
    gamma = mean(-zz[idx])
    lambda = rep(0, K)
    lambda[-idx] = zz[-idx] + gamma

    expect_true(all(lambda*res >= -1e-2))
  }
})

####################################################

## estimate_a works properly

test_that("estimate_a returns properly", {
  set.seed(10)
  theta_mat <- matrix(rnorm(25), 5, 5); theta_mat <- theta_mat + t(theta_mat)
  c_mat <- matrix(rnorm(100), 10, 10); c_mat <- c_mat + t(c_mat)
  pure_idx <- c(1,3,5,7,10)

  res <- estimate_a(theta_mat, c_mat, pure_idx)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(10, 5)))
  expect_true(all(res[1:5, 1:5] == diag(5)))
})

###################################################

## group_cluster works properly

test_that("group_cluster returns properly", {
  set.seed(10)
  a_mat <- matrix(0, 20, 10)
  a_mat[1:10, 1:10] <- diag(10)
  for(i in 1:50){
    i <- sample(11:20, 1); j <- sample(1:10, 1)
    a_mat[i,j] <- 1
  }
  a_mat <- t(apply(a_mat, 1, function(x){
    if(all(x == 0)) return(rep(1/10, 10)) else return(x/sum(x))}))

  res <- group_cluster(a_mat)

  expect_true(is.list(res))
  expect_true(length(res) == 10)
  expect_true(all(sapply(res, length)) >= 1)
  expect_true(all(sapply(res, function(x){all(x %% 1 == 0)})))
})

##################################################

## partition_cluster works properly

test_that("partition_cluster returns properly", {
  set.seed(10)
  a_mat <- matrix(0, 20, 10)
  a_mat[1:10, 1:10] <- diag(10)
  for(i in 1:50){
    i <- sample(11:20, 1); j <- sample(1:10, 1)
    a_mat[i,j] <- 1
  }
  a_mat <- t(apply(a_mat, 1, function(x){
    if(all(x == 0)) return(rep(1/10, 10)) else return(x/sum(x))}))

  res <- partition_cluster(a_mat)

  expect_true(is.list(res))
  expect_true(length(res) == 10)
  expect_true(all(sapply(res, length)) >= 1)
  expect_true(all(sapply(res, function(x){all(x %% 1 == 0)})))

  expect_true(all(sort(unlist(res)) == 1:20))
})
