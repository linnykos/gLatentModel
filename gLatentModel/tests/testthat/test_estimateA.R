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


#####################################

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

#######################################

## .KKT_grid_solver is correct

test_that(".KKT_grid_solver returns properly", {
  set.seed(10)
  gamma <- 1
  K <- 5
  vec <- stats::rnorm(K)
  mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)

  inv_mat <- solve(t(mat) %*% mat)

  res <- .KKT_grid_solver(gamma, vec, mat, inv_mat)

  expect_true(length(res) == K)
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
})

# test_that(".KKT_grid_solver enforces nonnegativety", {
#   for(trial in 1:50){
#     set.seed(10*trial)
#     gamma <- 1
#     K <- 5
#     vec <- stats::rnorm(K)
#     mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)
#
#     inv_mat <- solve(t(mat) %*% mat)
#
#     res <- .KKT_grid_solver(gamma, vec, mat, inv_mat)
#
#     expect_true(all(res >= -1e-4)) #FOR SOME REASON, NOT ALWAYS NON-NEGATIVE
#   }
# })
#
# test_that(".KKT_grid_solver enforces complimentary slackness",{
#   for(trial in 1:50){
#     set.seed(10*trial)
#     gamma <- 1
#     K <- 5
#     vec <- stats::rnorm(K)
#     mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)
#
#     inv_mat <- solve(t(mat) %*% mat)
#
#     res <- .KKT_grid_solver(gamma, vec, mat, inv_mat)
#
#     lambda <- -2*(vec - res %*% mat) %*% mat + gamma
#     expect_true(all(lambda * res >= -1e-4))
#   }
# })
#
# test_that(".KKT_grid_solver enforces dual feasibility",{
#   for(trial in 1:50){
#     set.seed(10*trial)
#     gamma <- 1
#     K <- 5
#     vec <- stats::rnorm(K)
#     mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)
#
#     inv_mat <- solve(t(mat) %*% mat)
#
#     res <- .KKT_grid_solver(gamma, vec, mat, inv_mat)
#
#     lambda <- -2*(vec - res %*% mat) %*% mat + gamma
#     expect_true(all(lambda >= -1e-4))
#   }
# })

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

# test_that(".optim_solver_constrainLS enforces primal feasibility",{
#   for(trial in 1:50){
#     set.seed(10*trial)
#     K <- 5
#     vec <- stats::rnorm(K)
#     mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)
#
#     res <- .optim_solver_constrainLS(vec, mat)
#
#     expect_true(all(res >= -1e-4))
#     expect_true(abs(sum(res) - 1) <= 1e-4)
#   }
# })
#
# test_that(".optim_solver_contrainLS gives the same value of gamma", {
#   for(trial in 1:50){
#     set.seed(10*trial)
#     K <- 5
#     vec <- stats::rnorm(K)
#     mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)
#
#     res <- .optim_solver_constrainLS(vec, mat)
#
#     zz = -2*(vec-res%*%mat)%*%mat
#     idx = which(res != 0)
#     gamma = mean(-zz[idx])
#
#     expect_true(all(abs(-zz[idx] - gamma) <= 1e-4))
#   }
# })
#
# test_that(".optim_solver_constrainLS enforces dual feasibility", {
#   for(trial in 1:50){
#     set.seed(10*trial)
#     K <- 5
#     vec <- stats::rnorm(K)
#     mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)
#
#     res <- .optim_solver_constrainLS(vec, mat)
#
#     zz = -2*(vec-res%*%mat)%*%mat
#     idx = which(res != 0)
#     gamma = mean(-zz[idx])
#     lambda = rep(0, K)
#     lambda[-idx] = zz[-idx] + gamma
#
#     expect_true(all(lambda >= -1e-4))
#   }
# })
#
# test_that(".optim_solver_constrainLS enforces complimentary slackness", {
#   for(trial in 1:50){
#     set.seed(10*trial)
#     K <- 5
#     vec <- stats::rnorm(K)
#     mat <- matrix(stats::rnorm(K^2), K, K); mat <- mat + t(mat)
#
#     res <- .optim_solver_constrainLS(vec, mat)
#
#     zz = -2*(vec-res%*%mat)%*%mat
#     idx = which(res != 0)
#     gamma = mean(-zz[idx])
#     lambda = rep(0, K)
#     lambda[-idx] = zz[-idx] + gamma
#
#     expect_true(all(lambda*res >= -1e-4))
#   }
# })
