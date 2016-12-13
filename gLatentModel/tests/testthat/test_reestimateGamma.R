context("Test reestimate gamma")

## reestimate_gamma is correct

test_that("reestimate_gamma returns properly", {
  set.seed(10)
  a_mat <- matrix(0, 20, 10)
  a_mat[1:10, 1:10] <- diag(10)
  for(i in 1:50){
    i <- sample(11:20, 1); j <- sample(1:10, 1)
    a_mat[i,j] <- 1
  }
  a_mat <- t(apply(a_mat, 1, function(x){
    if(all(x == 0)) return(rep(1/10, 10)) else return(x/sum(x))}))
  group_list <- group_cluster(a_mat)
  dat <- MASS::mvrnorm(20, rep(0, 10), diag(10))

  res <- reestimate_gamma(dat, group_list)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
  expect_true(all(res > 0))
})
