context("Test distance between partitions")

## .check_partition is correct

test_that(".check_partition works", {
  vec <- rep(1:4, times = 10)
  bool <- .check_partition(vec)

  expect_true(bool)
})

test_that(".check_partition fails when labels are missing", {
  vec <- rep(1:4, times = 10)
  vec[vec == 2] <- 1
  expect_error(.check_partition(vec))
})

test_that(".check_partition fails for non-vectors", {
  vec <- rep(1:4, times = 10)
  tmp <- as.matrix(vec)
  expect_error(.check_partition(tmp))

  tmp <- as.factor(vec)
  expect_error(.check_partition(tmp))

  tmp <- as.list(vec)
  expect_error(.check_partition(tmp))
})

#########

## .partition_variation is correct

test_that(".partition_variation works", {
  set.seed(10)
  partition1 <- sample(1:4, 40, replace = T)
  partition2 <- sample(1:4, 40, replace = T)

  res <- .partition_variation(partition1, partition2)
  expect_true(is.numeric(res))
  expect_true(res > 0 )
})

test_that(".partition_variation returns 0 when the partition is the same", {
  set.seed(10)
  partition1 <- sample(1:4, 40, replace = T)
  partition2 <- partition1
  for(i in 1:4){partition2[partition2 == i] <- -i}
  partition2 <- partition2 + 5

  res <- .partition_variation(partition1, partition2)
  expect_true(res == 0)
})

test_that(".partition_variation is smaller for smaller variation", {
  #this is a heuristic check I believe should work
  set.seed(10)
  partition1 <- sample(1:4, 40, replace = T)
  partition2 <- sample(1:4, 40, replace = T)
  partition3 <- partition1

  partition3[1] <- (partition3[1]+1)%%4
  if(partition3[1] == 0) partition3[1] <- 4

  res1 <- .partition_variation(partition1, partition2)
  res2 <- .partition_variation(partition1, partition3)

  expect_true(res1 >= res2)
})

##########

## .partition_hamming is correct

test_that(".partition_hamming works", {
  set.seed(10)
  partition1 <- sample(1:4, 40, replace = T)
  partition2 <- sample(1:4, 40, replace = T)

  res <- .partition_hamming(partition1, partition2)
  expect_true(is.numeric(res))
  expect_true(res > 0 )
})

test_that(".partition_hamming returns 0 when the partition is the same", {
  set.seed(10)
  partition1 <- sample(1:4, 40, replace = T)
  partition2 <- partition1
  for(i in 1:4){partition2[partition2 == i] <- -i}
  partition2 <- partition2 + 5

  res <- .partition_hamming(partition1, partition2)
  expect_true(res == 0)
})

test_that(".partition_hamming is smaller for smaller variation", {
  #this is a heuristic check I believe should work
  set.seed(10)
  partition1 <- sample(1:4, 40, replace = T)
  partition2 <- sample(1:4, 40, replace = T)
  partition3 <- partition1

  partition3[1] <- (partition3[1]+1)%%4
  if(partition3[1] == 0) partition3[1] <- 4

  res1 <- .partition_hamming(partition1, partition2)
  res2 <- .partition_hamming(partition1, partition3)

  expect_true(res1 >= res2)
})

#########

## .hamming_check is correct

test_that(".hamming_check works", {
  set.seed(10)
  pair <- c(1,10)
  partition1 <- sample(1:4, 40, replace = T)
  partition2 <- sample(1:4, 40, replace = T)

  res <- .hamming_check(pair, partition1, partition2)

  expect_true(is.logical(res))
})

test_that(".hamming_check returns TRUE if pair are same partitions", {
  pair <- c(1,10)
  partition1 <- rep(1:4, times = 10)
  partition2 <- rep(4:1, times = 10)
  res <- .hamming_check(pair, partition1, partition2)

  expect_true(res)
})

test_that(".hamming_check returns FALSE if pair are different partitions", {
  pair <- c(1,10)
  partition1 <- rep(1:4, times = 10)
  partition2 <- rep(4:1, times = 10)

  partition2[10] <- (partition2[10]+1)%%4
  if(partition2[10] == 0) partition2[10] <- 4

  res <- .hamming_check(pair, partition1, partition2)

  expect_true(!res)
})

##########

## .l2vec is correct

test_that(".l2vec works", {
  res <- .l2vec(1:5)
  expect_true(is.numeric(res))
})

############

## .partition_hamming_covariance is correct

test_that(".partition_hamming_covariance works", {
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

  partition1 <- c(1,1,2,2)
  partition2 <- c(1,2,2,1)

  res <- .partition_hamming_covariance(partition1, partition2, cov_mat)

  expect_true(is.numeric(res))
})

test_that(".partition_hamming_covariance returns 0 if partitions are the same", {
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

  partition1 <- rep(1:4, times = 10)
  partition2 <- rep(4:1, times = 10)

  res <- .partition_hamming_covariance(partition1, partition2, cov_mat)

  expect_true(res == 0)
})

test_that(".partition_hamming_covariance scales correctly", {
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

  cov_mat2 <- cov_mat*5

  partition1 <- c(1,1,2,2)
  partition2 <- c(1,2,2,1)

  res1 <- .partition_hamming_covariance(partition1, partition2, cov_mat)
  res2 <- .partition_hamming_covariance(partition1, partition2, cov_mat2)

  expect_true(res2 > res1)
})

test_that(".partition_hamming_covariance can handle uneven partitions", {
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

  partition1 <- rep(1,4)
  partition2 <- c(1,2,2,1)

  res <- .partition_hamming_covariance(partition1, partition2, cov_mat)

  expect_true(is.numeric(res))
})
