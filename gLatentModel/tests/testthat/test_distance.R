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
