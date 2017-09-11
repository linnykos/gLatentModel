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
