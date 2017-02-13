context("Test adjustment cluster")

## .adjustment_cluster is correct

test_that(".adjustment_cluster correctly removes one cluster", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster(clust, 5)

  expect_true(all(sort(unique(res)) == c(1:5)))
})

test_that(".adjustment_cluster correctly adds one cluster", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster(clust, 7)

  expect_true(all(sort(unique(res)) == c(1:7)))
})

test_that(".adjustment_cluster correctly removes many clusters", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster(clust, 4)

  expect_true(all(sort(unique(res)) == c(1:4)))
})

test_that(".adjustment_cluster correctly adds many clusters", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster(clust, 9)

  expect_true(all(sort(unique(res)) == c(1:9)))
})

############################

## .adjustment_cluster_add is correct

test_that(".adjustment_cluster_add adds", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster_add(clust)

  expect_true(all(sort(unique(res)) == c(1:7)))
})

test_that(".adjustment_cluster_add adds when there's data", {
  clust <- rep(1:6, each = 5)
  dat <- matrix(1:300, 10, 30)
  res <- .adjustment_cluster_add(clust, dat)

  expect_true(all(sort(unique(res)) == c(1:7)))
})

##############################

## .adjustment_cluster_remove is correct

test_that(".adjustment_cluster_remove adds", {
  clust <- rep(1:6, each = 5)
  res <- .adjustment_cluster_remove(clust)

  expect_true(all(sort(unique(res)) == c(1:5)))
})

test_that(".adjustment_cluster_remove adds when there's data", {
  clust <- rep(1:6, each = 5)
  dat <- matrix(1:300, 10, 30)
  res <- .adjustment_cluster_remove(clust, dat)
})


