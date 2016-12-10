context("Test estimate gamma")

test_that(".l2norm is correct", {
  vec <- c(1:5)
  res <- .l2norm(vec)

  expect_true(res == sqrt(1+4+9+16+25))
})
