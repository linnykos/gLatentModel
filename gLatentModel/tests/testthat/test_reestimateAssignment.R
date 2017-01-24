context("Test reestimate assignment")

# .reestimate_assignment is correct

test_that(".reestimate_assignment works", {
  cluster_vec <- c(1,4,2,3,1,4,1,3)
  res <- .reestimate_assignment(cluster_vec)

  true_res <- matrix(c(1,0,0,0,1,0,1,0,
                       0,0,1,0,0,0,0,0,
                       0,0,0,1,0,0,0,1,
                       0,1,0,0,0,1,0,0), 8, 4)
  expect_true(all(true_res == res))
})
