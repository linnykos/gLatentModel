context("Test clustering representations")

## .convert_list_to_vector is correct

test_that(".convert_list_to_vector works", {
  lis <- lapply(1:4, function(x){
    seq(x, 40, by = 4)
  })

  res <- .convert_list_to_vector(lis)

  expect_true(length(res) == 40)
  expect_true(all(res %in% c(1:4)))
  expect_true(all(res == rep(1:4, times = 10)))
})

test_that(".convert_list_to_vector works when the list is of length 1", {
  lis <- list(c(1:10))
  res <- .convert_list_to_vector(lis)
  expect_true(all(res == rep(1,10)))
})

## .convert_vector_to_list is correct

test_that(".convert_vector_to_list works", {
  vec <- rep(1:4, times = 10)
  res <- .convert_vector_to_list(vec)

  lis <- lapply(1:4, function(x){
    seq(x, 40, by = 4)
  })

  expect_true(length(res) == 4)
  for(i in 1:4){
    expect_true(all(res[[i]] == lis[[i]]))
  }
})
