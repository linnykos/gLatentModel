context("Test reshuffle")

# .reshuffle is correct

test_that(".reshuffle works", {
  set.seed(10)
  K <- 4; n <- 100; times <- 3
  L <- huge::huge.generator(n = n, d = K, graph = "hub", g = 3, verbose = F)
  latent_dat <- L$data

  a_mat <- diag(K)
  for(i in 1:(times-1)){a_mat <- rbind(a_mat, diag(K))}

  dat <- latent_dat%*%t(a_mat)
  dat <- dat + rnorm(prod(dim(dat)))

  res <- gLatentModel(dat, K, seed = 10)

  cluster_vec2 <- c(4,2,1,3)[res$cluster]
  res2 <- res
  res2$cluster <- cluster_vec2; res2$cov_latent <- res2$cov_latent[c(4,2,1,3), c(4,2,1,3)]

  res3 <- .reshuffle(res2, res$cluster)

  expect_true(all(res$cluster == res3$cluster))
  expect_true(any(res$cluster != res2$cluster))

  expect_true(sum(abs(res$cov_latent - res3$cov_latent))/K^2 < 1e-4)
  expect_true(sum(abs(res$cov_latent - res2$cov_latent))/K^2 > 1e-4)
})
