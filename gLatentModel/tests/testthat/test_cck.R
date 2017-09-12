context("Test cck")

## row_difference_closure is correct

test_that("row_difference_closure works", {
  d <- 10
  dat <- matrix(0, d, d)
  dat[,1] <- 5; dat[1,] <- 5
  dat[,2] <- 1; dat[2,] <- 1
  vec <- dat[lower.tri(dat)]

  g <- row_difference_closure(1,2,d)
  res <- g(vec)

  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(res == 4)
})

test_that("row_difference_closure gives the right value", {
  d <- 10
  dat <- matrix(0, d, d)
  dat[,5] <- 1:d; dat[1,] <- 1:d
  dat[,2] <- 2*(1:d); dat[2,] <- 2*(1:d)
  vec <- dat[lower.tri(dat)]

  g <- row_difference_closure(2,5,d)
  res <- g(vec)

  expect_true(res == 10)
})

test_that("row_difference_closure gives a function that does the right calculation",{
  set.seed(10)
  n <- 100; d <- 4
  cov_mat <- diag(d)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5
  dat <- MASS::mvrnorm(n, mu = rep(0,4), Sigma = cov_mat)

  e <- stats::rnorm(n)
  sigma_vec <- apply(dat, 2, stats::sd)
  psi <- cor_vec(dat, sigma_vec = sigma_vec)
  psi_boot <- cor_vec(dat, sigma_vec = sigma_vec, noise_vec = e)
  g <- row_difference_closure(1,2,d)
  res <- g(psi_boot, average_vec = psi*sum(e)/n)

  #manual calculation
  dat <- scale(dat, scale = F)
  combn_mat <- combn(d, 2)
  psi_boot2 <- sapply(1:ncol(combn_mat), function(x){
    denom <- sigma_vec[combn_mat[1,x]]*sigma_vec[combn_mat[2,x]]
    val <- sapply(1:n, function(y){
      (dat[y,combn_mat[1,x]]*dat[y,combn_mat[2,x]]/denom - psi[x])*e[y]
    })
    mean(val)
  })
  #hard-coded indices
  diff_vec <- c(psi_boot2[2]-psi_boot2[4], psi_boot2[3]-psi_boot2[5])
  res2 <- max(abs(diff_vec))

  expect_true(abs(res - res2) <= 1e-5)
})

########################

## cor_vec is correct

test_that("cor_vec works", {
  dat <- matrix(rnorm(40), 8, 5)
  res <- cor_vec(dat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 5*4/2)
})

test_that("cor_vec gives the correlation matrix when noise_vec is 1", {
  set.seed(10)
  dat <- matrix(rnorm(40), 8, 5)
  sigma_vec <- apply(dat, 2, stats::sd)

  res <- cor_vec(dat, sigma_vec = sigma_vec)
  res2 <- stats::cor(dat)
  res2 <- res2[lower.tri(res2)]*(7/8)
})

test_that("cor_vec does the calculation correctly with noise_vec", {
  set.seed(10)
  dat <- matrix(rnorm(40), 8, 5)
  sigma_vec <- apply(dat, 2, stats::sd)
  e <- rnorm(8)
  res <- cor_vec(dat, sigma_vec = sigma_vec, noise = e)

  combn_mat <- combn(5, 2)
  vec <- apply(combn_mat, 2, function(x){
    mean(dat[,x[1]]*dat[,x[2]]*e)/(sigma_vec[x[1]]*sigma_vec[x[2]])
  })

  expect_true(all(sort(res) == sort(vec)))
})

##########################

## cck is correct

test_that("cck works", {
  set.seed(10)
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

  dat <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = cov_mat)
  g <- row_difference_closure(1,4,4)
  res <- cck(dat, g = g)

  expect_true(length(res) == 3)
  expect_true(is.list(res))
  expect_true(res$pval <= 1)
  expect_true(res$pval >= 0)
})

test_that("cck has sensible p-values", {
  set.seed(10)
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5

  dat <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = cov_mat)
  set.seed(10)
  g1 <- row_difference_closure(1,4,4)
  res1 <- cck(dat, g = g1)

  set.seed(10)
  g2 <- row_difference_closure(1,2,4)
  res2 <- cck(dat, g = g2)

  expect_true(res2$pval < res1$pval)
  expect_true(res1$t0 <= res2$t0)
})

test_that("cck has the right null distribution", {
  trials <- 50
  res_vec <- rep(NA, trials); res_vec2 <- rep(NA, trials)
  cov_mat <- diag(4)
  cov_mat[c(2:3), c(1,4)] <- 0.5;  cov_mat[c(1,4), c(2:3)] <- 0.5
  g1 <- row_difference_closure(1,4,4)
  g2 <- row_difference_closure(1,2,4)

  for(i in 1:trials){
    set.seed(i)
    dat <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = cov_mat)
    set.seed(10)
    res_vec[i] <- cck(dat, g = g1, trials = 50)$pval
  }

  for(i in 1:trials){
    set.seed(i)
    dat <- MASS::mvrnorm(100, mu = rep(0,4), Sigma = cov_mat)
    set.seed(10)
    res_vec2[i] <- cck(dat, g = g2, trials = 50)$pval
  }

  expect_true(sum(abs(quantile(res_vec) - seq(0, 1, length.out = 5)))
              <= sum(abs(quantile(res_vec2) - seq(0, 1, length.out = 5))))

})

########################

