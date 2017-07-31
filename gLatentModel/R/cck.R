cck <- function(dat, g, translate = cor_vec, alpha = 0.05, trials = 100){
  n <- nrow(dat)

  sigma_vec <- apply(dat, 2, stats::sd)
  t0 <- g(translate(dat, sigma_vec = sigma_vec))

  t_boot <- rep(NA, trials)

  for(i in 1:trials){
    e <- stats::rnorm(n)
    t_boot[i] <- g(translate(dat, sigma_vec = sigma_vec, noise_vec = e))
  }

  pval <- length(which(n^(1/2)*abs(t_boot-t0) >= n^(1/2)*abs(t0)))/trials

  list(pval = pval, quant = quantile(t_boot, 1-alpha), t0 = t0)
}

cor_vec <- function(dat, sigma_vec = rep(1, ncol(dat)), noise_vec = rep(1, nrow(dat))){
  n <- nrow(dat)
  dat <- scale(dat, scale = F)

  mat <- diag(1/sigma_vec) %*% (t(dat) %*% diag(noise_vec) %*% dat) %*% diag(1/sigma_vec)
  mat[lower.tri(mat)]/n
}

row_difference_closure <- function(i,j,d){
  stopifnot(i <= d, j <= d, i >= 1, j >= 1, i%%1 == 0, j%%1 == 0)
  tmp <- matrix(0, d, d)
  tmp[i,] <- 1:d; tmp[,i] <- 1:d
  tmp[j,] <- -(1:d); tmp[,j] <- -(1:d)
  tmp[i,j] <- 0; tmp[j,i] <- 0

  vec <- tmp[lower.tri(tmp)]
  idx1 <- order(vec, decreasing = T)[1:(d-2)]
  idx2 <- order(vec, decreasing = F)[1:(d-2)]

  function(vec){
    max(abs(vec[idx1] - vec[idx2]))
  }
}
