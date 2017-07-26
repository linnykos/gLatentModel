cck <- function(dat, g, translate = cor_vec, alpha = 0.05, trials = 100){
  n <- nrow(dat)

  t0 <- g(translate(dat))

  t_boot <- rep(NA, trials)

  for(i in 1:trials){
    e <- stats::rnorm(n)
    dat_boot <- diag(e)%*%dat
    t_boot[i] <- g(translate(dat_boot))
  }

  pval <- length(which(t_boot >= t0))/trials

  list(pval = pval, quant = quantile(t_boot, 1-alpha), t0 = t0)
}

cor_vec <- function(dat){
  cor_mat <- stats::cor(dat)
  cor_mat[lower.tri(cor_mat)]
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
