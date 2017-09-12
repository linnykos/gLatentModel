#' Applys CCK to the correlations for a function g
#'
#' @param dat nxd numeric matrix
#' @param g function that takes in (d choose 2) numbers
#' @param translate translates each sample of length d into some other vector
#' @param alpha significance value, between 0 and 1
#' @param trials positive integer of bootstrap trials
#' @param cores number of cores to parallelize computation over
#'
#' @return
#' @export
cck <- function(dat, g, translate = cor_vec, alpha = 0.05, trials = 100,
                cores = 2){
  n <- nrow(dat)
  dat <- scale(dat, scale = F)
  doMC::registerDoMC(cores = cores)

  sigma_vec <- apply(dat, 2, stats::sd)
  psi <- translate(dat, sigma_vec = sigma_vec)
  theta <- g(psi)
  t0 <- max(abs(theta))

  func <- function(i){
    idx <- sample(1:n, n, replace = T)
    dat_tmp <- dat[idx,]
    dat_tmp <- scale(dat_tmp, scale = F)
    sigma_tmp <- apply(dat_tmp, 2, stats::sd)
    psi_tmp <- cor_vec(dat_tmp, sigma_vec = sigma_tmp)
    g(psi_tmp)
  }

  i <- 1 #debugging purposes
  theta_boot <- foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                                          func(i))
  t_boot <- sapply(theta_boot, function(x){max(abs(x - theta))})

  pval <- length(which(n^(1/2)*t_boot >= n^(1/2)*t0))/trials

  list(pval = pval, quant = stats::quantile(t_boot, 1-alpha), t0 = t0)
}

#' Converts the data into correlations
#'
#' Allows \code{noise_vec} which multiplies each sample with a certain amount of
#' noise. Also requires the estimated standard deviations (marginal) in
#' \code{sigma_vec}.
#'
#' @param dat nxd numeric matrix
#' @param sigma_vec length d numeric vector
#' @param noise_vec length d numeric vector
#'
#' @return a vector of length (d choose 2).
#' @export
cor_vec <- function(dat, sigma_vec = rep(1, ncol(dat)), noise_vec = rep(1, nrow(dat))){
  n <- nrow(dat)

  mat <- diag(1/sigma_vec) %*% (t(dat) %*% diag(noise_vec) %*% dat) %*% diag(1/sigma_vec)
  mat[lower.tri(mat)]/n
}

#' Generate function to evaluate difference between two rows
#'
#' @param i integer between 1 and d
#' @param j integer between 1 and d
#' @param d positive ingeter
#'
#' @return a function
#' @export
row_difference_closure <- function(i,j,d){
  stopifnot(i <= d, j <= d, i >= 1, j >= 1, i%%1 == 0, j%%1 == 0)
  tmp <- matrix(0, d, d)
  tmp[i,] <- 1:d; tmp[,i] <- 1:d
  tmp[j,] <- -(1:d); tmp[,j] <- -(1:d)
  tmp[i,j] <- 0; tmp[j,i] <- 0

  vec <- tmp[lower.tri(tmp)]
  idx1 <- order(vec, decreasing = T)[1:(d-2)]
  idx2 <- order(vec, decreasing = F)[1:(d-2)]

  function(vec, average_vec = rep(0,length(vec))){
    new_vec <- vec - average_vec
    new_vec[idx1] - new_vec[idx2]
  }
}
