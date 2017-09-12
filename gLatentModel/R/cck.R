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
  doMC::registerDoMC(cores = cores)

  sigma_vec <- apply(dat, 2, stats::sd)
  psi <- translate(dat, sigma_vec = sigma_vec)
  t0 <- g(psi)

  t_boot <- rep(NA, trials)

  func <- function(i){
    e <- stats::rnorm(n)
    g(translate(dat, sigma_vec = sigma_vec, noise_vec = e),
      average_vec = psi*sum(e)/n)
  }

  i <- 1 #debugging purposes
  t_boot <- as.numeric(foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                                          func(i)))

  pval <- length(which(n^(1/2)*abs(t_boot - t0) >= n^(1/2)*abs(t0)))/trials

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
  dat <- scale(dat, scale = F)

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
    max(abs(new_vec[idx1] - new_vec[idx2]))
  }
}
