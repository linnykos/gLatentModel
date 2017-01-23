##formula from Eq (19) of https://arxiv.org/pdf/1606.05100v1.pdf

.estimate_cov_noise <- function(dat, num_subsample = NA){
  stopifnot(is.matrix(dat), is.numeric(dat))
  dat <- scale(dat, center = T, scale = F)
  p <- ncol(dat); n <- nrow(dat)

  sapply(1:p, function(x){
    ne1 <- .neighbor_noise_estimation(dat, x, num_subsample)
    ne2 <- .neighbor_noise_estimation(dat, c(x, ne1), num_subsample)

    t(dat[,x] - dat[,ne1])%*%(dat[,x] - dat[,ne2])/n
  })
}

.l2norm <- function(vec){
  stopifnot(is.numeric(vec))
  sqrt(sum(vec^2))
}

.v_innerproduct <- function(dat, a, b, num_subsample = NA){
  stopifnot(is.numeric(dat), is.matrix(dat), all(c(a,b) <= ncol(dat)))

  p <- ncol(dat)
  idx <- c(1:p)[-c(a,b)]; idx_comb <- utils::combn(idx, 2)
  if(ncol(idx_comb) > 10 & !is.na(num_subsample))
    idx_comb <- idx_comb[,sample(1:ncol(idx_comb))[1:num_subsample]]
  vec_diff <- dat[,a] - dat[,b]

  res <- apply(idx_comb, 2, function(x){
    vec <- dat[,x[1]] - dat[,x[2]]
    if(all(vec  == 0)) 0 else t(vec_diff)%*%(vec/.l2norm(vec))
  })

  max(res, na.rm = T)
}

.neighbor_noise_estimation <- function(dat, a.vec, num_subsample = NA){
  stopifnot(is.numeric(dat), is.matrix(dat), all(a.vec <= ncol(dat)))

  p <- ncol(dat)
  idx <- c(1:p)[-a.vec]
  res <- sapply(idx, .v_innerproduct, dat = dat, a = a.vec[1],
    num_subsample = num_subsample)

  idx[which.min(res)]
}
