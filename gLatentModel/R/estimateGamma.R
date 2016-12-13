##formula from Eq (19) of https://arxiv.org/pdf/1606.05100v1.pdf

estimate_gamma <- function(dat){
  stopifnot(is.matrix(dat), is.numeric(dat))
  dat <- scale(dat, center = T, scale = F)
  p <- ncol(dat); n <- nrow(dat)

  sapply(1:p, function(x){
    ne1 <- .neighbor_gamma_estimation(dat, x)
    ne2 <- .neighbor_gamma_estimation(dat, c(x, ne1))

    t(dat[,x] - dat[,ne1])%*%(dat[,x] - dat[,ne2])/n
  })
}

.l2norm <- function(vec){
  stopifnot(is.numeric(vec))
  sqrt(sum(vec^2))
}

.v_innerproduct <- function(dat, a, b, sample = 10){
  stopifnot(is.numeric(dat), is.matrix(dat), all(c(a,b) <= ncol(dat)))

  p <- ncol(dat)
  idx <- c(1:p)[-c(a,b)]; idx_comb <- utils::combn(idx, 2)
  if(ncol(idx_comb) > 10) idx_comb <- idx_comb[,sample(1:ncol(idx_comb))[1:sample]]
  vec_diff <- dat[,a] - dat[,b]

  res <- apply(idx_comb, 2, function(x){
    vec <- dat[,x[1]] - dat[,x[2]]
    if(all(vec  == 0)) 0 else t(vec_diff)%*%(vec/.l2norm(vec))
  })

  max(res, na.rm = T)
}

.neighbor_gamma_estimation <- function(dat, a.vec){
  stopifnot(is.numeric(dat), is.matrix(dat), all(a.vec <= ncol(dat)))

  p <- ncol(dat)
  idx <- c(1:p)[-a.vec]
  res <- sapply(idx, .v_innerproduct, dat = dat, a = a.vec[1])

  idx[which.min(res)]
}
