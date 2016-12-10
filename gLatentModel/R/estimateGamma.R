##formula from Eq (19) of https://arxiv.org/pdf/1606.05100v1.pdf

estimate_gamma <- function(mat){
  stopifnot(is.matrix(mat), is.numeric(mat))
  mat <- scale(mat, center = T, scale = F)
  p <- ncol(mat); n <- nrow(mat)

  res <- sapply(1:p, function(x){
    ne1 <- .neighbor_gamma_estimation(mat, x)
    ne2 <- .neighbor_gamma_estimation(mat, c(x, ne1))

    t(mat[,x] - mat[,ne1])%*%(mat[,x] - mat[,ne2])/n
  })

  diag(res)
}

.l2norm <- function(vec){
  stopifnot(is.numeric(vec))
  sqrt(sum(vec^2))
}

.v_innerproduct <- function(mat, a, b){
  stopifnot(is.numeric(mat), is.matrix(mat), all(c(a,b) <= ncol(mat)))

  p <- ncol(mat)
  idx <- c(1:p)[-c(a,b)]; idx_comb <- utils::combn(idx, 2)
  vec_diff <- mat[,a] - mat[,b]

  res <- apply(idx_comb, 2, function(x){
    vec <- mat[,x[1]] - mat[,x[2]]
    if(all(vec  == 0)) 0 else t(vec_diff)%*%(vec/.l2norm(vec))
  })

  max(res, na.rm = T)
}

.neighbor_gamma_estimation <- function(mat, a.vec){
  stopifnot(is.numeric(mat), is.matrix(mat), all(a.vec) <= ncol(mat))

  p <- ncol(mat)
  idx <- c(1:p)[-a.vec]
  res <- sapply(idx, .v_innerproduct, mat = mat, a = a.vec[1])

  idx[which.min(res)]
}
