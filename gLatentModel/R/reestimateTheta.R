reestimate_theta <- function(cov_mat, gamma_mat, a_mat){
  stopifnot(ncol(gamma_mat) == nrow(gamma_mat))
  stopifnot(nrow(a_mat) == nrow(cov_mat), nrow(a_mat) == nrow(gamma_mat))
  stopifnot(all(as.numeric(a_mat) %in% c(0,1)))

  a_inv <- solve(t(a_mat) %*% a_mat)
  a_inv %*% t(a_mat) %*% (cov_mat - gamma_mat) %*% a_mat %*% a_inv
}
