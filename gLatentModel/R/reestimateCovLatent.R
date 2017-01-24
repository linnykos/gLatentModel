.reestimate_cov_latent <- function(cov_mat, gamma_mat, assignment_mat){
  stopifnot(ncol(gamma_mat) == nrow(gamma_mat))
  stopifnot(nrow(cov_mat) == nrow(gamma_mat))

  a_inv <- solve(t(assignment_mat) %*% assignment_mat)
  a_inv %*% t(assignment_mat) %*% (cov_mat - gamma_mat) %*% assignment_mat %*% a_inv
}
