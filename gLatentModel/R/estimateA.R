# estimate_a <- function(theta_mat, c_mat, pure_idx){
#   stopifnot(is.matrix(c_mat), is.numeric(c_mat), all(c_mat == t(c_mat)))
#   stopifnot(is.matrix(theta_mat), is.numeric(theta_mat), all(theta_mat == t(theta_mat)))
#   stopifnot(all(pure_idx <= ncol(c_mat)), all(pure_idx >= 1), all(pure_idx %%1 == 0))
#   stopifnot(all(sort(unique(pure_idx)) == sort(pure_idx)))
#
#   K <- length(pure_idx); d <- ncol(c_mat)
#   idx <- sort(idx); idx_rest <- c(1:d)[-pure_idx]
#
#   c_mat_shuf <- c_mat[c(idx, idx_rest), c(idx, idx_rest)]
#
#   u_mat <- t(sapply(1:(d-K)), function(x){
#     .optim_solver_constrainLS(c_mat_shuf[x+K, c(1:K)], theta_mat)
#   })
#
#   rbind(diag(K), u_mat)
# }
#
estimate_theta <- function(c_mat, pure_idx){
  stopifnot(is.matrix(c_mat), is.numeric(c_mat), all(c_mat == t(c_mat)))
  stopifnot(all(pure_idx <= ncol(c_mat)), all(pure_idx >= 1), all(pure_idx %%1 == 0))
  stopifnot(all(sort(unique(pure_idx)) == sort(pure_idx)))

  c_mat[pure_idx, pure_idx]
}

pure_nodes <- function(c_mat, k){
  order(abs(diag(c_mat)), decreasing = T)[1:k]
}

.optim_solver_constrainLS <- function(vec, mat, z_vec = rep(1/ncol(mat), ncol(mat)),
  lambda_vec = rep(0, ncol(mat)), rho = 1, tol = 1e-4, max_iter = 100){

  iter <- 1
  beta_vec_prev <- rep(0, ncol(mat))

  beta_mat <- matrix(0, ncol(mat), max_iter)
  z_mat <- matrix(0, ncol(mat), max_iter)
  lambda_mat <- matrix(0, ncol(mat), max_iter)

  while(iter < max_iter){
    beta_vec_new <- .ADMM_solver_primal(vec, mat, lambda_vec, z_vec, rho)
    z_vec <- .ADMM_solver_dual(beta_vec_new, rho, lambda_vec)
    lambda_vec <- lambda_vec + rho*(beta_vec_new - z_vec)

    if(.l2norm(beta_vec_new - beta_vec_prev) <= tol) break()
    beta_vec_old <- beta_vec_new
    iter <- iter + 1

    beta_mat[,iter] <- beta_vec_new
    z_mat[,iter] <- z_vec
    lambda_mat[,iter] <- lambda_vec
  }

  #list(beta = beta_mat, z = z_mat, lambda = lambda_mat)
  beta_vec_new
}

.ADMM_solver_primal <- function(vec, mat, lambda_vec, z_vec, rho){
  K <- ncol(mat)

  res <- solve(-2*t(mat)%*%mat + diag(K))%*%(-lambda_vec + rho * z_vec)
  as.numeric(res/sum(res))
}

.ADMM_solver_dual <- function(beta_vec, rho, lambda_vec){
  res <- 2/rho * (lambda_vec + rho*beta_vec/2)
  res[res < 0] <- 0
  as.numeric(res)
}
