estimate_a <- function(theta_mat, c_mat, pure_idx){
  stopifnot(is.matrix(c_mat), is.numeric(c_mat), all(c_mat == t(c_mat)))
  stopifnot(is.matrix(theta_mat), is.numeric(theta_mat), all(theta_mat == t(theta_mat)))
  stopifnot(all(pure_idx <= ncol(c_mat)), all(pure_idx >= 1), all(pure_idx %%1 == 0))
  stopifnot(all(sort(unique(pure_idx)) == sort(pure_idx)))

  K <- length(pure_idx); d <- ncol(c_mat)
  pure_idx <- sort(pure_idx); idx_rest <- c(1:d)[-pure_idx]

  c_mat_shuf <- c_mat[c(pure_idx, idx_rest), c(pure_idx, idx_rest)]

  u_mat <- t(sapply(1:(d-K), function(x){
    .optim_solver_constrainLS(c_mat_shuf[x+K, c(1:K)], theta_mat)
  }))

  rbind(diag(K), u_mat)
}

estimate_theta <- function(c_mat, pure_idx){
  stopifnot(is.matrix(c_mat), is.numeric(c_mat), all(c_mat == t(c_mat)))
  stopifnot(all(pure_idx <= ncol(c_mat)), all(pure_idx >= 1), all(pure_idx %%1 == 0))
  stopifnot(all(sort(unique(pure_idx)) == sort(pure_idx)))

  c_mat[pure_idx, pure_idx]
}

pure_nodes <- function(c_mat, k){
  sort(order(abs(diag(c_mat)), decreasing = T)[1:k])
}

group_cluster <- function(a_mat){
  K <- ncol(a_mat)

  as.list(sapply(1:K, function(x){which(a_mat[,x] != 0)}))
}

partition_cluster <- function(a_mat){
  K <- ncol(a_mat)

  a_mat <- t(apply(a_mat, 1, function(x){
    idx <- rep(0, K)
    val <- max(x); loc <- which(abs(x - val) < 1e-4)
    if(length(loc) > 1) loc <- loc[sample(1:length(loc), 1)]
    idx[loc] <- 1
    idx
  }))

  plyr::alply(a_mat, 2, function(x){which(x != 0)})
}

# using slide 27 from https://web.stanford.edu/~boyd/papers/pdf/admm_slides.pdf
# thanks!
.optim_solver_constrainLS <- function(vec, mat, z_vec = rep(1/ncol(mat), ncol(mat)),
  lambda_vec = rep(0, ncol(mat)), rho = 10, tol = 1e-4, max_iter = 100){

  iter <- 1
  beta_vec_prev <- rep(0, ncol(mat))

  while(iter < max_iter){
    beta_vec_new <- .ADMM_solver_primal(vec, mat, lambda_vec, z_vec, rho)
    z_vec <- .ADMM_solver_dual(beta_vec_new + lambda_vec)
    lambda_vec <- lambda_vec + beta_vec_new - z_vec

    if(.l2norm(beta_vec_new - beta_vec_prev) <= tol) break()
    beta_vec_old <- beta_vec_new
    iter <- iter + 1
  }

  z_vec
}

.ADMM_solver_primal <- function(vec, mat, lambda_vec, z_vec, rho){
  K <- ncol(mat)

  res <- solve(2*t(mat)%*%mat + rho*diag(K)) %*% (2*t(mat)%*%vec +
      rho*z_vec - rho*lambda_vec)
  as.numeric(res)
}

.ADMM_solver_dual <- function(y_vec, tol = 1e-4){
  K <- length(y_vec)
  y_sort <- sort(y_vec, decreasing = T)

  res_dif <- sapply(1:K, function(x){
    (sum(y_sort[1:x]) - 1)/x
  })

  res_sum <- sapply(res_dif, function(x){
     sum(sapply(y_vec - x, function(y){max(0, y)}))
  })

  if(all(abs(res_sum - 1) > tol)){
    stop("Error in ADMM Dual to solve A")
  } else {
    idx <- which.min(abs(res_sum - 1))
    sapply(y_vec - res_dif[idx], function(y){max(0, y)})
  }
}
