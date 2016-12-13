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

  #beta_mat <- matrix(0, ncol(mat), max_iter)
  #z_mat <- matrix(0, ncol(mat), max_iter)
  #lambda_mat <- matrix(0, ncol(mat), max_iter)

  while(iter < max_iter){
    beta_vec_new <- .ADMM_solver_primal(vec, mat, lambda_vec, z_vec, rho)
    z_vec <- .ADMM_solver_dual(beta_vec_new, rho, lambda_vec)
    lambda_vec <- lambda_vec + rho*(beta_vec_new - z_vec)

    if(.l2norm(beta_vec_new - beta_vec_prev) <= tol) break()
    beta_vec_old <- beta_vec_new
    iter <- iter + 1

    #beta_mat[,iter] <- beta_vec_new
    #z_mat[,iter] <- z_vec
    #lambda_mat[,iter] <- lambda_vec
  }

  #list(beta = beta_mat, z = z_mat, lambda = lambda_mat)
  z_vec
}

.ADMM_solver_primal <- function(vec, mat, lambda_vec, z_vec, rho){
  K <- ncol(mat)

  res <- solve(2*t(mat)%*%mat + rho*diag(K))%*%(-lambda_vec + rho * z_vec +
      2*t(mat)%*%vec)
  #res
  as.numeric(res/sum(res))
}

.ADMM_solver_dual <- function(beta_vec, rho, lambda_vec, depth = 1,
  max_depth = 3, grid_length = 100, start = -10, end = 10,
  warnings = 0, tol = 1e-4, multiplier_range = 10){

  if(warnings > max_depth) stop("Error in solving for A dual. Cannot calibrate.")

  gamma_grid <- seq(start, end, length.out = grid_length)

  dual_mat <- sapply(gamma_grid, function(x){
    sapply(beta_vec + 1/rho*(lambda_vec + x), function(y){max(y,0)})
  })
  dual_sum <- apply(dual_mat, 2, sum)

  res <- .select_minimum_index(dual_sum, gamma_grid)

  if(max(dual_sum) < 1-tol | min(dual_sum) > 1+tol){
    #cannot find right u_sum
    dif <- multiplier_range*abs(end - start)
    .ADMM_solver_dual(beta_vec, rho, lambda_vec, start = start-dif,
      end = end+dif, warnings = warnings+1)
  } else if(depth < max_depth & abs(dual_sum[res$idx] - 1) > tol){
    #go to next depth
    .ADMM_solver_dual(beta_vec, rho, lambda_vec, start = res$start,
      end = res$end, depth = depth+1, warnings = 0)
  } else{
    #return
    as.numeric(dual_mat[,res$idx])
  }
}

.select_minimum_index <- function(vec, gamma_grid){
  stopifnot(length(vec) == length(gamma_grid))

  len <- length(gamma_grid)
  idx <- which.min(abs(vec - 1))
  val <- vec[idx]
  dif <- gamma_grid[2] - gamma_grid[1]

  uniq_vec <- sort(unique(vec))
  uniq_idx <- which(uniq_vec == val)

  tmp <- which(vec == uniq_vec[max(1, uniq_idx - 1)])
  start <- gamma_grid[tmp[length(tmp)]] - dif
  tmp <- which(vec == uniq_vec[min(len, uniq_idx + 1)])
  end <- gamma_grid[tmp[1]] + dif

  list(idx = idx, start = start, end = end)
}
