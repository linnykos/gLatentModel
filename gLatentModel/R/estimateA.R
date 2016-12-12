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

# .optim_solver_constrainLS <- function(vec, mat, depth = 1, max_depth = 3,
#   grid_length = 100, multiplier_range = 10, shift_idx = 1, start = -10,
#   end = 10, tol = 1e-4, warnings = 0){
#
#   if(warnings > max_depth) stop("Error in solving for A. Cannot calibrate.")
#
#   gamma_grid <- seq(start, end, length.out = grid_length)
#   inv_mat <- solve(t(mat)%*%mat)
#
#   u_mat <- sapply(gamma_grid, .KKT_grid_solver,
#     vec = vec, mat = mat, inv_mat = inv_mat)
#   u_sum <- apply(u_mat, 2, sum)
#   idx <- which.min(abs(u_sum - 1))
#
#   if(max(u_sum) < 1-tol | min(u_sum) > 1+tol){
#     #cannot find right u_sum
#     dif <- multiplier_range*abs(end - start)
#     .optim_solver_constrainLS(vec, mat, start = start-dif,
#       end = end+dif, warnings = warnings+1)
#   } else if(depth < max_depth){
#     #go to next depth
#     dif <- gamma_grid[2] - gamma_grid[1]
#     next_start <- gamma_grid[max(1, idx - shift_idx)] - dif
#     next_end <- gamma_grid[min(grid_length, idx + shift_idx)] + dif
#     .optim_solver_constrainLS(vec, mat, start = next_start,
#       end = next_end, depth = depth+1)
#   } else{
#     #return
#     as.numeric(u_mat[,idx])
#   }
# }

.KKT_grid_solver <- function(gamma, vec, mat, inv_mat){
  K <- ncol(mat)

  res <- ((2*t(vec)%*%mat - gamma)%*%inv_mat)/2
  idx <- which(res > 0)
  if(length(idx) == 0) return(rep(0, K))

  tmp <- (t(vec)%*%mat[,idx] - gamma/2) %*% solve(t(mat[,idx])%*%mat[,idx])
  res <- rep(0, K)
  res[idx] = tmp
  as.numeric(res)
}
