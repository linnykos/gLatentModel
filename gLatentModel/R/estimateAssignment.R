estimate_assignment <- function(cov_latent, cov_denoise, pure_idx){
  stopifnot(is.matrix(cov_denoise), is.numeric(cov_denoise), all(cov_denoise == t(cov_denoise)))
  stopifnot(is.matrix(cov_latent), is.numeric(cov_latent), all(cov_latent == t(cov_latent)))
  stopifnot(all(pure_idx <= ncol(cov_denoise)), all(pure_idx >= 1), all(pure_idx %%1 == 0))
  stopifnot(all(sort(unique(pure_idx)) == sort(pure_idx)))

  K <- length(pure_idx); d <- ncol(cov_denoise)
  pure_idx <- sort(pure_idx); idx_rest <- c(1:d)[-pure_idx]

  cov_denoise_shuf <- cov_denoise[c(pure_idx, idx_rest), c(pure_idx, idx_rest)]

  assignment_mat <- t(sapply(1:(d-K), function(x){
    .optim_solver_constrainLS(cov_denoise_shuf[x+K, c(1:K)], cov_latent)
  }))

  rbind(diag(K), assignment_mat)
}

estimate_cov_latent <- function(cov_denoise, pure_idx){
  stopifnot(is.matrix(cov_denoise), is.numeric(cov_denoise), all(cov_denoise == t(cov_denoise)))
  stopifnot(all(pure_idx <= ncol(cov_denoise)), all(pure_idx >= 1), all(pure_idx %%1 == 0))
  stopifnot(all(sort(unique(pure_idx)) == sort(pure_idx)))

  cov_denoise[pure_idx, pure_idx]
}

group_cluster <- function(assignment_mat){
  K <- ncol(assignment_mat)

  as.list(sapply(1:K, function(x){which(assignment_mat[,x] != 0)}))
}

partition_cluster <- function(assignment_mat){
  K <- ncol(assignment_mat)

  assignment_mat <- t(apply(assignment_mat, 1, function(x){
    idx <- rep(0, K)
    val <- max(x); loc <- which(abs(x - val) < 1e-4)
    if(length(loc) > 1) loc <- loc[sample(1:length(loc), 1)]
    idx[loc] <- 1
    idx
  }))

  plyr::alply(assignment_mat, 2, function(x){which(x != 0)})
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
