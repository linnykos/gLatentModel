n <- 500
d <- 100; d4 <- d/4; stopifnot(d %% 4 == 0)
cov_mat <- lapply(1:4, function(x){
  if(x == 1){
    matrix(rep(rep(c(.6, 0, 0, 0.7), each = d4), times = d4), ncol = d4, nrow = d)
  } else if (x == 2){
    matrix(rep(rep(c(0, .9, .5, 0), each = d4), times = d4), ncol = d4, nrow = d)
  } else if (x == 3){
    matrix(rep(rep(c(0, .5, .4, .7), each = d4), times = d4), ncol = d4, nrow = d)
  } else {
    matrix(rep(rep(c(.7, 0, .7, .9), each = d4), times = d4), ncol = d4, nrow = d)
  }
})
cov_mat <- do.call(cbind, cov_mat)

#symmetrize
diag(cov_mat) <- 1
cov_mat <- (cov_mat + t(cov_mat))/2
eig <- eigen(cov_mat)
eig$values[eig$values < 0] <- 0
cov_mat <- eig$vectors%*%diag(eig$values)%*%t(eig$vectors)
cov_mat <- diag(1/diag(cov_mat)^(1/2))%*%cov_mat%*%diag(1/diag(cov_mat)^(1/2))

# cov_mat_simplified <- matrix(sapply(c(0:3)*d4+1, function(x){
#   cov_mat[c(0:3)*d4+2,x]
# }), ncol = 4, nrow = 4)

set.seed(10)
dat <- MASS::mvrnorm(n, mu = rep(0,d), Sigma = cov_mat)

# clockwise90 = function(a) { t(a[nrow(a):1,]) }
# image(clockwise90(cov(dat)))

## use cck
combn_mat <- combn(d, 2)
g_list <- lapply(1:ncol(combn_mat), function(x){
  gLatentModel::row_difference_closure(combn_mat[1,x], combn_mat[2,x], d)})
cck_idx <- gLatentModel::stepdown(dat, g_list)
cck_res <- gLatentModel::connected_components(d, combn_mat[,idx])

## use hierarchical clust
dd <- as.dist((1 - cor(dat))/2)
obj <- hclust(dd, method = "single")
hier_res <- cutree(obj, 4)

save.image("results.RData")

