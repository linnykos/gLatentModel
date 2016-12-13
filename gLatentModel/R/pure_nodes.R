pure_nodes <- function(c_mat, k){
  sort(order(abs(diag(c_mat)), decreasing = T)[1:k])
}
