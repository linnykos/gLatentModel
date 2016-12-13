reestimate_gamma <- function(dat, group_list){
  n <- ncol(dat)

  sapply(group_list, function(x){
    var_vec <- apply(dat[x,], 2, var)
    sum(var_vec)/(n * length(x))
  })
}
