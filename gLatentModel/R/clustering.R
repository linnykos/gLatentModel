.convert_list_to_vector <- function(list){
  n <- max(sapply(list, max))
  vec <- rep(0, n)
  for(i in 1:length(list)){
    vec[list[[i]]] <- i
  }

  vec
}

.convert_vector_to_list <- function(vec){
  lapply(unique(vec), function(x){which(vec == x)})
}
