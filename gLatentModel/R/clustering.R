#' Converts a partition representation from list to vector
#'
#' @param list a list
#'
#' @return a vector
#' @export
.convert_list_to_vector <- function(list){
  stopifnot(is.list(list))

  n <- max(sapply(list, max))
  vec <- rep(0, n)
  for(i in 1:length(list)){
    vec[list[[i]]] <- i
  }

  vec
}

#' Converts a partition representation from vector to a list
#'
#' @param vec a vector
#'
#' @return a list
#' @export
.convert_vector_to_list <- function(vec){
  lapply(unique(vec), function(x){which(vec == x)})
}
