#' Computes the distance between two partitions.
#'
#' \code{partition1} and \code{partition2} must be valid partitions. See
#' \code{.check_partition} to see what this means.
#' You can also pass in a covariance matrix (to take the distance with
#' respect to), but this is only required for certain \code{method}s.
#'
#' @param partition1 a partition vector
#' @param partition2 a partition vector
#' @param method a string
#' @param covariance a PSD matrix
#'
#' @return a numeric
#' @export
partition_distance <- function(partition1, partition2, method = "variation",
                               covariance = NA){
  stopifnot(length(partition1) == length(partition2))
  .check_partition(partition1); .check_partition(partition2)

  if(method == "variation") {.partition_variation(partition1, partition2)
  } else if(method == "hamming") {
    .partition_hamming(partition1, partition2)
  } else if(method == "hamming_covariance") {
    .partition_hamming_covariance(partition1, partition2, covariance)
  } else {
    stop("Method not found.")
  }
}

#' Checks if a partition is valid
#'
#' A partition is valid if it contains only positive integers. A parition
#' of \code{n} samples will be a vector of length \code{n}. If there are
#' \code{k} partitions, then the elements of this vector must range
#' between \code{1} and \code{k} (inclusive).
#'
#' This function will throw an error if \code{vec} is not a valid partition.
#' Otherwise, it will return \code{TRUE}.
#'
#' @param vec a partition vector
#'
#' @return a boolean
#' @export
.check_partition <- function(vec){
  stopifnot(is.numeric(vec), !is.list(vec), !is.matrix(vec), !is.factor(vec))
  n <- max(vec)
  stopifnot(all(vec %% 1 == 0), all(vec > 0))
  stopifnot(all(1:n %in% vec))

  TRUE
}

#' Partition distance with variation of information
#'
#' See \url{https://en.wikipedia.org/wiki/Variation_of_information} for the
#' formula used here.
#'
#' @param partition1 a partition
#' @param partition2 a partition
#'
#' @return a numeric
#' @export
.partition_variation <- function(partition1, partition2){
  n <- length(partition1)
  vec1 <- table(partition1)/n; vec2 <- table(partition2)/n
  num1 <- length(vec1); num2 <- length(vec2)
  tab <- table(partition1, partition2)/n

  grid <- expand.grid(1:num1, 1:num2)

  vec <- apply(grid, 1, function(x){
    if(tab[x[1],x[2]] == 0) return(0)
    -tab[x[1],x[2]]*(log(tab[x[1],x[2]] / vec1[x[1]]) +
                        log(tab[x[1],x[2]] / vec2[x[2]]))
  })

  sum(vec)
}

.partition_hamming <- function(partition1, partition2){
  n <- length(partition1)
  combn_mat <- utils::combn(n, 2)

  vec <- apply(combn_mat, 2, function(x){
    bool <- .hamming_check(x, partition1, partition2)
    if(bool) 0 else 1
  })

  sum(vec)/ncol(combn_mat)
}

#' Hamming partition distance with respect to a covariance matrix
#'
#' This partition first calculates the hamming distance between every
#' pair of samples \code{n}.If the hamming distance is \code{TRUE}, then
#' a distance of \code{0} is assigned to this pair. Otherwise, the
#' L2 distance squared of the distance between the respective rows of the
#' covariance matrix is assigned to this pair. (Note: each vector does not
#' contain the two samples being used.) Then the total distance is the
#' square-root of the sum of all the individual distances.
#'
#' @param partition1 a partition vector
#' @param partition2 a parititon vector
#' @param covariance a PSD matrix
#'
#' @return a numeric
#' @export
.partition_hamming_covariance <- function(partition1, partition2, covariance){
  stopifnot(is.matrix(covariance), all(dim(covariance) > 2), nrow(covariance) == ncol(covariance))

  n <- length(partition1)
  combn_mat <- utils::combn(n, 2)

  vec <- apply(combn_mat, 2, function(x){
    bool <- .hamming_check(x, partition1, partition2)
    if(bool) 0 else .l2vec(covariance[-x,x[1]] - covariance[-x,x[2]])^2
  })

  sqrt(sum(vec))
}

#' Hamming check for two samples
#'
#' If both samples are assigned to the same cluster or different clusters in both partitions,
#' then return \code{TRUE}. Otherwise, return \code{FALSE}.
#'
#' @param pair a pair of integers, from 1 to \code{n}
#' @param partition1 a partition vector
#' @param partition2 a partition vector
#'
#' @return a boolean
#' @export
.hamming_check <- function(pair, partition1, partition2){
  bool1 <- ifelse(partition1[pair[1]] == partition1[pair[2]], TRUE, FALSE)
  bool2 <- ifelse(partition2[pair[1]] == partition2[pair[2]], TRUE, FALSE)

  if(bool1 == bool2) TRUE else FALSE
}

#' The L2 distance
#'
#' @param vec a numeric vector
#'
#' @return a numeric
#' @export
.l2vec <- function(vec){
  sqrt(sum((vec)^2))
}
