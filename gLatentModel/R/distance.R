partition_distance <- function(partition1, partition2, method = "variation",
                               covariance = NA){
  stopifnot(length(partition1) == length(partition2))
  .check_partition(partition1); .check_partition(partition2)

  if(method == "variation") {.partition_variation(partition1, partition2)
  } else if(method == "hamming_covariance") {
    .partition_hamming_covariance(partition1, partition2, covariance)
  } else {
    stop("Method not found.")
  }
}

.check_partition <- function(vec){
  stopifnot(is.numeric(vec), !is.list(vec), !is.matrix(vec), !is.factor(vec))
  n <- max(vec)
  stopifnot(all(vec %% 1 == 0), all(vec > 0))
  stopifnot(all(1:n %in% vec))

  TRUE
}

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

.hamming_check <- function(pair, partition1, partition2){
  bool1 <- ifelse(partition1[pair[1]] == partition1[pair[2]], TRUE, FALSE)
  bool2 <- ifelse(partition2[pair[1]] == partition2[pair[2]], TRUE, FALSE)

  if(bool1 == bool2) TRUE else FALSE
}

.l2vec <- function(vec){
  sqrt(sum((vec)^2))
}
