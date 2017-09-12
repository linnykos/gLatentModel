#' Stepdown procedure
#'
#' @param dat nxd numeric matrix
#' @param g_list a list of function that takes in (d choose 2) numbers
#' @param translate translates each sample of length d into some other vector
#' @param alpha significance value, between 0 and 1
#' @param trials positive integer of bootstrap trials
#' @param cores number of cores to parallelize computation over
#'
#' @return a set of indices of integers between 1 to \code{length(g_list)}
#' to denote all the hypothesis that passed (i.e., failed to reject the null)
#' @export
stepdown <- function(dat, g_list, translate = cor_vec, alpha = 0.05, trials = 100,
                     cores = 2){
  n <- nrow(dat); len <- length(g_list)
  dat <- scale(dat, scale = F)
  doMC::registerDoMC(cores = cores)

  sigma_vec <- apply(dat, 2, stats::sd)
  psi <- translate(dat, sigma_vec = sigma_vec)
  theta_list <- lapply(1:len, function(x){g_list[[x]](psi)})
  t0_vec <- sapply(theta_list, function(x){max(abs(x))})
  idx_in <- rep(TRUE, length(g_list))

  func <- function(i){
    set.seed(round*10*i)
    idx <- sample(1:n, n, replace = T)
    dat_tmp <- dat[idx,]
    dat_tmp <- scale(dat_tmp, scale = F)
    sigma_tmp <- apply(dat_tmp, 2, stats::sd)

    t_boot <- sapply(1:len, function(x){
      if(idx_in[x]) {
        theta_boot <- g_list[[x]](translate(dat_tmp, sigma_vec = sigma_tmp))
        max(abs(theta_boot - theta_list[[x]]))
      } else NA})

    if(all(is.na(t_boot))) return(numeric(0)) else n^(1/2)*max(t_boot, na.rm = T)
  }

  round <- 1
  while(TRUE){
    t_boot <- rep(NA, trials)
    i <- 1 #debugging purposes
    t_boot <- as.numeric(foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                                            func(i)))

    cutoff <- stats::quantile(t_boot, 1-alpha, na.rm = T)
    idx <- intersect(which(n^(1/2)*abs(t0_vec) >= cutoff), which(idx_in))

    if(length(idx) == 0) break()

    idx_in[idx] <- FALSE
    round <- round+1
  }

  which(idx_in)
}

#' Returns the number of connected components from a list of edges
#'
#' @param d number of nodes
#' @param edges 2x(number of edges) matrix of integers
#'
#' @return list of the connected components
#' @export
connected_components <- function(d, edges){
  stopifnot(nrow(edges) == 2)
  stopifnot(all(as.numeric(edges) %% 1 == 0), all(as.numeric(edges) <= d),
            all(as.numeric(edges) >= 1))
  stopifnot(d %% 1 == 0, d > 1)

  g <- igraph::graph.empty(n = d, directed = F)
  g <- igraph::add_edges(g, edges = edges)
  res <- igraph::components(g)

  lapply(1:res$no, function(x){which(res$membership == x)})
}
