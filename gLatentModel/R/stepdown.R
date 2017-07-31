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
#' to denote all the hypothesis that passed
#' @export
stepdown <- function(dat, g_list, translate = cor_vec, alpha = 0.05, trials = 100,
                     cores = 2){
  n <- nrow(dat); len <- length(g_list)
  doMC::registerDoMC(cores = cores)

  sigma_vec <- apply(dat, 2, stats::sd)
  psi <- translate(dat, sigma_vec = sigma_vec)
  t0_vec <- sapply(1:len, function(x){g_list[[x]](psi)})
  idx_in <- rep(TRUE, length(g_list))

  func <- function(i){
    set.seed(round*10*i)
    e <- stats::rnorm(n)
    t_boot_vec <- sapply(1:len, function(x){
      if(idx_in[x]) {
        res <- g_list[[x]](translate(dat, sigma_vec = sigma_vec, noise_vec = e),
                           average_vec = psi*sum(e)/n)
      } else NA})
    if(all(is.na(t_boot_vec))) return(numeric(0)) else t_boot[i] <- n^(1/2)*max(t_boot_vec, na.rm = T)
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
