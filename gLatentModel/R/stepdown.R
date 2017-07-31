stepdown <- function(dat, g_list, translate = cor_vec, alpha = 0.05, trials = 100){
  n <- nrow(dat); len <- length(g_list)

  sigma_vec <- apply(dat, 2, stats::sd)
  t0_vec <- sapply(1:len, function(x){g_list[[x]](translate(dat, sigma_vec = sigma_vec))})
  idx_in <- rep(TRUE, length(g_list))

  round <- 1
  while(TRUE){

    t_boot <- rep(NA, trials)
    for(i in 1:trials){
      set.seed(round*10*i)
      e <- stats::rnorm(n)
      t_boot_vec <- sapply(1:len, function(x){
        if(idx_in[x]) {
          res <- g_list[[x]](translate(dat, sigma_vec = sigma_vec, noise_vec = e))
          n^(1/2)*abs(res - t0_vec[x])
        } else NA})
      if(all(is.na(t_boot_vec))) return(numeric(0)) else t_boot[i] <- max(t_boot_vec, na.rm = T)
    }

    cutoff <- quantile(t_boot, 1-alpha, na.rm = T)
    idx <- intersect(which(n^(1/2)*abs(t0_vec) >= cutoff), which(idx_in))

    if(length(idx) == 0) break()

    idx_in[idx] <- FALSE
    round <- round+1
  }

  which(idx_in)
}
