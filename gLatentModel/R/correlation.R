correlation <- function(dat, g, translate = cor_vec, alpha = 0.05, trials = 100,
                cores = 2){
  n <- nrow(dat)
  doMC::registerDoMC(cores = cores)

  sigma_vec <- apply(dat, 2, stats::sd)
  t0 <- g(translate(dat, sigma_vec = sigma_vec))

  t_boot <- rep(NA, trials)

  func <- function(i){
    e <- stats::rnorm(n)
    g(translate(dat, sigma_vec = sigma_vec, noise_vec = e))
  }

  i <- 1 #debugging purposes
  t_boot <- as.numeric(foreach::"%dopar%"(foreach::foreach(i = 1:trials),
                                          func(i)))

  pval <- length(which(n^(1/2)*abs(t_boot-t0) >= n^(1/2)*abs(t0)))/trials

  list(pval = pval, quant = stats::quantile(t_boot, 1-alpha), t0 = t0)
}

.average_squared <- function(){

}

.average_absolute <- function(){

}
