simulationGenerator <- function(rule, paramMat, criterion, trials,
                                cores = NA, filename = NA){

  if(!is.na(cores)) doMC::registerDoMC(cores = cores)

  res <- vector(list, nrow(paramMat))
  names(res) <- sapply(1:nrow(paramMat), function(x){
    paste0(paramMat[x,], collapse = "-")})

  for(i in 1:nrow(paramMat)){
    cat(paste0("Row ", x, " started!\n"))

    fun <- function(y, verbose = F){if(verbose) print(y); set.seed(y); criterion(rule(paramMat[x,]), paramMat[x,])}
    if(is.na(cores)){
      res[[i]] <- sapply(1:trials, fun, verbose = T)
    } else {
      res[[i]] <- .adjustFormat(foreach::"%dopar%"(foreach::foreach(trial = 1:trials), fun(trial)))
    }

    if(!is.na(filename)) save(res, file = filename)
  }

  res
}

.adjustFormat <- function(lis){
  len <- sapply(lis, length)
  if(length(unique(len)) != 1) return(lis)

  ncol <- unique(len)
  if(length(ncol) == 1 && ncol == 1) return(as.numeric(unlist(lis)))
  do.call(cbind, lis)
}
