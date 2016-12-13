load("results.RData")

glatent_mat <- sapply(glatent_res, function(x){apply(x, 1, mean)})

hclust_mat <- sapply(hclust_res, function(x){apply(x, 1, mean)})

sbm_mat <- sapply(sbm_res, function(x){apply(x, 1, mean)})

glatent_sd <- sapply(glatent_res, function(x){apply(x, 1, sd)})

hclust_sd <- sapply(hclust_res, function(x){apply(x, 1, sd)})

sbm_sd <- sapply(sbm_res, function(x){apply(x, 1, sd)})

name.vec <- c("Max", "Forbenius", "Spectral", "L1")

mean_list <- list(glatent_mat, hclust_mat, sbm_mat)
sd_list <- list(glatent_sd, hclust_sd, sbm_sd)

par(mfrow = c(1, 4))
for(i in 1:4){
  xlim = c(min(paramMat[,4]), max(paramMat[,4]))
  ylim = c(min(sapply(mean_list, function(x){min(x[i,])})),
    max(sapply(mean_list, function(x){max(x[i,])})))

  plot(NA, xlim = xlim, ylim = ylim)
  for(j in 1:3){
    lines(paramMat[,4], mean_list[[j]][i,], col = j+1, lwd = 3)
  }
  lines()
}

par(mfrow = c(1,1))
xlim = c(min(paramMat[,4]), max(paramMat[,4]))
ylim = c(min(sapply(mean_list, function(x){min(x[5,])})),
  max(sapply(mean_list, function(x){max(x[5,])})))

plot(NA, xlim = xlim, ylim = ylim)
for(j in 1:3){
  lines(paramMat[,4], mean_list[[j]][5,], col = j+1, lwd = 3)
}


