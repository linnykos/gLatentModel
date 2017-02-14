rm(list=ls())
load("../results/results_standard.RData")
#load("../results/results_fragment.RData")

glatent_mat <- sapply(glatent_res, function(x){apply(x, 1, mean)})
glatent_sbm_mat <- sapply(glatent_sbm_res, function(x){apply(x, 1, mean)})
hclust_mat <- sapply(hclust_res, function(x){apply(x, 1, mean)})
sbm_mat <- sapply(sbm_res, function(x){apply(x, 1, mean)})

glatent_sd <- sapply(glatent_res, function(x){apply(x, 1, sd)})
glatent_sbm_sd <- sapply(glatent_res, function(x){apply(x, 1, sd)})
hclust_sd <- sapply(hclust_res, function(x){apply(x, 1, sd)})
sbm_sd <- sapply(sbm_res, function(x){apply(x, 1, sd)})

name.vec <- c("Max", "Forbenius", "Spectral", "L1", "Graph Edge Distance",
              "Partnership Difference", "Jaccard Difference")

mean_list <- list(glatent_mat, glatent_sbm_mat, hclust_mat, sbm_mat)
sd_list <- list(glatent_sd, glatent_sbm_sd, hclust_sd, sbm_sd)

par(mfrow = c(1, 4))
for(i in 1:4){
  xlim = c(min(paramMat[,3]), max(paramMat[,3]))
  ylim = c(min(sapply(1:3, function(x){min(mean_list[[x]][i,] - sd_list[[x]][i,])})),
    max(sapply(1:3, function(x){max(mean_list[[x]][i,] + sd_list[[x]][i,])})))

  plot(NA, xlim = xlim, ylim = ylim, xlab = "n", ylab = "Error",
       main = name.vec[i])
  for(j in 1:4){
    offset <- rnorm(1)
    arrows(paramMat[,3] + offset, mean_list[[j]][i,]-sd_list[[j]][i,],
           paramMat[,3] + offset, mean_list[[j]][i,]+sd_list[[j]][i,],
           code = 3, angle = 90, length = 0.02, col = j, lty = 2)
  }

  for(j in c(4,3,1,2)){
    lines(paramMat[,3], mean_list[[j]][i,], col = j, lwd = 4)
  }

  legend("topright", c("GLatent", "GLatent SBM", "Hier.", "SBM"), fill = 1:4)
}

par(mfrow = c(1, 2))
for(i in 6:7){
  xlim = c(min(paramMat[,3]), max(paramMat[,3]))
  ylim = c(min(sapply(1:3, function(x){min(mean_list[[x]][i,] - sd_list[[x]][i,])})),
           max(sapply(1:3, function(x){max(mean_list[[x]][i,] + sd_list[[x]][i,])})))

  plot(NA, xlim = xlim, ylim = ylim, xlab = "Error", ylab = "n",
       main = name.vec[i])
  for(j in 1:4){
    offset <- rnorm(1)
    arrows(paramMat[,3] + offset, mean_list[[j]][i,]-sd_list[[j]][i,],
           paramMat[,3] + offset, mean_list[[j]][i,]+sd_list[[j]][i,],
           code = 3, angle = 90, length = 0.02, col = j, lty = 2)
  }

  for(j in c(4,3,1,2)){
    lines(paramMat[,3], mean_list[[j]][i,], col = j, lwd = 4)
  }

  legend("topright", c("GLatent", "GLatent SBM", "Hier.", "SBM"), fill = 1:4)
}
