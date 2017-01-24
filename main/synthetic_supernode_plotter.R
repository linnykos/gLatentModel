rm(list=ls())
load("../main/results_standard.RData")
#load("../main/results_fragment.RData")

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
  xlim = c(min(paramMat[,3]), max(paramMat[,3]))
  ylim = c(min(sapply(mean_list, function(x){min(x[i,])})),
    max(sapply(mean_list, function(x){max(x[i,])})))

  plot(NA, xlim = xlim, ylim = ylim, xlab = "Error", ylab = "n",
       main = name.vec[i])
  for(j in 1:3){
    lines(paramMat[,3], mean_list[[j]][i,], col = j+1, lwd = 3)
    arrows(paramMat[,3], mean_list[[j]][i,]-sd_list[[j]][i,],
           paramMat[,3], mean_list[[j]][i,]+sd_list[[j]][i,],
           code = 3, angle = 90, length = 0.02)
  }

  legend("topright", c("GLatent", "Hier.", "SBM"), fill = 2:4)
}

par(mfrow = c(1,1))
xlim = c(min(paramMat[,3]), max(paramMat[,3]))
ylim = c(min(sapply(mean_list, function(x){min(x[5,])})),
  max(sapply(mean_list, function(x){max(x[5,])})))

plot(NA, xlim = xlim, ylim = ylim, xlab = "n", ylab = "Clustering Dissimilarity")
for(j in 1:3){
  lines(paramMat[,3], mean_list[[j]][5,], col = j+1, lwd = 3)
  arrows(paramMat[,3], mean_list[[j]][5,]-sd_list[[j]][5,],
         paramMat[,3], mean_list[[j]][5,]+sd_list[[j]][5,],
         code = 3, angle = 90, length = 0.02)
}
legend("topright", c("GLatent", "Hier.", "SBM"), fill = 2:4)

