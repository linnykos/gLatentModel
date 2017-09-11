source("../main/simulation_helper.R")
n <- 100; d <- 100
cov_mat <- generate_matrix(n = n, d = d, strength = 0)
cov_mat <- clean_matrix(cov_mat)

clockwise90 <- function(a) { t(a[nrow(a):1,]) }
redgreen <- colorRampPalette( c("green", "black", "red"), space="rgb")(32)

set.seed(10)
dat <- MASS::mvrnorm(n, mu = rep(0,d), Sigma = cov_mat)
shuff_idx <- sample(1:d)
rev_shuff_idx <- order(shuff_idx)

pdf("../main/data.pdf", height = 6, width = 4)
image(clockwise90(dat), col = redgreen, xaxt = "n", yaxt = "n", ann = F,
      bty = "n")
graphics.off()

pdf("../main/cor.pdf", height = 4, width = 4)
image(clockwise90(cor(dat)), col = rev(heat.colors(50)), xaxt = "n", yaxt = "n", ann = F,
      bty = "n", asp = T, zlim = c(0,1))
graphics.off()

pdf("../main/data_shuff.pdf", height = 6, width = 4)
image(clockwise90(dat[,shuff_idx]), col = redgreen, xaxt = "n", yaxt = "n", ann = F,
      bty = "n")
graphics.off()

pdf("../main/cor_shuff.pdf", height = 4, width = 4)
image(clockwise90(cor(dat[,shuff_idx])), col = rev(heat.colors(50)), xaxt = "n", yaxt = "n", ann = F,
      bty = "n", asp = T, zlim = c(0,1))
graphics.off()


strength_vec <- c(0, 0.3, 0.6, 0.9)
for(i in 1:length(strength_vec)){
  set.seed(10)
  n <- 100; d <- 100
  cov_mat <- generate_matrix(n = n, d = d, strength = strength_vec[i])
  cov_mat <- clean_matrix(cov_mat)
  dat <- MASS::mvrnorm(n, mu = rep(0,d), Sigma = cov_mat)

  pdf(paste0("../main/cor_", i, ".pdf"), height = 4, width = 4)
  image(clockwise90(cor(dat)), col = rev(heat.colors(50)), xaxt = "n", yaxt = "n", ann = F,
        bty = "n", asp = T, zlim = c(0,1))
  graphics.off()
}

clust_list_thres <- list(c(3.8, 2.6, 1.3, 1), c(4, 3.6, 1.6, 1),
                         c(4, 4, 2.7, 1), c(4, 4, 3.9, 1))
error_list_thres <- list(c(0.25, 0.42, 0.69, 0.76), c(0.04, 0.29, 0.66, 0.76),
                         c(0, 0.01, 0.33, 0.76), c(0, 0, 0.07, 0.76))
clust_list_test <- list(c(3.7, 2.9, 1.5, 1), c(4, 3.5, 1.9, 1),
                         c(4, 4, 3.1, 1.8), c(4, 4, 4, 2.5))
error_list_test <- list(c(0.21, 0.49, 0.65, 0.76), c(0.05, 0.3, 0.71, 0.76),
                         c(0, 0.01, 0.21, 0.52), c(0, 0, 0.01, 0.33))

xvec <- c(0.6, 0.42, 0.24, 0.06)

pdf("../main/comparison1.pdf", height = 4, width = 7)
par(mfrow = c(1,2))
plot(NA, xlim = c(0.06, 0.6), ylim = c(1, 4), main = "Average number of\nestimated clusters",
     xlab = "Correlation gap", ylab = "Num. clust")
for(i in 1:1){
  lines(x = xvec, y = clust_list_thres[[i]], col = i, lty = 2, lwd = 2)
  lines(x = xvec, y = clust_list_test[[i]], col = i, lwd = 2)
}

plot(NA, xlim = c(0.06, 0.6), ylim = c(0,1), main = "Average error\nin clustering",
     xlab = "Correlation gap", ylab = "Error in clustering")
for(i in 1:1){
  lines(x = xvec, y = error_list_thres[[i]], col = i, lty = 2, lwd = 2)
  lines(x = xvec, y = error_list_test[[i]], col = i, lwd = 2)
}
dev.off()


pdf("../main/comparison2.pdf", height = 4, width = 7)
par(mfrow = c(1,2))
plot(NA, xlim = c(0.06, 0.6), ylim = c(1, 4), main = "Average number of\nestimated clusters",
     xlab = "Correlation gap", ylab = "Num. clust")
lines(x = xvec, y = clust_list_thres[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = clust_list_test[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = clust_list_thres[[2]], col = rgb(1, 0, 0), lty = 2, lwd = 2)
lines(x = xvec, y = clust_list_test[[2]], col = rgb(1, 0, 0), lwd = 2)

plot(NA, xlim = c(0.06, 0.6), ylim = c(0,1), main = "Average error\nin clustering",
     xlab = "Correlation gap", ylab = "Error in clustering")
lines(x = xvec, y = error_list_thres[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = error_list_test[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = error_list_thres[[2]], col = rgb(1, 0, 0), lty = 2, lwd = 2)
lines(x = xvec, y = error_list_test[[2]], col = rgb(1, 0, 0), lwd = 2)
dev.off()

pdf("../main/comparison3.pdf", height = 4, width = 7)
par(mfrow = c(1,2))
plot(NA, xlim = c(0.06, 0.6), ylim = c(1, 4), main = "Average number of\nestimated clusters",
     xlab = "Correlation gap", ylab = "Num. clust")
lines(x = xvec, y = clust_list_thres[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = clust_list_test[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = clust_list_thres[[2]], col = rgb(1, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = clust_list_test[[2]], col = rgb(1, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = clust_list_thres[[3]], col = rgb(0, 0.7, 0), lty = 2, lwd = 2)
lines(x = xvec, y = clust_list_test[[3]], col = rgb(0, 0.7, 0), lwd = 2)

plot(NA, xlim = c(0.06, 0.6), ylim = c(0,1), main = "Average error\nin clustering",
     xlab = "Correlation gap", ylab = "Error in clustering")
lines(x = xvec, y = error_list_thres[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = error_list_test[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = error_list_thres[[2]], col = rgb(1, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = error_list_test[[2]], col = rgb(1, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = error_list_thres[[3]], col = rgb(0, 0.7, 0), lty = 2, lwd = 2)
lines(x = xvec, y = error_list_test[[3]], col = rgb(0, 0.7, 0), lwd = 2)
dev.off()

pdf("../main/comparison4.pdf", height = 4, width = 7)
par(mfrow = c(1,2))
plot(NA, xlim = c(0.06, 0.6), ylim = c(1, 4), main = "Average number of\nestimated clusters",
     xlab = "Correlation gap", ylab = "Num. clust")
lines(x = xvec, y = clust_list_thres[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = clust_list_test[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = clust_list_thres[[2]], col = rgb(1, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = clust_list_test[[2]], col = rgb(1, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = clust_list_thres[[3]], col = rgb(0.4, 0.7, 0.4, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = clust_list_test[[3]], col = rgb(0.4, 0.7, 0.4, 0.7), lwd = 2)
lines(x = xvec, y = clust_list_thres[[4]], col = rgb(0, 0, 1), lty = 2, lwd = 2)
lines(x = xvec, y = clust_list_test[[4]], col = rgb(0, 0, 1), lwd = 2)

plot(NA, xlim = c(0.06, 0.6), ylim = c(0,1), main = "Average error\nin clustering",
     xlab = "Correlation gap", ylab = "Error in clustering")
lines(x = xvec, y = error_list_thres[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = error_list_test[[1]], col = rgb(0.7, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = error_list_thres[[2]], col = rgb(1, 0.7, 0.7, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = error_list_test[[2]], col = rgb(1, 0.7, 0.7, 0.7), lwd = 2)
lines(x = xvec, y = error_list_thres[[3]], col = rgb(0.4, 0.7, 0.4, 0.7), lty = 2, lwd = 2)
lines(x = xvec, y = error_list_test[[3]], col = rgb(0.4, 0.7, 0.4, 0.7), lwd = 2)
lines(x = xvec, y = error_list_thres[[4]], col = rgb(0, 0, 1), lty = 2, lwd = 2)
lines(x = xvec, y = error_list_test[[4]], col = rgb(0, 0, 1), lwd = 2)
dev.off()


