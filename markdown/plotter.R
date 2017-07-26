create_matrix <- function(glatent_res, glatent_sbm_res, hclust_res, sbm_res){
  glatent_mat <- sapply(glatent_res, function(x){apply(x, 1, mean)})
  glatent_sbm_mat <- sapply(glatent_sbm_res, function(x){apply(x, 1, mean)})
  hclust_mat <- sapply(hclust_res, function(x){apply(x, 1, mean)})
  sbm_mat <- sapply(sbm_res, function(x){apply(x, 1, mean)})

  glatent_sd <- sapply(glatent_res, function(x){apply(x, 1, sd)})
  glatent_sbm_sd <- sapply(glatent_res, function(x){apply(x, 1, sd)})
  hclust_sd <- sapply(hclust_res, function(x){apply(x, 1, sd)})
  sbm_sd <- sapply(sbm_res, function(x){apply(x, 1, sd)})

  name_vec <- c("Max", "Forbenius", "Spectral", "L1", "Graph Edge Distance",
                "Partnership Difference", "Jaccard Difference")

  mean_list <- list(glatent_mat, glatent_sbm_mat, hclust_mat, sbm_mat)
  sd_list <- list(glatent_sd, glatent_sbm_sd, hclust_sd, sbm_sd)

  list(mean_list = mean_list, sd_list = sd_list, name_vec = name_vec)
}

plot_data <- function(res, i_vec, n_vec){
  par(mfrow = c(1, length(i_vec)))

  for(i in i_vec){
    xlim = c(min(n_vec), max(n_vec))
    ylim = c(min(sapply(1:4, function(x){min(res$mean_list[[x]][i,] -
                                               res$sd_list[[x]][i,])})),
             max(sapply(1:4, function(x){max(res$mean_list[[x]][i,] +
                                               res$sd_list[[x]][i,])})))

    plot(NA, xlim = xlim, ylim = ylim, xlab = "n", ylab = "Error",
         main = res$name_vec[[i]])

    for(j in 1:4){
      offset <- rnorm(1)
      arrows(n_vec + offset, res$mean_list[[j]][i,] - res$sd_list[[j]][i,],
             n_vec + offset, res$mean_list[[j]][i,] + res$sd_list[[j]][i,],
             code = 3, angle = 90, length = 0.02, col = j, lty = 2)
    }

    for(j in c(4,3,1,2)){
      lines(n_vec, res$mean_list[[j]][i,], col = j, lwd = 4)
    }

    legend("topright", c("Cord", "Spectral", "Hier.", "SBM"), fill = 1:4)
  }

  invisible()
}

load("../results/results_standard.RData")

res <- create_matrix(glatent_res, glatent_sbm_res, hclust_res, sbm_res)

pdf("../markdown/fig/standard_estimation.pdf", width = 8, height = 3)
par(mar = c(4,4,1.5,0.5))
plot_data(res, 1:4, paramMat[,3])
graphics.off()

pdf("../markdown/fig/standard_clustering.pdf", width = 8, height = 4)
par(mar = c(4,4,1.5,0.5))
suppressWarnings(plot_data(res, 6:7, paramMat[,3]))
graphics.off()

###################

load("../results/results_mixing.RData")

res <- create_matrix(glatent_res, glatent_sbm_res, hclust_res, sbm_res)

pdf("../markdown/fig/mixing_estimation.pdf", width = 8, height = 3)
par(mar = c(4,4,1.5,0.5))
plot_data(res, 1:4, paramMat[,3])
graphics.off()

pdf("../markdown/fig/mixing_clustering.pdf", width = 8, height = 4)
par(mar = c(4,4,1.5,0.5))
suppressWarnings(plot_data(res, 6:7, paramMat[,3]))
graphics.off()

