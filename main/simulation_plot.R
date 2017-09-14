rm(list=ls())

load("../main/results.RData")
source("../main/simulation_helper.R")

calculate_error <- function(res, method, true_partition, d = 20){
  lapply(sort(unique(param_mat[,2])), function(x){
    true_cov <- generate_matrix(d = d, strength = x)
    true_cov <- clean_matrix(true_cov)

    idx <- which(param_mat[,2] == x)

    sapply(idx, function(x){
      sapply(1:nrow(res[[x]]), function(y){
        vec <- sapply(1:ncol(res[[x]]), function(z){
          partition_distance(.convert_list_to_vector(res[[x]][y,z][[1]]),
                             true_partition, method, true_cov)
        })

        mean(vec)
      })
    })
  })
}

true_partition <- rep(1:4, each = 5)
res_list_variation <- calculate_error(res, "variation", true_partition)
res_list_hamming_covariance <- calculate_error(res, "hamming_covariance", true_partition)
res_list_hamming <- calculate_error(res, "hamming", true_partition)

res <- list(res_list_variation, res_list_hamming_covariance)
label_vec <- c("Variation distance", "Distance w.r.t. covariance")

par(mfrow = c(2,4), mar = c(3,3,1,1))
for(i in 1:2){
  labrange <- c(0, max(as.numeric(unlist(res[[i]]))))
  for(k in 1:4){
    plot(NA, xlim = c(0, max(n_seq)), ylim = labrange, main = paste0(label_vec[i], ": SNR ", 1-strength_seq[k]))
    for(j in 1:2){
      lines(y = res[[i]][[k]][j,], x = n_seq, col = j, lwd = 2)
    }
  }
}



