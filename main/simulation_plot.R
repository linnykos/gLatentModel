rm(list=ls())

load("../results/results.RData")
source("../main/simulation_helper.R")

calculate_error <- function(res, method, true_partition, d = 40){
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

true_partition <- rep(1:4, each = 10)
res_list_variation <- calculate_error(res, "variation", true_partition)
res_list_hamming_covariance <- calculate_error(res, "hamming_covariance", true_partition)
res_list_hamming <- calculate_error(res, "hamming", true_partition)
