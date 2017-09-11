rm(list=ls())
devtools::install_github("linnylin92/gLatentModel", ref = "cck",
                         subdir = "gLatentModel", force = T)

source("../main/simulation_base.R")
source("../main/simulation_helper.R")

library(gLatentModel)
library(cord)
set.seed(10)
