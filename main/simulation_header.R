rm(list=ls())

library(devtools)
install_github("linnylin92/gLatentModel", ref = "kevin", subdir = "gLatentModel",
  force = T)

source("simulation_base.R")
