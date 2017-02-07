rm(list=ls())

library(devtools)
install_github("linnylin92/gLatentModel", ref = "kevin", subdir = "gLatentModel",
  force = T)

library(gLatentModel); library(psych)
library(foreach); library(doMC)

source("../main/simulation_base.R")
source("../main/simulation_functions.R")
