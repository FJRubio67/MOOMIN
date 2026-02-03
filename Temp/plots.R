rm(list = ls())
library(twopiece)
library(numDeriv)

source("C:/Users/Javier/Documents/GitHub/MOOMIN/routines/routines_tp.R")

eps <- seq(-3,3,by = 0.1)

stpn <- sdiscrepancy_min_tpn(eps)
ptpn <- unprior_min_tpn(eps)

stplogis <- sdiscrepancy_min_tplogis(eps)
ptplpgis <- unprior_min_tplogis(eps)

stplap <- sdiscrepancy_min_tplap(eps)
ptplap <- unprior_min_tplap(eps)

