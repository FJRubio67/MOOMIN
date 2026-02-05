rm(list=ls())

library(sn)
library(numDeriv)

source("~/Documents/GitHub/MOOMIN/routines/routines_slogis.R")

#######################################################################
# Prior based on minimum discrepancy from baseline
#######################################################################

pdf("dis_slogis.pdf", width = 8, height = 6)
curve(discrepancy_min_slogis,-20,20,n=201, lwd = 2,
      xlab = expression(lambda), ylab = "Discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

pdf("sdis_slogis.pdf", width = 8, height = 6)
curve(sdiscrepancy_min_slogis,-20,20,n=201, lwd = 2,
      xlab = expression(lambda), ylab = "Signed discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

pdf("prior_slogis.pdf", width = 8, height = 6)
curve(unprior_min_slogis,-20,20,n=201, lwd = 2,
      xlab = expression(lambda), ylab = "Prior", cex.axis = 1.5, cex.lab = 1.5)
dev.off()


