rm(list=ls())

library(sn)
library(numDeriv)

source("~/Documents/GitHub/OMOM/R/routines/routines_sn.R")

#source("C:/Users/Javier/Documents/GitHub/OMOM/R/routines/routines_sn.R")

#######################################################################
# Prior based on minimum discrepancy from baseline
#######################################################################

pdf("dis_sn_min.pdf", width = 8, height = 6)
curve(discrepancy_min,-20,20,n=200, lwd = 2,
      xlab = expression(lambda), ylab = "Discrepancy", cex.axis = 1.5, cex.lab = 1.5)
abline(v = -1, col = "red", lwd = 2)
abline(v = 1, col = "red", lwd = 2)
dev.off()

pdf("sdis_sn_min.pdf", width = 8, height = 6)
curve(sdiscrepancy_min,-20,20,n=200, lwd = 2,
      xlab = expression(lambda), ylab = "Signed discrepancy", cex.axis = 1.5, cex.lab = 1.5)
abline(v = -1, col = "red", lwd = 2)
abline(v = 1, col = "red", lwd = 2)
dev.off()

pdf("prior_sn_min.pdf", width = 8, height = 6)
curve(prior_min,-20,20,n=300, lwd = 2,
      xlab = expression(lambda), ylab = "Prior", cex.axis = 1.5, cex.lab = 1.5)
dev.off()

# Mode of the prior
optimize(prior_min, interval =  c(0,3), maximum = TRUE)


pdf("prior_sn_min.pdf", width = 8, height = 6)
curve(prior_min,-20,20,n=200, lwd = 2,
      xlab = expression(lambda), ylab = "Prior", cex.axis = 1.5, cex.lab = 1.5)
curve(tprior_app_min,-20,20, n = 200, lwd = 2, lty = 2, col = "red", add = TRUE)
legend("topright", legend = c("Prior", "Approximation"), lty = c(1,2), col = c("black","red"), lwd = c(2,2))
dev.off()
