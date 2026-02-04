rm(list = ls())
library(twopiece)
library(numDeriv)

source("C:/Users/Javier/Documents/GitHub/MOOMIN/routines/routines_tp.R")
#source("~/Documents/GitHub/MOOMIN/routines/routines_tp.R")

eps <- round(sort(c(0,seq(-3,3,length.out = 100))), digits = 2)

disctpn <- discrepancy_min_tpn(eps)
stpn <- sdiscrepancy_min_tpn(eps)
ptpn <- unprior_min_tpn(eps)

disctplogis <- discrepancy_min_tplogis(eps)
stplogis <- sdiscrepancy_min_tplogis(eps)
ptplpgis <- unprior_min_tplogis(eps)

disctplap <- discrepancy_min_tplap(eps)
stplap <- sdiscrepancy_min_tplap(eps)
ptplap <- unprior_min_tplap(eps)



pdf("disc_min_tpn.pdf", width = 8, height = 6)
plot(eps,disctpn, lwd = 2, main = "two-piece normal", type = "l",
      xlab = expression(epsilon), ylab = "Discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("sdisc_min_tpn.pdf", width = 8, height = 6)
plot(eps,stpn, lwd = 2, main = "two-piece normal", type = "l",
     xlab = expression(epsilon), ylab = "Signed discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("prior_min_tpn.pdf", width = 8, height = 6)
plot(eps,ptpn, lwd = 2, main = "two-piece normal", type = "l",
     xlab = expression(epsilon), ylab = "Prior", cex.axis = 1.5, cex.lab = 1.5)
dev.off()





pdf("disc_min_tplogis.pdf", width = 8, height = 6)
plot(eps,disctplogis, lwd = 2, main = "two-piece logistic", type = "l",
     xlab = expression(epsilon), ylab = "Discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("sdisc_min_tplogis.pdf", width = 8, height = 6)
plot(eps,stplogis, lwd = 2, main = "two-piece logistic", type = "l",
     xlab = expression(epsilon), ylab = "Signed discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("prior_min_tplogis.pdf", width = 8, height = 6)
plot(eps,ptplogis, lwd = 2, main = "two-piece logistic", type = "l",
     xlab = expression(epsilon), ylab = "Prior", cex.axis = 1.5, cex.lab = 1.5)
dev.off()



pdf("disc_min_tplap.pdf", width = 8, height = 6)
plot(eps,disctplap, lwd = 2, main = "two-piece Laplace", type = "l",
     xlab = expression(epsilon), ylab = "Discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("sdisc_min_tplap.pdf", width = 8, height = 6)
plot(eps,stplap, lwd = 2, main = "two-piece Laplace", type = "l",
     xlab = expression(epsilon), ylab = "Signed discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("prior_min_tplap.pdf", width = 8, height = 6)
plot(eps,ptplap, lwd = 2, main = "two-piece Laplace", type = "l",
     xlab = expression(epsilon), ylab = "Prior", cex.axis = 1.5, cex.lab = 1.5)
dev.off()