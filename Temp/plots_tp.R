rm(list = ls())

# Required packages
library(twopiece)
library(numDeriv)

# Routines
#source("C:/Users/Javier/Documents/GitHub/MOOMIN/routines/routines_tp.R")
source("~/Documents/GitHub/MOOMIN/routines/routines_tp.R")

################################################################################
# Grid of points on (-3,0) for evaluations
################################################################################
eps <- round(seq(-3,0,length.out = 90), digits = 2)

################################################################################
# Evaluations
################################################################################

# twopiece normal
disctpn <- discrepancy_min_tpn(eps)
ptpn <- unprior_min_tpn(eps)

# twopiece logistic
disctplogis <- discrepancy_min_tplogis(eps)
ptplogis <- unprior_min_tplogis(eps)

# twopiece Laplace
disctplap <- discrepancy_min_tplap(eps)
ptplap <- unprior_min_tplap(eps)

################################################################################
# Data preparation
################################################################################

# Grid
eps_neg  <- eps
eps_pos  <- -rev(eps_neg[-length(eps_neg)])

# twopiece normal
disctpn_neg <- disctpn
disctpn_pos <-  rev(disctpn_neg[-length(disctpn_neg)])

sdisctpn_neg <- -disctpn
sdisctpn_pos <-  -rev(sdisctpn_neg[-length(sdisctpn_neg)])

ptpn_neg <- ptpn
ptpn_pos <-  rev(ptpn[-length(ptpn)])

# twopiece logistic
disctplogis_neg <- disctplogis
disctplogis_pos <-  rev(disctplogis_neg[-length(disctplogis_neg)])

sdisctplogis_neg <- -disctplogis
sdisctplogis_pos <-  -rev(sdisctplogis_neg[-length(sdisctplogis_neg)])

ptplogis_neg <- ptplogis
ptplogis_pos <-  rev(ptplogis[-length(ptplogis)])

# twopiece Laplace
disctplap_neg <- disctplap
disctplap_pos <-  rev(disctplap_neg[-length(disctplap_neg)])

sdisctplap_neg <- -disctplap
sdisctplap_pos <-  -rev(sdisctplap_neg[-length(sdisctplap_neg)])

ptplap_neg <- ptplap
ptplap_pos <-  rev(ptplap[-length(ptplap)])

#---------------------------
# full grid and values
#---------------------------

# Grid
eps_full  <- c(eps_neg, eps_pos)

# twopiece normal
disctpn_full <- c(disctpn_neg, disctpn_pos)
sdisctpn_full <- c(sdisctpn_neg, sdisctpn_pos)
ptpn_full <- c(ptpn_neg, ptpn_pos)

# twopiece logistic
disctplogis_full <- c(disctplogis_neg, disctplogis_pos)
sdisctplogis_full <- c(sdisctplogis_neg, sdisctplogis_pos)
ptplogis_full <- c(ptplogis_neg, ptplogis_pos)

# twopiece Laplace
disctplap_full <- c(disctplap_neg, disctplap_pos)
sdisctplap_full <- c(sdisctplap_neg, sdisctplap_pos)
ptplap_full <- c(ptplap_neg, ptplap_pos)

################################################################################
# Plots
################################################################################

indbad_tpn <- c(which(ptpn_full>0.1),which(abs(eps_full) == 2.87))
# twopiece normal
pdf("disc_min_tpn.pdf", width = 8, height = 6)
plot(eps_full,disctpn_full, lwd = 2, main = "two-piece normal", type = "l",
      xlab = expression(epsilon), ylab = "Discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("sdisc_min_tpn.pdf", width = 8, height = 6)
plot(eps_full,sdisctpn_full, lwd = 2, main = "two-piece normal", type = "l",
     xlab = expression(epsilon), ylab = "Signed discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("prior_min_tpn.pdf", width = 8, height = 6)
plot(eps_full[-indbad_tpn],ptpn_full[-indbad_tpn], lwd = 2, main = "two-piece normal", type = "l",
     xlab = expression(epsilon), ylab = "Prior", cex.axis = 1.5, cex.lab = 1.5)
dev.off()



# twopiece logistic
ind_tplogis_bad <- which(ptplogis_full>0.1)
pdf("disc_min_tplogis.pdf", width = 8, height = 6)
plot(eps_full,disctplogis_full, lwd = 2, main = "two-piece logistic", type = "l",
     xlab = expression(epsilon), ylab = "Discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("sdisc_min_tplogis.pdf", width = 8, height = 6)
plot(eps_full,sdisctplogis_full, lwd = 2, main = "two-piece logistic", type = "l",
     xlab = expression(epsilon), ylab = "Signed discrepancy", cex.axis = 1.5, cex.lab = 1.5)
dev.off()
pdf("prior_min_tplogis.pdf", width = 8, height = 6)
plot(eps_full[-ind_tplogis_bad],ptplogis_full[-ind_tplogis_bad], lwd = 2, main = "two-piece logistic", type = "l",
     xlab = expression(epsilon), ylab = "Prior", cex.axis = 1.5, cex.lab = 1.5)
dev.off()


# twopiece Laplace
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