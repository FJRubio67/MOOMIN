## ----setup, include=FALSE------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------------------------------
rm(list=ls())

# Required packages
library(sn)
library(numDeriv)
library(ggplot2)
library(tidyr)
library(dplyr)
library(knitr)
library(kableExtra)

# Routines
source("~/Documents/GitHub/OMOM/R/routines/routines_sn.R")
#source("C:/Users/Javier/Documents/GitHub/OMOM/R/routines/routines_sn.R")


#  Data: BMI Australian athletes
data(ais)
data <- ais$BMI[ais$sex == "female"]
length(data)

# histogram
hist(data, breaks = 15, xlab = "BMI", ylab = "Density", probability = TRUE,
     cex.lab = 1.5, cex.axis = 1.5, main = "Full data") 
box()

boxplot(data, ylab = "BMI",
     cex.lab = 1.5, cex.axis = 1.5, main = "Full data")


## ------------------------------------------------------------------------------------------------
curve(tprior_app_min, -20, 20, lwd = 2, n = 500, xlab = expression(lambda), ylab = "Density",
      cex.lab = 1.5, cex.axis = 1.5, ylim = c(0,0.2))
curve(d_mom, -20, 20, add = TRUE, col="red", lwd = 2, n = 500, lty = 2)
curve(da_jeff,-20,20, add=TRUE, col="blue", lwd = 2, n = 500, lty = 3)
legend( "topright", legend = c("MOOMIN", "DIMOM", "Jeffreys"), lwd = 2, col = c("black","red","blue"), lty = c(1,2,3))


## ------------------------------------------------------------------------------------------------
# Normality tests
NormTest <- snNormalityTest(data, method = "ILA")

# Posterior probability of sn model based on Jeffreys prior
probsn_J <- NormTest$probsn_J
print(probsn_J)

# Posterior probability of sn model based on DIMOM prior
probsn_MOM <- NormTest$probsn_MOM
print(probsn_MOM)

# Posterior probability of sn model based on MOOMIN prior
probsn_MOOMIN <- NormTest$probsn_MOOMIN
print(probsn_MOOMIN)

# BIC sn model
BICSN <- NormTest$BICSN
print(BICSN)

# BIC normal model
BICN <- NormTest$BICN
print(BICN)

# p-value
psw <- NormTest$pval
print(psw)

# Estimators

EST <- rbind(c(NormTest$MLEN,NA), NormTest$MLESN,
             c(NormTest$MAPN,NA),
             NormTest$MAPSN_J, 
             NormTest$MAPSN_MOM, NormTest$MAPSN_MOOMIN)

rownames(EST) <- c("MLE_N", "MLE_SN",
                  "MAP_N", 
                  "MAPSN_J", 
                  "MAPSN_DIMOM", "MAPSN_MOOMIN")

kable(EST, digits = 3, format = "html") %>%
  kable_styling(full_width = FALSE) %>%
  column_spec(1:ncol(EST), extra_css = "padding-left: 15px; padding-right: 15px;")


## ------------------------------------------------------------------------------------------------
##################################################################################################
# Robust outlier detection
# Location estimated with the median
# Scale estimated with the normalised mean absolute deviation
# https://rpubs.com/FJRubio/Outlier
##################################################################################################


# Normalised Median Absolute Deviation
b <- 1/qnorm(0.75) # Normalisation constant
NMAD <- function(data)  b*median( abs(data - median(data)))

# Number of NMADS from the Median
a <- 3

# NMAD and Median
NMAD.data <- NMAD(data)
med <- median(data)

c(med,NMAD.data)

# Robust outlier detection
Rob.Out.detect <- Vectorize(function(x){
  val <- abs((x - med)/NMAD.data)
  return(ifelse(val>a, TRUE, FALSE))
}  )

# Identifying outliers in the data with the robust method
which(Rob.Out.detect(data))

# Removing outlier
datao <- data[!Rob.Out.detect(data)]

length(datao)

hist(datao, breaks = 15, xlab = "BMI", ylab = "Density", probability = TRUE,
     cex.lab = 1.5, cex.axis = 1.5, main = "Removing outlier") 
box()

boxplot(data, ylab = "BMI",
     cex.lab = 1.5, cex.axis = 1.5, main = "Removing outlier")


## ------------------------------------------------------------------------------------------------
# Normality tests
NormTestO <- snNormalityTest(datao, method = "ILA")

# Posterior probability of sn model based on Jeffreys prior
Oprobsn_J <- NormTestO$probsn_J
print(Oprobsn_J)

# Posterior probability of sn model based on DIMOM prior
Oprobsn_MOM <- NormTestO$probsn_MOM
print(Oprobsn_MOM)

# Posterior probability of sn model based on MOOMIN prior
Oprobsn_MOOMIN <- NormTestO$probsn_MOOMIN
print(Oprobsn_MOOMIN)

# BIC sn model
BICSNO <- NormTestO$BICSN
print(BICSNO)

# BIC normal model
BICNO <- NormTestO$BICN
print(BICNO)

# p-value
pswo <- NormTestO$pval
print(pswo)


# Estimators

ESTO <- rbind(c(NormTestO$MLEN,NA), NormTestO$MLESN,
             c(NormTestO$MAPN,NA),
             NormTestO$MAPSN_J, 
             NormTestO$MAPSN_MOM, NormTestO$MAPSN_MOOMIN)

rownames(ESTO) <- c("MLE_N", "MLE_SN",
                  "MAP_N", 
                  "MAPSN_J",
                  "MAPSN_DIMOM",  "MAPSN_MOOMIN")


kable(ESTO, digits = 3, format = "html") %>%
  kable_styling(full_width = FALSE) %>%
  column_spec(1:ncol(EST), extra_css = "padding-left: 15px; padding-right: 15px;")


