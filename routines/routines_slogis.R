
################################################################################
# Minimum Discrepancy measure: logistic(mu,sigma) vs skew logistic
# lambda: real number, skewness parameter
################################################################################
# lambda: real parameter

discrepancy_min_slogis <- Vectorize(function(lambda){
  # Absolute value using the symmetry of the discrepancy measure
  if(lambda ==0 ) val = 0
  if(lambda != 0){
    lambda <- abs(lambda)
    # Discrepancy
    disc <- function(par){
      # Integrand
      tempf <- Vectorize(function(x){
        num <- dlogis(x,par[1],exp(par[2]))^2
        den <- dlogis(x,par[1],exp(par[2])) + 2*dlogis(x)*plogis(lambda*x)
        out <- num/den
        return(out)
      })
      # Integral (heuristic choice of integration range)
      int <- integrate(tempf,-10,10)$value
      return(int)
    }
    
    # Initial value (mean and sd)
    init <- c(0.5,0)
    
    # Minimum translated discrepancy
    val <- optim(init,disc)$value-0.5
  }
  return(val)
})

################################################################################
# Signed Minimum Discrepancy measure: logistic(mu,sigma) vs skew logistic
# lambda: real number, skewness parameter
################################################################################
# lambda: real parameter

sdiscrepancy_min_slogis <- Vectorize(function(lambda){
  if(lambda ==0 ){ 
        sgn = 0
        val = 0
  }
  if(lambda != 0){
    # Absolute value using the symmetry of the discrepancy measure
    sgn <- sign(lambda)
    lambda <- abs(lambda)
    # Discrepancy
    disc <- function(par){
      # Integrand
      tempf <- Vectorize(function(x){
        num <- dlogis(x,par[1],exp(par[2]))^2
        den <- dlogis(x,par[1],exp(par[2])) + 2*dlogis(x)*plogis(lambda*x)
        out <- num/den
        return(out)
      })
      # Integral
      int <- integrate(tempf,-10,10)$value
      return(int)
    }
    
    # Initial value (mean and sd)
    init <- c(0.5,0)
    
    # Minimum translated signed discrepancy
    val <- nlminb(init,disc)$objective-0.5
  }
  return(sgn*abs(val))
})


################################################################################
# MOOMIN Prior on the skewness parameter
# Based on assuming a uniform prior on the signed minimum discrepancy measure
################################################################################
# lambda: real parameter


# Un-normalised prior
# Based on numerical integration of the expression obtained in a Proposition
unprior_min_slogis <- Vectorize(function(lambda){
  if(lambda ==0 ) prior = 0
  if(lambda != 0){
    # Absolute value using the symmetry of the discrepancy measure
    lambda <- abs(lambda)
    # Discrepancy
    disc <- function(par){
      # Integrand of discrepancy
      tempf <- Vectorize(function(x){
        num <- dlogis(x,par[1],exp(par[2]))^2
        den <- dlogis(x,par[1],exp(par[2])) + 2*dlogis(x)*plogis(lambda*x)
        out <- num/den
        return(out)
      })
      # Integral
      int <- integrate(tempf,-10,10)$value
      return(int)
    }
    
    # Initial value (mean and sd)
    init <- c(0.5,0)
    
    # Minimum discrepancy
    OPT <- optim(init,disc)
    
    mu_opt <- OPT$par[1]
    sd_opt <- exp(OPT$par[2])
    
    
    # Integrand of prior
    tempf2 <- Vectorize(function(x){
      num <- -(dlogis(x,mu_opt,sd_opt)^2)*2*dlogis(x)*dlogis(lambda*x)*x
      den <- (dlogis(x,mu_opt,sd_opt) + 2*dlogis(x)*plogis(lambda*x))^2
      out <- num/den
      return(out)
    })
    # prior
    prior <- abs(integrate(tempf2,-10,10)$value)
  }
  return(prior)
  
  
})






