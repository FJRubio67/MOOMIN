
################################################################################
# Minimum Discrepancy measure: exp(lambda) vs weibull(lambda,nu)
# nu: real number, log(shape)
################################################################################
# lambda: real parameter

discrepancy_min_w <- Vectorize(function(nu){
  # Absolute value using the symmetry of the discrepancy measure
  if(exp(nu) == 1 ) val = 0
  if(exp(nu) != 1){
    # Discrepancy
    disc <- Vectorize(function(par){
      # Integrand
      tempf <- Vectorize(function(x){
        num <- dexp(x, rate = exp(par))^2
        den <- dexp(x, rate = exp(par)) + dweibull(x, scale = 1, shape = exp(nu))
        out <- num/den
        return(out)
      })
      # Integral (heuristic choice of integration range)
      int <- integrate(tempf,0,5)$value
      return(int)
    })
    
    # Initial value (mean and sd)
    init <- c(0.5)
    
    # Minimum translated discrepancy
    val <- optimize(f = disc, interval = c(-10,10))$objective-0.5
  }
  return(val)
})

################################################################################
# Signed Minimum Discrepancy measure: sech(mu,sigma) vs skew hyperbolic secant
# lambda: real number, skewness parameter
################################################################################
# lambda: real parameter

sdiscrepancy_min_shs <- Vectorize(function(lambda){
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
        num <- dsech(x,par[1],exp(par[2]))^2
        den <- dsech(x,par[1],exp(par[2])) + 2*dsech(x)*psech(lambda*x)
        out <- num/den
        return(out)
      })
      # Integral
      int <- integrate(tempf,-5,5)$value
      return(int)
    }
    
    # Initial value (mean and sd)
    init <- c(0.5,0)
    
    # Minimum translated signed discrepancy
    val <- optim(init,disc)$value-0.5
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
unprior_min_shs <- Vectorize(function(lambda){
  if(lambda ==0 ) prior = 0
  if(lambda != 0){
    # Absolute value using the symmetry of the discrepancy measure
    lambda <- abs(lambda)
    # Discrepancy
    disc <- function(par){
      # Integrand of discrepancy
      tempf <- Vectorize(function(x){
        num <- dsech(x,par[1],exp(par[2]))^2
        den <- dsech(x,par[1],exp(par[2])) + 2*dsech(x)*psech(lambda*x)
        out <- num/den
        return(out)
      })
      # Integral
      int <- integrate(tempf,-5,5)$value
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
      num <- -(dsech(x,mu_opt,sd_opt)^2)*2*dsech(x)*dsech(lambda*x)*x
      den <- (dsech(x,mu_opt,sd_opt) + 2*dsech(x)*psech(lambda*x))^2
      out <- num/den
      return(out)
    })
    # prior
    prior <- abs(integrate(tempf2,-5,5)$value)
  }
  return(prior)
  
  
})






