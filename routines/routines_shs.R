
#------------------------------------------------------------------
# PDF
#------------------------------------------------------------------

dsech <- Vectorize(function(x,mu=0,sigma=1,log = FALSE){
  logden <-  -log(2) - log(sigma) - log( cosh( 0.5*pi*(x-mu)/sigma ) ) 
  val <- ifelse(log, logden, exp(logden)) 
  return(val)
})

#------------------------------------------------------------------
# CDF
#------------------------------------------------------------------

psech <- Vectorize(function(x,mu=0,sigma=1,log.p = FALSE){
  logcdf <-  log(2) - log(pi) + log( atan( exp( 0.5*pi*(x-mu)/sigma ) ) )
  val <- ifelse(log.p, logcdf, exp(logcdf))
  return(val)
})

#------------------------------------------------------------------
# Quantile function
#------------------------------------------------------------------

qsech <- Vectorize(function(p,mu=0,sigma=1){
  val <- sigma*2*log( tan( 0.5*pi*p ) )/pi + mu
  return(val)
})

#------------------------------------------------------------------
# Random number generation
#------------------------------------------------------------------

rsech <- function(n,mu=0,sigma=1){
  u <- runif(n)
  val <-  sigma*2*log( tan( 0.5*pi*u ) )/pi + mu
  return(val)
}


################################################################################
# Minimum Discrepancy measure: sech(mu,sigma) vs skew hyperbolic secant
# lambda: real number, skewness parameter
################################################################################
# lambda: real parameter

discrepancy_min_shs <- Vectorize(function(lambda){
  # Absolute value using the symmetry of the discrepancy measure
  if(lambda ==0 ) val = 0
  if(lambda != 0){
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
      # Integral (heuristic choice of integration range)
      int <- integrate(tempf,-5,5)$value
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






