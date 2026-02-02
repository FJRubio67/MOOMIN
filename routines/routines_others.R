logit <- Vectorize(function(x) qlogis(x))
expit <- Vectorize(function(x) plogis(x))

################################################################################
# Signed Minimum Discrepancy measure: normal(mu,sigma) vs twopiece normal
# par: real number, skewness parameter
################################################################################
# lambda: real parameter

sdiscrepancy_min_tp <- Vectorize(function(eps){
  # Absolute value using the symmetry of the discrepancy measure
  sgn <- sign(eps)
  eps <- tanh(eps)
  eps = abs(eps)
  dtp <- Vectorize(function(x) dtp3(x,0,1,eps, param = "eps", FUN = dnorm))
  # Discrepancy
  disc <- function(par){
    # Integrand
    tempf <- Vectorize(function(x){
      num <- dnorm(x,par[1],exp(par[2]))^2
      den <- dnorm(x,par[1],exp(par[2])) + dtp(x)
      out <- num/den
      return(out)
    })
    # Integral
    int <- integrate(tempf,-15,15)$value
    return(int)
  }
  
  # Initial value (mean and sd)
  init <- c(0,0)
  
  # Minimum translated signed discrepancy
  val <- optim(init,disc)$value-0.5
  return(sgn*abs(val))
})


################################################################################
# MOOMIN Prior on the skewness parameter
# Based on assuming a uniform prior on the signed minimum discrepancy measure
################################################################################
# lambda: real parameter


# Un-normalised prior
# Based on numerical integration of the expression obtained in a Proposition
unprior_mintp <- Vectorize(function(eps){
  eps <- tanh(eps)
  eps = abs(eps)
  # prior
  prior <- grad(sdiscrepancy_min_tp, x = eps)
  return(prior)
})

curve(unprior_mintp,-2,2)





dsas <- function(x,mu,sigma,epsilon,delta,log=FALSE){
  ifelse(sigma>0 & delta>0,
         logPDF <- dnorm(sinh(delta*asinh((x-mu)/sigma)-epsilon),0,1,log=T) + log(delta) + log(cosh(delta*asinh((x-mu)/sigma)-epsilon)) -0.5*log(1+(x-mu)^2/sigma^2) - log(sigma),
         logPDF <- 'parameters out of range')
  ifelse( is.numeric(logPDF),ifelse( log, return(logPDF), return(exp(logPDF)) ), logPDF )
}

rsas <- function(n,mu,sigma,epsilon,delta){
  ifelse(sigma>0 & delta>0,
         sample <- mu+sigma*sinh((asinh(rnorm(n,0,1))+epsilon)/delta),
         sample <- 'parameters out of range')
  return(sample)
}

################################################################################
# Signed Minimum Discrepancy measure: normal(mu,sigma) vs twopiece normal
# par: real number, skewness parameter
################################################################################
# lambda: real parameter

sdiscrepancy_min_sas <- Vectorize(function(eps){
  # Absolute value using the symmetry of the discrepancy measure
  sgn <- sign(eps)
  eps = abs(eps)
  ds <- Vectorize(function(x) dsas(x,0,1,eps,1))
  # Discrepancy
  disc <- function(par){
    # Integrand
    tempf <- Vectorize(function(x){
      num <- dnorm(x,par[1],exp(par[2]))^2
      den <- dnorm(x,par[1],exp(par[2])) + ds(x)
      out <- num/den
      return(out)
    })
    # Integral
    if(abs(eps)<=1)  int <- integrate(tempf,-5,5)$value
    if(abs(eps)>1)  int <- integrate(tempf,-50,50)$value
    return(int)
  }
  
  # Initial value (mean and sd)
  sim <- rsas(1000,0,1,eps,1)
  init <- c(mean(sim),log(sd(sim)))
  
  # Minimum translated signed discrepancy
  val <- nlminb(init,disc)$objective-0.5
  return(sgn*abs(val))
})


################################################################################
# MOOMIN Prior on the skewness parameter
# Based on assuming a uniform prior on the signed minimum discrepancy measure
################################################################################
# lambda: real parameter


# Un-normalised prior
# Based on numerical integration of the expression obtained in a Proposition
unprior_min_sas <- Vectorize(function(eps){
  # prior
  prior <- grad(sdiscrepancy_min_sas, x = eps)
  return(prior)
  
  
})

curve(unprior_min_sas,-2,2)

