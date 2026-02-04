################################################################################
# Minimum Discrepancy measure: normal(mu,sigma) vs twopiece normal
# eps: real number, skewness parameter
################################################################################

discrepancy_min_tpn <- Vectorize(function(eps){
  # Absolute value using the symmetry of the discrepancy measure
  eps <- tanh(eps)
  eps = abs(eps)
  dtp <- Vectorize(function(x) dtp3(x,0,1,eps, param = "eps", FUN = dnorm))
  # Discrepancy
  disc <- function(par){
    # Integrand
    tempf <- Vectorize(function(x){
      num <- dnorm(x,par[1],exp(par[2]))^2
      den <- dnorm(x,par[1],exp(par[2])) + dtp(x)
      if(den <= 1e-8) out <- 0
      if(den > 1e-8) out <- num/den
      return(out)
    })
    # Integral
    int <- integrate(tempf,-12.5,0)$value +  integrate(tempf,0,12.5)$value
    return(int)
  }
  
  # Initial value (mean and sd)
  sim <- rtp3(10000,0,1,eps, param = "eps", FUN = rnorm)
  init <- c(mean(sim),sd(sim))
  
  # Minimum translated signed discrepancy
  val <- nlminb(init,disc)$objective-0.5
  return(abs(val))
})


################################################################################
# MOOMIN Prior on the skewness parameter of the twopiece normal
# Based on assuming a uniform prior on the signed minimum discrepancy measure
################################################################################
# par: real parameter


# Un-normalised prior
# Based on numerical integration of the expression obtained in a Proposition
unprior_min_tpn <- Vectorize(function(par){
  if(par==0) prior = 0
  if(par!=0){
    par = abs(par)
  # prior
  prior <- grad(discrepancy_min_tpn, x = par,  method.args=list(eps=1e-8) )
  }
  return(abs(prior))
})




################################################################################
# Minimum Discrepancy measure: logistic(mu,sigma) vs twopiece logistic
# eps: real number, skewness parameter
################################################################################

discrepancy_min_tplogis <- Vectorize(function(eps){
  # Absolute value using the symmetry of the discrepancy measure
  eps <- tanh(eps)
  eps = abs(eps)
  dtp <- Vectorize(function(x) dtp3(x,0,1,eps, param = "eps", FUN = dlogis))
  # Discrepancy
  disc <- function(par){
    # Integrand
    tempf <- Vectorize(function(x){
      num <- dlogis(x,par[1],exp(par[2]))^2
      den <- dlogis(x,par[1],exp(par[2])) + dtp(x)
      if(den <= 1e-8) out <- 0
      if(den > 1e-8) out <- num/den
      return(out)
    })
    # Integral
    int <- integrate(tempf,-75,0)$value +  integrate(tempf,0,75)$value
    return(int)
  }
  
  # Initial value (mean and sd)
  sim <- rtp3(10000,0,1,eps, param = "eps", FUN = rlogis)
  init <- c(mean(sim),sd(sim))
  
  # Minimum translated signed discrepancy
  val <- optim(init,disc)$value-0.5
  return(abs(val))
})


################################################################################
# MOOMIN Prior on the skewness parameter of the twopiece logistic
# Based on assuming a uniform prior on the signed minimum discrepancy measure
################################################################################
# par: real parameter


# Un-normalised prior
# Based on numerical integration of the expression obtained in a Proposition
unprior_min_tplogis <- Vectorize(function(par){
  if(par==0) prior = 0
  if(par!=0){
    par = abs(par)
  # prior
  prior <- grad(discrepancy_min_tplogis, x = par,  method.args=list(eps=1e-8) )
  }
  return(abs(prior))
})







# Probability Density Function: location parameter mu and scale parameter sigma
dlap = function(x,mu=0,sigma=1,log=FALSE){
  log.den = - log(2) - log(sigma) - abs(x-mu)/sigma
  if(log) log.den
  else exp(log.den)
}

# Random number generation 
rlap = function(n,mu=0,sigma=1){
  u = runif(n)
  rgen = ifelse(u<0.5,mu + sigma*log(2*u),mu - sigma*log(2*(1-u)))
  return(rgen)
}

################################################################################
# Signed Minimum Discrepancy measure: Laplace(mu,sigma) vs twopiece Laplace
# eps: real number, skewness parameter
################################################################################

discrepancy_min_tplap <- Vectorize(function(eps){
  # Absolute value using the symmetry of the discrepancy measure
  eps <- tanh(eps)
  eps = abs(eps)
  dtp <- Vectorize(function(x) dtp3(x,0,1,eps, param = "eps", FUN = dlap))
  # Discrepancy
  disc <- function(par){
    # Integrand
    tempf <- Vectorize(function(x){
      num <- dlap(x,par[1],exp(par[2]))^2
      den <- dlap(x,par[1],exp(par[2])) + dtp(x)
      if(den <= 1e-12) out <- 0
      if(den > 1e-12) out <- num/den
      return(out)
    })
    # Integral
    int <- integrate(tempf,-60,0)$value +  integrate(tempf,0,60)$value
    return(int)
  }
  
  # Initial value (mean and sd)
  sim <- rtp3(10000,0,1,eps, param = "eps", FUN = rlap)
  init <- c(mean(sim),sd(sim))
  
  # Minimum translated signed discrepancy
  val <- nlminb(init,disc)$objective-0.5
  return(abs(val))
})


################################################################################
# MOOMIN Prior on the skewness parameter of the twopiece Laplace
# Based on assuming a uniform prior on the signed minimum discrepancy measure
################################################################################
# lambda: real parameter


# Un-normalised prior
# Based on numerical integration of the expression obtained in a Proposition
unprior_min_tplap <- Vectorize(function(par){
  if(par==0) prior = 0
  if(par!=0){
    par = abs(par)
  # prior
  prior <- grad(sdiscrepancy_min_tplap, x = par,  method.args=list(eps=1e-8) )
  }
  return(abs(prior))
})








