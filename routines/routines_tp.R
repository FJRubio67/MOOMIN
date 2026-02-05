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
      if(den <= 1e-12) out <- 0
      if(den > 1e-12) out <- num/den
      return(out)
    })
    # Integral
    int <- integrate(tempf,-12.5,0)$value +  integrate(tempf,0,12.5)$value
    return(int)
  }
  
  # Initial value (mean and sd)
  set.seed(123)
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
 # prior <- grad(discrepancy_min_tpn, x = par, 
#                method = "Richardson",
#                method.args = list(eps = 1e-6, d = 0.01, r = 6))
  
prior <- fderiv(discrepancy_min_tpn, x = par, h = 1e-12)
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
      out <- ifelse(den > 1e-15, num / den, 0)
      out[!is.finite(out)] <- 0
      return(out)
    })
    # Integral
#    int <- integrate(tempf,-75,0)$value +  integrate(tempf,0,75)$value
    # Integral with error handling
    int1 <- tryCatch(
      integrate(tempf, -75, 0)$value,
      error = function(e) return(0)
    )
    int2 <- tryCatch(
      integrate(tempf, 0, 75)$value,
      error = function(e) return(0)
    )
    int <- int1 + int2
    return(int)
  }


  # Initial value (mean and sd)
  set.seed(123)
  sim <- rtp3(10000,0,1,eps, param = "eps", FUN = rlap)
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
 # prior <- grad(discrepancy_min_tplogis, x = par, 
#                method = "Richardson",
#                method.args = list(eps = 1e-6, d = 0.01, r = 6))
    prior <- fderiv(discrepancy_min_tplogis, x = par, h = 1e-8)
    
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
      out <- ifelse(den > 1e-15, num / den, 0)
      out[!is.finite(out)] <- 0
      return(out)
    })
    # Integral
    # Integral with error handling
    int1 <- tryCatch(
      integrate(tempf, -60, 0)$value,
      error = function(e) return(0)
    )
    int2 <- tryCatch(
      integrate(tempf, 0, 60)$value,
      error = function(e) return(0)
    )
    int <- int1 + int2
    return(int)
  }
  
  # Initial value (mean and sd)
  set.seed(123)
  sim <- rtp3(10000,0,1,eps, param = "eps", FUN = rlap)
  init <- c(mean(sim),sd(sim))
  
  # Minimum translated signed discrepancy
  val <- optim(init,disc)$value-0.5
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
 # prior <- grad(discrepancy_min_tplap, x = par )
 prior <- fderiv(discrepancy_min_tplap, x = par, h = 1e-12)
  }
  return(abs(prior))
})





################################################################################
# Function to detect spikes
################################################################################

detect_spikes_robust <- function(x, y, k = 3) {
  # First derivative
  dy <- diff(y) / diff(x)
  
  # Second derivative (curvature)
  d2y <- diff(dy) / diff(x[-1])
  
  # Detect high curvature points
  threshold <- median(abs(d2y)) + k * mad(abs(d2y))
  spike_indices <- which(abs(d2y) > threshold) + 1
  
  # Filter out the zero point if needed
  spike_indices <- spike_indices[abs(x[spike_indices]) > 0.1]
  
  return(list(
    indices = spike_indices,
    x_values = x[spike_indices],
    y_values = y[spike_indices]
  ))
}


