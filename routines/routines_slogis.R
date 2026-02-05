
################################################################################
# Initial point for skew normal
################################################################################
# x: data set
sn_start <- function(x) {
  # sample moments
  m <- mean(x)
  s <- sd(x)
  g1 <- mean(((x - m) / s)^3)  # sample skewness
  
  # constant
  c <- (4 - pi) / 2
  
  if (abs(g1) < 1e-6) {
    # essentially symmetric → normal
    lambda <- 0
    delta <- 0
  } else {
    # method-of-moments inversion
    k <- (abs(g1) / c)^(2/3)
    omega <- sign(g1) * sqrt(k / (1 + k))
    delta <- omega * sqrt(pi / 2)
    # numerical safeguard
    delta <- max(min(delta, 0.995), -0.995)
    lambda <- delta / sqrt(1 - delta^2)
  }
  
  sigma <- s / sqrt(1 - (2 * delta^2) / pi)
  mu <- m - sigma * delta * sqrt(2 / pi)
  
  out =  c(mu = mu, sigma = sigma, lambda = lambda)
  
  return(out)
}


################################################################################
# Minimum Discrepancy measure: normal(mu,sigma) vs skew logistic
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
# Signed Minimum Discrepancy measure: normal(mu,sigma) vs skew normal
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
    my_dp <- c(xi = 0, omega = 1, alpha = lambda)
    my_cp <- dp2cp(my_dp, family = "SN")
    init <- c(my_cp[1],log(my_cp[2]))
    
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
    my_dp <- c(xi = 0, omega = 1, alpha = lambda)
    my_cp <- dp2cp(my_dp, family = "SN")
    init <- c(my_cp[1],log(my_cp[2]))
    
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






