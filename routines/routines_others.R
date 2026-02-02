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
unprior_min <- Vectorize(function(lambda){
  # Absolute value using the symmetry of the discrepancy measure
  lambda <- abs(lambda)
  # Discrepancy
  disc <- function(par){
    # Integrand of discrepancy
    tempf <- Vectorize(function(x){
      num <- dnorm(x,par[1],exp(par[2]))^2
      den <- dnorm(x,par[1],exp(par[2])) + 2*dnorm(x)*pnorm(lambda*x)
      out <- num/den
      return(out)
    })
    # Integral
    int <- integrate(tempf,-5,7.5)$value
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
    num <- -(dnorm(x,mu_opt,sd_opt)^2)*2*dnorm(x)*dnorm(lambda*x)*x
    den <- (dnorm(x,mu_opt,sd_opt) + 2*dnorm(x)*pnorm(lambda*x))^2
    out <- num/den
    return(out)
  })
  # prior
  prior <- abs(integrate(tempf2,-5,7.5)$value)
  return(prior)
  
  
})

# Normalising constant
#nc_prior <- integrate(unprior_min,-100,100)$value
#nc_prior <- 0.07327272
nc_prior <- 0.08034964

# Normalised prior based on minimum discrepancy (based on normalising prior nc_prior)
# Based on numerical integration of the expression obtained in a Proposition
prior_min <- Vectorize(function(lambda){ unprior_min(lambda)/nc_prior })

# Approximated normalised prior based on minimum discrepancy
alpha_moomin = 2
k_moomin = 4
a_moomin <- 0.28
# alpha_moomin = 2
# k_moomin = 4
# a_moomin <- 0.275
m_moomin <- (k_moomin + alpha_moomin )/2
nc_moomin <- (a_moomin ^(-0.5 - k_moomin /2.)*
                gamma((1 + k_moomin )/2.)*
                gamma(-0.5 - k_moomin /2. + m_moomin ))/gamma(m_moomin )


tprior_app_min <- Vectorize(function(lambda){
  
  out <- abs(lambda)^k_moomin /(nc_moomin *(1+a_moomin *lambda^2)^m_moomin )
  
  return(out)
  
})


ltprior_app_min <- Vectorize(function(lambda){
  
  out <- k_moomin*log(abs(lambda)) -
    log(nc_moomin) - m_moomin*log(1+a_moomin *lambda^2)
  
  return(out)
  
})
