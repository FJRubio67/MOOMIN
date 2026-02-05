
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
