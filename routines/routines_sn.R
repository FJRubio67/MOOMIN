################################################################################
# t-approximation to the Jeffreys prior
################################################################################
cJ <-pi/2 # scale parameter

# Density
da_jeff <- Vectorize(function(lambda) dt(lambda/cJ,df=1/2)/cJ)

# log-density
lda_jeff <- Vectorize(function(lambda){
 out <-   dt(lambda/cJ,df=1/2, log = TRUE) -log(cJ) 
})


################################################################################
# t-approximation to the TV prior
################################################################################
c_TV <- 0.92 # scale parameter

# Density
da_tv <- Vectorize(function(lambda) dt(lambda/c_TV,df=1)/c_TV)

# log-density
lda_tv <- Vectorize(function(lambda){
  out <-   dt(lambda/c_TV,df=1, log = TRUE) -log(c_TV) 
})

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
# DIMOM prior
################################################################################

s_mom = 1.69 # scale parameter (5% probability mass on |lambda| < 1)

# Density
d_mom <- Vectorize(function(lambda) (lambda^2/s_mom^2)*dnorm(lambda,mean = 0,sd = s_mom) )

# Log-density
ld_mom <- Vectorize(function(lambda) 2*log(abs(lambda)) -
                      2*log(s_mom) + dnorm(lambda,mean = 0,sd = s_mom, log = TRUE) )



################################################################################
# Minimum Discrepancy measure: normal(mu,sigma) vs skew normal
# lambda: real number, skewness parameter
################################################################################
# lambda: real parameter

discrepancy_min <- Vectorize(function(lambda){
  # Absolute value using the symmetry of the discrepancy measure
  lambda <- abs(lambda)
  # Discrepancy
  disc <- function(par){
    # Integrand
    tempf <- Vectorize(function(x){
      num <- dnorm(x,par[1],exp(par[2]))^2
      den <- dnorm(x,par[1],exp(par[2])) + 2*dnorm(x)*pnorm(lambda*x)
      out <- num/den
      return(out)
    })
    # Integral (heuristic choice of integration range)
    int <- integrate(tempf,-5,7.5)$value
    return(int)
  }
  
  # Initial value (mean and sd)
  my_dp <- c(xi = 0, omega = 1, alpha = lambda)
  my_cp <- dp2cp(my_dp, family = "SN")
  init <- c(my_cp[1],log(my_cp[2]))
  
  # Minimum translated discrepancy
  val <- optim(init,disc)$value-0.5
  return(val)
})

################################################################################
# Signed Minimum Discrepancy measure: normal(mu,sigma) vs skew normal
# lambda: real number, skewness parameter
################################################################################
# lambda: real parameter

sdiscrepancy_min <- Vectorize(function(lambda){
  # Absolute value using the symmetry of the discrepancy measure
  sgn <- sign(lambda)
  lambda <- abs(lambda)
  # Discrepancy
  disc <- function(par){
    # Integrand
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


################################################################################
# Marginal likelihood for normal model under reference prior
# Direct integration
################################################################################
# data: a numerical vector containing the data set

logmlN <-  function(data, log = FALSE) {
  if (!is.numeric(data)) {
    stop("data must be a numeric vector")
  }
  
  n <- length(data)
  if (n < 2) {
    stop("Need at least n >= 2 observations")
  }
  
  xbar <- mean(data)
  S <- sum((data - xbar)^2)
  
  if (S <= 0) {
    stop("Sum of squares is zero; marginal likelihood is not finite")
  }
  
  log_m <- 
    -((n - 1) / 2) * log(2 * pi) -
    0.5 * log(n) +
    lgamma((n - 1) / 2) +
    ((n - 1) / 2) * (log(2) - log(S))
  
  if (log) {
    return(log_m)
  } else {
    return(exp(log_m))
  }
}

################################################################################
# Marginal likelihood for normal model under reference prior
# Integrated Laplace approximation
################################################################################
# data: a numerical vector containing the data set
# log_pi_lambda: log-prior on lambda

laplace_mu_sigma <- function(data, lambda) {
  
  start <- c(mean(data), log(sd(data)))
  
  neg_logpost <- function(par) {
    mu = par[1]
    sigma = exp(par[2])
    loglik <-  sum(dsn(data, xi = mu, omega = sigma, alpha = lambda, log = TRUE))
    lprior <- - log(sigma)
    ljac <- par[2]
    out <- -loglik - lprior - ljac
    
    return(out)
  }
  
  opt <- optim(
    start,
    neg_logpost,
    hessian = TRUE,
    method = "BFGS"
  )
  
  if (opt$convergence != 0)
    stop("Optimization failed for lambda = ", lambda)
  
  log_det_hess <- determinant(opt$hessian, logarithm = TRUE)$modulus
  
  log_m_lambda <-  0.5*2*log(2*pi) - opt$value - 
    0.5*log_det_hess
  
  
  
  return(as.numeric(log_m_lambda))
}


ml_sn_integrated <- function(data, log_pi_lambda,
                             lambda_range = c(-10, 10)) {
  
  integrand <- Vectorize(function(lambda) {
    out <- laplace_mu_sigma(data, lambda) + log_pi_lambda(lambda)
    return(out)
    
  })  
  
  out <- integrate(
    function(l) exp(integrand(l)),
    lower = lambda_range[1],
    upper = lambda_range[2],
    rel.tol = 1e-8
  )
  
 ml <- log(out$value)
 
 return(ml)
}

################################################################################
# Normality tests against skew-normal alternatives 
# It also calculates the Shapiro-Wilks test
################################################################################
# data: a numerical vector containing the data set
# method : Integrated Laplace Approximation (ILA) or Laplace Approximation (LA)

snNormalityTest <- function(data, method = "ILA"){
  
  data <- as.vector(data)
  
  n <- length(data)
  sd_pop <- sd(data) * sqrt((n - 1) / n)
  m_pop <- mean(data)
  signl <- sign(mean(data)-median(data))
  initsn = sn_start(data)
  
  #---------------------------------------
  # Calculations for normal model
  # Closed-form calculations
  #---------------------------------------
  
  # log-likelihood
  loglikn <- function(par){
    
    mu = par[1]
    sigma = exp(par[2])
    
    loglik <-  sum(dnorm(data, mean = mu, sd = sigma, log = TRUE))
    
    out <- -loglik 
    
    return(out)
    
  }
  
  # log-posterior
  logpostn <- function(par){
    
    mu = par[1]
    sigma = exp(par[2])
    
    loglik <-  sum(dnorm(data, mean = mu, sd = sigma, log = TRUE))
    lprior <- - log(sigma)
    ljac <- par[2]
    
    out <- -loglik - lprior - ljac
    
    return(out)
    
  }
  
  # MLE
  MLEN = c(m_pop, log(sd_pop))
  # BIC
  BICN = 2*loglikn(MLEN) + log(n)*length(MLEN)
  
  # MAP and marginal likelihood
  OPTPN = nlminb(start = MLEN, logpostn, control = list(iter.max = 10000))
  MAPN <- OPTPN$par

  mln <- logmlN(data, log = TRUE)
  mln <- as.numeric(mln)
  
  #-------------------------------------------------------------
  # Calculations for skew-normal model: MLE
  #-------------------------------------------------------------
  
  # log-likelihood
  logliksn <- function(par){
    
    mu = par[1]
    sigma = exp(par[2])
    lambda = par[3]
    
    loglik <-  sum(dsn(data, xi = mu, omega = sigma, alpha = lambda, log = TRUE))
    
    out <- -loglik
    
    return(out)
    
  }
  
  # MLE
  OPTLSN <- nlminb(start = initsn, logliksn, control = list(iter.max = 10000))
  MLESN <- OPTLSN$par
  # BIC
  BICSN <- 2*OPTLSN$objective + log(n)*length(OPTLSN$par)
  
  
  #-------------------------------------------------------------
  # Calculations for skew-normal model: Jeffrey's prior
  #-------------------------------------------------------------
  
  
  # log-posterior
  logpostsnJ <- function(par){
    
    mu = par[1]
    sigma = exp(par[2])
    lambda = par[3]
    
    loglik <-  sum(dsn(data, xi = mu, omega = sigma, alpha = lambda, log = TRUE))
    lprior <- lda_jeff(lambda) - log(sigma)
    ljac <- par[2]
    
    out <- -loglik - lprior - ljac
    
    return(out)
    
  }
  
  # MAP and marginal likelihood
  OPTPSNJ = nlminb(start = initsn, logpostsnJ, control = list(iter.max = 10000))
  MAPSNJ <- OPTPSNJ$par
  
  if(method == "LA"){
  
  LPSNJ <- -OPTPSNJ$objective
  HessSNJ <- hessian(logpostsnJ, x = MAPSNJ)
  npsnJ <- length(MAPSNJ)
  
  mlsnJ <- 0.5*npsnJ*log(2*pi) + LPSNJ - 
    0.5*as.numeric(determinant(HessSNJ, logarithm = TRUE)$modulus)
  mlsnJ <- as.numeric(mlsnJ)
  }
  
  if(method == "ILA"){
    
    mlsnJ <- as.numeric(ml_sn_integrated(data, lda_jeff))
  }
  
  #-------------------------------------------------------------
  # Calculations for skew-normal model: TV prior
  #-------------------------------------------------------------
  
  
  # log-posterior
  logpostsnTV <- function(par){
    
    mu = par[1]
    sigma = exp(par[2])
    lambda = par[3]
    
    loglik <-  sum(dsn(data, xi = mu, omega = sigma, alpha = lambda, log = TRUE))
    lprior <- lda_tv(lambda) - log(sigma)
    ljac <- par[2]
    
    out <- -loglik - lprior - ljac
    
    return(out)
    
  }
  
  # MAP and marginal likelihood
  OPTPSNTV = nlminb(start = initsn, logpostsnTV, control = list(iter.max = 10000))
  MAPSNTV <- OPTPSNTV$par
  
  if(method == "LA"){
  LPSNTV <- -OPTPSNTV$objective
  HessSNTV <- hessian(logpostsnTV, x = MAPSNTV)
  npsnTV <- length(MAPSNTV)
  
  mlsnTV <- 0.5*npsnTV*log(2*pi) + LPSNTV - 
    0.5*as.numeric(determinant(HessSNTV, logarithm = TRUE)$modulus)
  mlsnTV <- as.numeric(mlsnTV)
  }
  
  if(method == "ILA"){
    
    mlsnTV <-  as.numeric(ml_sn_integrated(data, lda_tv))
    
  }
  
  
  #-------------------------------------------------------------
  # Calculations for skew-normal model: MOM prior
  #-------------------------------------------------------------
  
  
  # log-posterior
  logpostsnM <- function(par){
    
    mu = par[1]
    sigma = exp(par[2])
    lambda = par[3]
    
    loglik <-  sum(dsn(data, xi = mu, omega = sigma, alpha = lambda, log = TRUE))
    lprior <- ld_mom(lambda) - log(sigma)
    ljac <- par[2]
    
    out <- -loglik - lprior - ljac
    
    return(out)
    
  }
  
  
  # MAP and marginal likelihood
  OPTPSNM = nlminb(start = initsn, logpostsnM, control = list(iter.max = 10000))
  MAPSNM <- OPTPSNM$par
  
  
  if(method == "LA"){
  
  LPSNM <- -OPTPSNM$objective
  HessSNM <- hessian(logpostsnM, x = MAPSNM)
  npsnM <- length(MAPSNM)
  
  mlsnM <- 0.5*npsnM*log(2*pi) + LPSNM - 
    0.5*as.numeric(determinant(HessSNM, logarithm = TRUE)$modulus)
  mlsnM <- as.numeric(mlsnM)
  }
  
  if(method =="ILA"){
    
    mlsnM <- as.numeric(ml_sn_integrated(data, ld_mom))
  }
  
  #-------------------------------------------------------------
  # Calculations for skew-normal model: MOOMIN prior
  # Based on the approximate MOOMIN prior
  #-------------------------------------------------------------
  
  # log-posterior
  logpostsn_min <- function(par){
    
    mu = par[1]
    sigma = exp(par[2])
    lambda = par[3]
    
    loglik <-  sum(dsn(data, xi = mu, omega = sigma, alpha = lambda, log = TRUE))
    lprior <- ltprior_app_min(lambda) - log(sigma)
    ljac <- par[2]
    
    out <- -loglik - lprior - ljac
    
    return(out)
    
  }
  
  # MAP and marginal likelihood
  OPTPSN_min = nlminb(start = initsn, logpostsn_min, control = list(iter.max = 10000))
  MAPSN_min <- OPTPSN_min$par
  
  if(method == "LA"){
  LPSN_min <- -OPTPSN_min$objective
  HessSN_min <- hessian(logpostsn_min, x = MAPSN_min)
  npsn_min <- length(MAPSN_min)
  
  mlsn_min <- 0.5*npsn_min*log(2*pi) + LPSN_min - 
    0.5*as.numeric(determinant(HessSN_min, logarithm = TRUE)$modulus)
  mlsn_min <- as.numeric(mlsn_min)
  }
  
  if(method == "ILA"){
    mlsn_min <- as.numeric(ml_sn_integrated(data, ltprior_app_min))
    
  }
  
  #-------------------------------------------------------------
  # Processing output
  #-------------------------------------------------------------

  BF10J <- exp(mlsnJ - mln)
  
  BF10TV <- exp(mlsnTV - mln)
  
  BF10M <- exp(mlsnM - mln)
  
  BF10_min <- exp(mlsn_min - mln)
  
  psnJ <- 1/(1 + exp(-mlsnJ + mln))
  
  psnTV <- 1/(1 + exp(-mlsnTV + mln))
  
  psnM <- 1/(1 + exp(-mlsnM + mln))
  
  psn_min <- 1/(1 + exp(-mlsn_min + mln))

  
  if(n<=5000){
    pval = shapiro.test(data)$p.value
    
    output <- list(MLEN = MLEN, MLESN = MLESN,
      MAPN = MAPN,  MAPSN_J = MAPSNJ, MAPSN_TV = MAPSNTV, MAPSN_MOM = MAPSNM, MAPSN_MOOMIN = MAPSN_min, 
      lmarglikn = mln, lmargliksn_J = mlsnJ, lmargliksn_TV = mlsnTV, lmargliksn_MOM = mlsnM, lmargliksn_MOOMIN = mlsn_min,
      BF10_J = BF10J, BF10_TV = BF10TV, BF10_MOM = BF10M, BF10_MOOMIN = BF10_min,
      probsn_J = psnJ, probsn_TV = psnTV, probsn_MOM = psnM,  probsn_MOOMIN = psn_min,
      BICN = BICN, BICSN = BICSN,
      pval = pval)
  }
  if(n>5000){

    
    output <- list(MLEN = MLEN, MLESN = MLESN,
                   MAPN = MAPN,  MAPSN_J = MAPSNJ, MAPSN_TV = MAPSNTV, MAPSN_MOM = MAPSNM, MAPSN_MOOMIN = MAPSN_min, 
                   lmarglikn = mln, lmargliksn_J = mlsnJ, lmargliksn_TV = mlsnTV, lmargliksn_MOM = mlsnM, lmargliksn_MOOMIN = mlsn_min,
                   BF10_J = BF10J, BF10_TV = BF10TV, BF10_MOM = BF10M, BF10_MOOMIN = BF10_min,
                   probsn_J = psnJ, probsn_TV = psnTV, probsn_MOM = psnM,  probsn_MOOMIN = psn_min,
                   BICN = BICN, BICSN = BICSN)
  }
  
  
return(output)
    
}










