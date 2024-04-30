#!/usr/bin/Rscript

# library(ggplot2) # not used?

# acf and pacf ----
plotAcf = function(x, n, lag.max = NULL, pacf=F){
  # Function to plot acf and pacf
  
  if (is.null(lag.max)){
    lag.max = as.integer(length(x)-1)
  }
  lag.max = min(lag.max, length(x)-1) 
  
  # lags = 0:lag.max
  CI <- c(1.96/sqrt(n), -1.96/sqrt(n))
  
  if (pacf) {
    ylim = c(min(CI[2],x),max(CI[1], x))
    ylab = "Partial ACF"
    lags = 1:(lag.max+1)
  }
  else{
    ylim = c(min(CI[2],x),1)
    ylab = "ACF"
    lags = 0:lag.max
  }
  
  plot(lags, x, type="h", ylab = ylab, xlab = "Lag", 
       ylim = ylim) 
  abline(h=0) + abline(h=CI, col=c("blue","blue"), lty=c(2,2))
}

# One lag ac
sample_acf <- function(h,x){
  #Function to calculate the sample autocorrelation function.
  #h: lag to evaluate ACF at
  #x: time series data
  
  x_mean = mean(x)
  t <- length(x)
  sample_acvf <- 0
  for (i in 1:(t-h)){
    sample_acvf <- sample_acvf + (x[i+h]-x_mean)*(x[i]-x_mean)
  }
  sample_acvf <- sample_acvf * (1/t)
  
  if (h == 0){
    return (sample_acvf)
  }
  else{
    return (sample_acvf / sample_acf(0,x))
  }
}

# acf
sampleAcf = function(x, lag.max = NULL, plot=T){
  # Computes sample acf for all lags until lag.max
  # x: Time series
  
  # Assert lag.max not too large
  if (is.null(lag.max)){
    lag.max = as.integer(length(x)/4)
  }
  lag.max = min(lag.max, length(x)-1) 
  
  lags <- 0:lag.max
  result <- sapply(lags,sample_acf,x = x)
  result[1] <- result[1] / sample_acf(0,x)
  if (plot) {
    plotAcf(result, length(x), lag.max)
  }
  
  return(result)
}

# Durbin Levinson
durbin_levinson <- function(h,x){
  #Simplified version of Durbin-Levinson.
  phi <- matrix(nrow = h+1,ncol = h+1)
  phi[1,1] <- 0
  phi[2,2] <- sample_acf(1,x)
  
  for (n in 2:(h)){
    numer <- sample_acf(n,x)
    denom <- 1
    for (k in 1:(n-1)){
      if ((n-1) != k){
        phi[n,k+1] <- phi[n-1,k+1] - phi[n,n]*phi[n-1,n-k]
        phi[k+1,n] <- phi[n,k+1]
      }
      
      numer <- numer - phi[n,k+1]*sample_acf(n-k,x)
      denom <- denom - phi[n,k+1]*sample_acf(k,x)
    }
    phi[n+1,n+1] <- as.numeric(numer/denom)
  }
  return (diag(phi)[-1])
}

# pacf
samplePacf <- function(x, lag.max=NULL, plot=T){
  #Computes and plots the partial autocorrelation function with a 95% confidence interval.
  #h: lag to compute PACF up to.
  #x: time series
  
  # Assert lag.max not too large
  if (is.null(lag.max)){
    lag.max = as.integer(length(x)/4)
  }
  lag.max = min(lag.max, length(x)-1) 
  
  coeffs <- durbin_levinson(lag.max,x)
  
  if (plot) {
    plotAcf(coeffs, length(x), lag.max, pacf = T)
  }
  
  return(coeffs)
}


# GARCH ----
# log likelihood (objective) function
critFunc = function(params, r, p=1, q=1){
  # LogLike of alpha and beta given the first r_t values t=(1,..., max(p,q))(Eq. 5.46)
  # Defaults to GARCH(1,1) and uses only the first three params (alpha0, alpha1, beta1)
  # params: vector of log(parameters). First p+1 are alphas, rest are betas
  # r: vector of returns (Eq. 5.34)
  # p,q: GARCH(p,q) model parameters
  
  # Don't allow for parameters to be smaller than zero by raising (optim might wanna do this)
  n = length(r)
  m = max(p,q)
  # If else statements in case it is not "generalized" (ARCH(p))
  if (p == 0){alpha = c(0); p=1}
  else {alpha = exp(params[1:(p+1)])}
  
  if (q == 0){beta = 0; q = 1}
  else {beta = exp(params[-(1:(p+1))])}
  
  if (min(c(alpha, beta))<0) {
    print("Some params smaller than zero")
    print(params)
  }
  
  # Initialize first m conditional variances.
  # TODO: This is not correct initialization of the first sigma[1:m], but it works for now.
  sigma_sq = numeric(n)
  sigma_sq[1:m] = t(alpha) %*% c(1,r[p:1]^2)
  # sigma[1:m] = alpha[1]/(1-sum(alpha[-1] + beta))
  
  ll = 0 # log likelihood (objective is to minimize this)
  for (t in (m+1):n){
    sigma_sq[t] = t(alpha) %*% c(1,r[(t-1):(t-p)]^2) + t(beta) %*% sigma_sq[(t-1):(t-q)]
    # ll = ll + log(sigma[t]) + r[t]^2/sigma[t] # can compute as below
  }
  ll = sum(log(sigma_sq) + r^2/sigma_sq)
  return(ll)
}

# One-step-ahead forecast of sigma
sigmaForecast = function(alpha, beta, r){
  # One-step-ahead forecasts of the volatility sigma (Eq. 5.52)
  # alpha = (alpha_0,..., alpha_p): Estimated parameters by num. optim.
  # beta = (beta_1,..., beta_q): Estimated parameters by num. optim.
  # r: vector of returns.
  
  if (is.null(beta)) {beta=0; q = 1}
  p = length(alpha) - 1
  q = length(beta)
  m = max(p,q)
  n = length(r)
  
  # Initialize first m conditional variances.
  # TODO: This is not correct initialization of the first sigma[1:m], but it works for now.
  sigma = numeric(n)
  sigma[1:m] = rep(t(alpha) %*% c(1,r[p:1]^2),m)
  for (t in (m+1):n){
    # sigma[t] = t(alpha) %*% c(1,r[(t-1):(t-p)]^2) + t(beta) %*% sigma[(t-1):(t-q)]
    sigma[t] = t(alpha) %*% c(1,r[(t-1):(t-p)]^2) + t(beta) %*% sigma[(t-1):(t-q)]^2
  }
  return(sigma)
}

# GARCH(p,q)
garch = function(p, q, alpha, beta, r){
  if (is.null(beta)){params=alpha}
  else {params = c(alpha, beta)}
  
  params = optim(par=params, fn = critFunc, r=r, p=p, q=q, method = "BFGS")$par
  sigma = sigmaForecast(alpha, beta, r)
  
  return(list(params = params, sigma = sigma))
}


# Convenience functions ----
extractAlphaBeta = function(params, p, q=0){
  return(list(alpha = params[1:(p+1)], beta = params[-(1:(p+1))]))
}

# Debug options ----
if (F){
  options(error=recover)
  options(error=NULL)
}

# runs only when script is run by itself
# Test ----
if (sys.nframe() != 0){
  # Load and find well behaved TS
  df <- read.csv('projectdata.csv', header = T, sep=";", dec=",", stringsAsFactors=FALSE)
  infl = df$Inflation
  fitVals = lm(Inflation~Unemployed+Consumption+InterestRate, data = df)$fitted.values
  res = infl - fitVals
  # par(mfrow = c(1,3))
  # qqnorm(infl)
  # qqnorm(diff(infl))
  # qqnorm(res)
  # par(mfrow=c(1,1))

  # Fit GARCH
  # garch10 = garch(1,0,alpha = c(0.1,0.1), beta=NULL, r = res)
  params = exp(optim(par=log(c(0.1,0.1)), fn = critFunc, r = res, p=1, q=0, method="BFGS")$par)
  print(c(params = params))
  sigma = sigmaForecast(params, NULL, res)
  plot(res)
  lines(sigma)
}


# Test: Manufactured data ----
## ARCH(1) ----
ar1_a = c(0.2,0.4)
ar1_p = length(ar1_a)-1
ar1_n = 1e3
ar1_rs = numeric(ar1_n) # r sample
for (t in (ar1_p+1):ar1_n){
  rs[t] = rnorm(1, sd=sqrt(t(ar1_a)%*%c(1,ar1_rs[(t-1):(t-ar1_p)]^2)))
}
ar1_init = c(0.1,0.1)
ar1_params = exp(optim(par=log(ar1_init), fn=critFunc, r=ar1_rs, p=ar1_p, q=0)$par)
rbind(real = ar1_a, estim = ar1_params)

## GARCH(2,2) ----
garch22 = list(
  p = 2, 
  q= 2,
  m= 2,
  n = 1e4,
  alpha_init = rep(0.1,3),
  beta_init = rep(0.1,2),
  # Simulation values
  alpha_sim = 
)




# Initialize vectors to store simulated data
gar22_rs <- numeric(garch22$n)

# Set initial values
gar22_rs[1:2] <- rnorm(2)
# conditional_variances[1:2] <- omega / (1 - alpha1 - alpha2 - beta1 - beta2)  # Ensure that it's a valid GARCH process

# Simulate GARCH(2,2) process
for (t in (garch22$m+1):garch22$n) {
  gar22_rs[t] = rnorm(1, sd = sqrt(t(garch22$alpha_sim)%*%c(1,gar22_rs[(t-1):(t-garch22$p)]^2) + 
                                     t(garch22$beta_sim) %*% gar22_rs[(t-1):(t-garch22$q)]^2))
  if (is.nan(gar22_rs[t])){
    browser()
  }
}

qqnorm(gar22_rs)



params = exp(optim(par=log(c(garch22$alpha_init, garch22$beta_init)), 
                   fn=critFunc, r=gar22_rs, p=garch22$p, q=garch22$q, method = "BFGS")$par)
round(rbind(real = c(garch22$alpha_sim,garch22$beta_sim), 
            estim = params, 
            init = c(garch22$alpha_init, garch22$beta_init)),2)

ae = extractAlphaBeta(params, 2)$alpha
be = extractAlphaBeta(params, 2)$beta

sigma = sigmaForecast(ae, be, rs)
plot(rs, type="l")
lines(sigma, col="cyan", lty=2)


