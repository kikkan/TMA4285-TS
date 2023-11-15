#!/usr/bin/Rscript

library(ggplot2) # not used?

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
  # params: vector of parameters. First p+1 are alphas, rest are betas
  # r: vector of returns (Eq. 5.34)
  # p,q: GARCH(p,q) model parameters
  
  n = length(r)
  m = max(p,q)
  
  # If else statements in case it is not "generalized" (ARCH(p))
  if (p == 0){
    # Prolly not needed, never gonna remove the AR part?
    alpha = c(0)
    p=1
  }
  else {alpha = params[1:(p+1)]}
  if (q == 0){
    beta = c(0)
    q = 1
  }
  else {beta = params[(p+2):(1+p+q)]}
  
  # Initialize first m conditional variances.
  # TODO: This is not correct initialization of the first sigma[1:m], but it works for now.
  sigma = numeric(n)
  sigma[1:m] = (r[1:m]^2)
  
  ll = 0 # log likelihood (objective is to minimize this)
  for (t in (m+1):n){
    sigma[t] = t(alpha) %*% c(1,r[(t-1):(t-p)]^2) + t(beta) %*% sigma[(t-1):(t-q)]
    ll = ll + log(sigma[t]) + r[t]^2/(sigma[t])
  }
  return(ll)
}

# One-step-ahead forecast of sigma
sigmaForecast = function(alpha, beta, r){
  # One-step-ahead forecasts of the volatility sigma (Eq. 5.52)
  # alpha = (alpha_0,..., alpha_p): Estimated parameters by num. optim.
  # beta = (beta_1,..., beta_q): Estimated parameters by num. optim.
  # r: vector of returns.
  
  p = length(alpha) - 1
  q = length(beta)
  m = max(p,q)
  
  # Initialize first m conditional variances.
  # TODO: This is not correct initialization of the first sigma[1:m], but it works for now.
  sigma = numeric(n)
  sigma[1:m] = (r[1:m]^2)
  for (t in (m+1):n){
    sigma[t] = t(alpha) %*% c(1,r[(t-1):(t-p)]^2) + t(beta) %*% sigma[(t-1):(t-q)]
  }
  return(sigma)
}




# runs only when script is run by itself
# Test ----
if (sys.nframe() == 0){
  df <- read.csv('projectdata.csv', header = T, sep=";", dec=",", stringsAsFactors=FALSE)
  # Comparison
  d = diff(df$Inflation)
  sAcf = sampleAcf(d, 22)
  acf(d, lag.max = 22)
  
  samplePacf(d)
  pacf(d)
}
