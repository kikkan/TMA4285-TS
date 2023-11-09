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

plot_acf <- function(h,x){
  #Function to plot the sample autocorrelation.
  #h: lag to evaluate ACF up to
  #x: time series
  lags <- 0:h
  result <- sapply(lags,sample_acf,x = x)
  result[1] <- result[1] / sample_acf(0,x)
  CI <- c(1.96/sqrt(length(x)), -1.96/sqrt(length(x)))
  plot(lags, result, type="h", main="Sample autocorrelation", ylab = "ACF", xlab = "Lag", ylim = c(min(CI[2],result),1)) +
    abline(h=0) + abline(h=CI, col=c("blue","blue"), lty=c(2,2))

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


PACF <- function(h,x){
  #Computes and plots the partial autocorrelation function with a 95% confidence interval.
  #h: lag to compute PACF up to.
  #x: time series

  coeffs <- durbin_levinson(h,x)
  lower_lim <- -1.96/sqrt(length(x))
  upper_lim <- 1.96/sqrt(length(x))

  plot(1:h,coeffs,xlab = "LAG",ylab = "PACF",type="h")
  points(1:h,coeffs, xlab="LAG", pch=16,col="red", type = "p")
  abline(h = lower_lim,col = "blue",lty = "dashed")
  abline(h = upper_lim,col = "blue", lty = "dashed")
}

# denne kan fjernes?
ACF_plot <- function(x, h){  # h = number of lags wanted
  sample_corr <- numeric(h)
  h = c(0:h)
  x <- na.omit(x)
  x_mean <- mean(x)
  n <- length(x)

  for (i in h) {
    cov <- 0
    for (t in 1:(n-i)) {
      cov <- cov + ( (x[t+i]-x_mean) * (x[t]-x_mean) )
    }
    corr <- cov/n
    sample_corr[i+1] <- corr
  }

  sample_corr <- sample_corr/sample_corr[1]
  CI <- c(1.96/sqrt(n), -1.96/sqrt(n))

  p <- plot(h, sample_corr, type="h", main="Sample autocorrelation", ylab = "ACF", xlab = "Lag", ylim = c(min(CI[2],sample_corr),1)) +
    abline(h=0) + abline(h=CI, col=c("blue","blue"), lty=c(2,2))

  print(p)
}


difference <- function(data, lag = 12){
  len <- length(data)
  diff_data <- vector("numeric",length=len-lag)

  for (i in lag:len) {
    diff_data[i-lag] <- data[i] - data[i-lag]
  }
  return(diff_data)
}


