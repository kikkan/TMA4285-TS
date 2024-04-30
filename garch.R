library(viridis)

figPath = "./figures/"

# Garch(p,q) implementation ----
garchLL = function(params, r, p=1, q=0){
  # LogLike of alpha and beta given the first r_t values t=(1,..., max(p,q))(Eq. 5.46)
  # Defaults to GARCH(1,1) and uses only the first three params (alpha0, alpha1, beta1)
  # params: vector of log(parameters). First p+1 are alphas, rest are betas
  # r: vector of returns (Eq. 5.34)
  # p,q: GARCH(p,q) model parameters
  # -> returns -log likelihood
  
  n = length(r)
  m = max(p,q)
  
  alpha = exp(params[1:(p+1)])
  if (q==0){beta = 0}
  else {beta = exp(params[-(1:(p+1))]) }
  
  # initialize first m values of the variance
  sigma_sq = numeric(n)
  sigma_sq[1:m] = t(alpha) %*% c(1,r[p:1]^2) # should they be zero?
  
  # Iteratively compute each variance
  for (t in (m+1):n){
    sigma_sq[t] = sum(alpha * c(1,r[(t-1):(t-p)]^2)) + sum(beta*sigma_sq[(t-1):(t-q)]) 
  }

  return(sum(log(sigma_sq) + r^2/sigma_sq))
}


sigmaForecast = function(mod){
  # One-step-ahead forecasts of the volatility sigma (Eq. 5.52)
  # mod: list of estimated parameters and the data it was fitted on.
  # -> return one-step-ahead forecasts
  
  # Extract values for readability.
  n = mod$n
  m = mod$m
  p = mod$p
  r = mod$r
  
  alpha = mod$estim[1:(p+1)]
  if (mod$q==0){beta = 0; q=1} # If ARCH(p) model
  else{beta = mod$estim[-(1:(p+1))]; q = mod$q}
  
  sf = numeric(n) # sigma squared forecasts
  sf[1:m] = sum(alpha*c(1,r[1:p]^2)) # first m values
  
  for (t in (m+1):n){
    sf[t] = sum(alpha * c(1,r[(t-1):(t-p)]^2)) + sum(beta*sf[(t-1):(t-q)])
  }
  return(sf)
}


Garch = function(r, p=1, q=0, init=c(0.1,0.1)){
  # Make, estimate and forecast volatility of a GARCH(p,q) given the parameters
  # r: the returns
  # p,q: GARCH(p,q)
  # init: initial parameters for alphas and betas
  
  n = length(r)
  rq = q
  
  alpha = init[1:(p+1)]
  if (q==0){beta = 0; q=1} # If ARCH(p) model
  else{beta = init[-(1:(p+1))]; q = q}
  
  estim = exp(optim(par = log(init), fn = garchLL, r = r, p=p, q=q, method = "BFGS")$par)
  
  mod = list(
    p = p,
    q = rq,
    n = n,
    m = max(p,q),
    init = init,
    estim = estim,
    r=r
  )
  
  mod$sigmaForecast = sigmaForecast(mod)
  return(mod)
  
}


# Test functions ----
testModel = function(mod){
  # Tests model compared to sampled data given real values of the alphas and the
  # betas.
  
  n = mod$n
  m = mod$m
  p = mod$p
  
  alpha = mod$paramSim[1:(p+1)]
  if (mod$q==0){beta = 0; q=1}
  else{beta = mod$paramSim[-(1:(p+1))]; q = mod$q}
  
  r = numeric(n)
  r[1:m] = rnorm(m)
  sigma_sq = numeric(n)
  sigma_sq[1:m] = sum(alpha*c(1,r[1:p]^2))
  
  for (t in (m + 1):n){
    sigma_sq[t] = sum(alpha*c(1,r[(t-1):(t-p)]^2)) + sum(beta*sigma_sq[(t-1):(t-q)])
    r[t] = rnorm(1, sd=sqrt(sigma_sq[t]))
  }
  
  # estimation
  estim = exp(optim(par = log(mod$init), fn = garchLL, r = r, p=p, q=q, method = "BFGS")$par)
  comparison = (rbind(real = mod$paramSim, estim = estim, init = mod$init))
  return(list(r = r, estim = estim, comparison = comparison))
}


# Convenience functions ----
modResults = function(mod, main="", plot=T){
  if (plot){
    plot(mod$r, type="l", main=main, xlab = "Month", 
         ylab = "Volatility [%]", cex.main=0.8, cex.axis = 0.8, cex.lab = 0.8)
    lines(sqrt(mod$sigmaForecast), col="cyan")
  }
}

# Tests ----
set.seed(420)
# ARCH(1) Test
arch1 = list(
  p = 1, 
  q= 0,
  m= 1,
  n = 1e4,
  init = rep(0.1,2),
  paramSim = c(0.01, 0.2)
)


testResult = testModel(arch1)
testResult$comparison

# GARCH(1,2) test
garch12 = list(
  p = 1, 
  q= 2,
  m= 2,
  n = 1e4,
  init = rep(0.1,4),
  paramSim = c(0.01, 0.2, 0.1, 0.5)
)


testResult = testModel(garch12)
garch12$estim = testResult$estim
round(testResult$comparison,3)


# Real data ----
df <- read.csv('projectdata.csv', header = T, sep=";", dec=",", stringsAsFactors=FALSE)

# fit GARCH(1,2) and plot results
fit = lm(Inflation~Unemployed+Consumption+InterestRate, data = df)
r = df$Inflation - fit$fitted.values
d = diff(df$Inflation)
par(mfrow=c(2,1))

## inflation - regression
garch12regr = Garch(r, p=1, q=2, init=rep(0.1, 1+2+1))
modResults(garch12regr, main="garch(1,2) on regression")

## diff inflation
garch12 = Garch(d, p=1, q=2, init=rep(0.1, 1+2+1))
modResults(garch12, "garch(1,2) on diff")
par(mfrow=c(1,1))

# fit several ARCH(p)
arch1regr = Garch(r, p=1, q=0, init=rep(0.1, 2))
arch1d = Garch(d, p=1, q=0, init=rep(0.1, 2))
arch2regr = Garch(r, p=2, q=0, init=rep(0.1, 3))
arch2d = Garch(d, p=2, q=0, init=rep(0.1, 3))
arch10regr = Garch(r, p=10, q=0, init=rep(0.1, 11))
arch10d = Garch(d, p=10, q=0, init=rep(0.1, 11))

## plot
par(mfrow=c(3,2))
modResults(arch1regr, main="arch(1) on inflation - regression")
modResults(arch1d, main = "arch(1) on diff(inflation)")
modResults(arch2regr, main="arch(2) on inflation - regression")
modResults(arch2d, "arch(2) on diff(inflation)")
modResults(arch10regr, main="arch(10) on inflation - regression")
modResults(arch10d, "arch(10) on diff(inflation)")
par(mfrow=c(1,1))

if (T){
  pdf(file = paste0(figPath, "garchVolatility.pdf"))
  par(mfrow=c(3,2))
  modResults(arch1regr, main="")
  modResults(arch1d, main = "")
  modResults(arch2regr, main="")
  modResults(arch2d, "")
  modResults(arch10regr, main="")
  modResults(arch10d, "")
  par(mfrow=c(1,1))
  dev.off()
}

## compare parameter estimates
list(
  arch1regr = arch1regr$estim,
  arch1d = arch1d$estim,
  arch2regr = arch2regr$estim,
  arch2d = arch2d$estim,
  arch10regr = arch10regr$estim,
  arch10d = arch10d$estim
)

models = list(
  arch1d = arch1d,
  arch2d = arch2d,
  arch10d = arch10d,
  arch1regr = arch1regr,
  arch2regr = arch2regr,
  arch10regr = arch10regr,
  garch12regr = garch12regr,
  garch12 = garch12
)


# Model selection ----
logLike = function(mod, m=1){
  i = (m+1):(mod$n) # can burn m first values
  ll = sum(-0.5*(log(mod$sigmaForecast[i]) + mod$r[i]^2/mod$sigmaForecast[i]))
  return(ll)
}


aic = function(mod, m=1){
  return(2*(mod$p + mod$q) -2*logLike(mod, m))
}

loglikes = numeric(length(models))
aics = numeric(length(models))
loglikesM = numeric(length(models))
aicsM = numeric(length(models))

for (i in 1:length(models)){
  loglikes[i] = logLike(models[[i]])
  aics[i] = aic(models[[i]])
  loglikesM[i] = logLike(models[[i]], models[[i]]$m)
  aicsM[i] = aic(models[[i]], models[[i]]$m)
  # print(sum(models[[i]]$r^2/models[[i]]$sigmaForecast)) # Interesting!
}

modSel = cbind(ll = loglikes, llM = loglikesM, aic = aics, aicM = aicsM)
rownames(modSel) = names(models)
modSel

plot(arch1d$r^2, type = "l")
lines(arch1d$sigmaForecast, col="cyan")
sum(arch1d$r^2/arch1d$sigmaForecast)
length(arch1d$sigmaForecast)


# return forecast ----
forecast = function(mod, n){
  m = mod$m
  p = mod$p
  r = mod$r
  nTot = mod$n + n
  
  alpha = mod$estim[1:(p+1)]
  if (mod$q==0){beta = 0; q=1} # If ARCH(p) model
  else{beta = mod$estim[-(1:(p+1))]; q = mod$q}
  
  rf = r
  length(rf) = nTot
  sf = mod$sigmaForecast
  length(sf) = nTot
  for (t in mod$n:(nTot)){
    sf[t+1] = sum(alpha * c(1,rf[(t+1-1):(t+1-p)]^2)) + sum(beta*sf[(t+1-1):(t+1-q)])
    rf[t+1] = rnorm(1, sd = sqrt(sf[t+1]))
  }
  return(rf)
}

plotForecast = function(models){
  c = viridis_pal(option = "D")(length(models))
  
  for (i in 1:length(models)){
    nTot = length(models[[i]]$forecast)
    w = models[[i]]$w
    lines((nTot-without):nTot, models[[i]]$forecast[(nTot - without):nTot], col=c[i], lty = 2)
  }
  return(c)
}

set.seed(420)
without = 12
dTot = diff(df$Inflation)
nTot = length(dTot)
d = diff(dTot[1:(nTot-without)])


arch1df = Garch(d, p=1, q=0, init=rep(0.1, 2))

arch2df = Garch(d, p=2, q=0, init=rep(0.1, 3))

arch10df = Garch(d, p=10, q=0, init=rep(0.1, 11))

modelsf = list(
  arch1df = arch1df,
  arch2df = arch2df,
  arch10df = arch10df
)



# for (mod in modelsf){
for (i in 1:length(modelsf)){
  # browser()
  w = 12
  modelsf[[i]]$w = w
  # mod$w = w
  modelsf[[i]]$forecast = forecast(modelsf[[i]], w)
}


back = 50
plot((nTot-back):nTot,dTot[(nTot-back):nTot], type = "l", xlab = "Month", ylab = "Volatility [%]", cex.main=0.8, cex.axis = 0.8, cex.lab = 0.8)
c = plotForecast(modelsf)
legend(160, -1, legend=c("data", "ARCH(1)", "ARCH(2)", "ARCH(10)"),
       col=c("#000000", c), lty=1:1:1, cex=0.5)


pdf(file = paste0(figPath, "garchForecasts.pdf"), height = 4)
plot((nTot-back):nTot,dTot[(nTot-back):nTot], type = "l", xlab = "Month", ylab = "Volatility [%]", cex.main=0.8, cex.axis = 0.8, cex.lab = 0.8)
c = plotForecast(modelsf)
legend(155, -1.05, legend=c("data", "ARCH(1)", "ARCH(2)", "ARCH(10)"),
       col=c("#000000", c), lty=1:1:1, cex=0.5)
dev.off()

forecastRes = numeric(length(modelsf))
for (i in 1:length(modelsf)){
  w = modelsf[[i]]$w
  forecastRes[i] = sum(dTot[(nTot-w):nTot] - modelsf[[i]]$forecast[(nTot-w):nTot])
  
}
forecastRes^2
# 9.590208  2.636010 11.154583