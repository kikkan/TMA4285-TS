df <- read.csv('projectdata.csv', header = T, sep=";", dec=",", stringsAsFactors=FALSE)
head(df)
df2 <- read.csv('projectdata(1).csv', header = T, sep=";", dec=",", stringsAsFactors=FALSE)
source("functions.r")

df$Inflation
model <- lm(Inflation~InterestRate+Consumption+Bankrupt, data=df)
x_series <- df$Inflation - model$fitted.values

# Plot residuals
year_divider <- c(0:15)*12
plot(c(1:length(x_series)), x_series, type="l", main="residuals", xlab="month", ylab="error") + abline(v=year_divider, lty=2, lw=1, col="gray")

# Need to look at ACF and PACF, to identify ARMA
plot_acf(100,x_series)
PACF(100,x_series)

# Difference
diff_x <- difference(x_series,lag=12)
plot(c(1:length(diff_x)), diff_x, type="l", main="residuals", xlab="month", ylab="error") + abline(v=year_divider, lty=2, lw=1, col="gray")
plot_acf(20,diff_x)
PACF(20,diff_x)
qqnorm(diff_x)


kalman_filter <- function(x,y,T,A,H,B,D){
  #Initial conditions
  s <- length(y)
  states <- vector("numeric",length = T)
  sds <- vector("numeric",length = T)
  innovations <- vector("numeric",length = s)

  S <- diag(x = dim(A)[1])*var(y) #Initial covariance matrix

  #Filtering
  for (t in 1:(s)){
    states[t] <- B %*% x
    #Innovations
    I <- y[t] - states[t]
    I <- as.numeric(I)
    innovations[t] <- I

    #Projection matrix (Kalman gain)
    C <- (B %*% S) %*% t(B) + D %*% t(D)
    C <- as.numeric(C)
    sds[t] <- sqrt(C)

    M <- (S %*% t(B)) %*% solve(C)

    #Condition on innovation/observation
    x_tt <- x + M %*% I

    #Update x,S
    x <- A %*% x_tt

    S <- (diag(dim(A)[1]) - M %*% B) %*% S
    S <- A %*% S %*% t(A) + H %*% t(H)
  }

  #Prediction
  if (T > s){
    for (t in s:(T)){
      x <- A %*% x
      S <- A %*% S %*% t(A) + H %*% t(H)
      states[t] <- B %*% x

      C <- (B %*% S) %*% t(B) + D %*% t(D)
      C <- as.numeric(C)
      sds[t] <- sqrt(C)
    }
  }

  return(data.frame("states" = states,"sds" = sds,"innovations" = innovations))
}


find_parameters <- function(x,y,T,A,H,B,D,initial_params){
  #Function to estimate unknown parameters in SARIMA-model.
  likelihood <- function(params){
    #Initial conditions
    #Change according to model
    A[1,dim(A)[2]] <- params[1]

    s <- length(y)
    S <- diag(x = dim(A)[1])*var(y) #Initial covariance matrix

    ll <- 0
    #Filtering
    for (t in 1:(s)){
      #Innovations
      I <- y[t] - B %*% x
      I <- as.numeric(I)

      #Projection matrix (Kalman gain)
      C <- (B %*% S) %*% t(B) + D %*% t(D)
      C <- as.numeric(C)

      M <- (S %*% t(B)) %*% solve(C)

      #Condition on innovation/observation
      x_tt <- x + M %*% I

      #Update x,S
      x <- A %*% x_tt

      S <- (diag(dim(A)[1]) - M %*% B) %*% S
      S <- A %*% S %*% t(A) + H %*% t(H)

      ll <- ll + 0.5*(log(abs(C)) + I^2/C)
    }
    return (ll)
  }

  new_params <- optim(initial_params,likelihood,method = "BFGS",control = list(maxit = 300))[1]
  return (new_params)
}

initial_params <- c(1) #alpha

A <- diag(nrow = 11)
firstrow <- rep(0,11)
A <- rbind(firstrow,A)
lastcol <- rep(0,12)
lastcol[1] <- 1 #alpha
A <- cbind(A,lastcol)

H <- matrix(0,nrow = 12)
H[1] <- 1

B <- matrix(0,ncol = 12)
B[1] <- 1

D <- 0
x <- matrix(0,nrow = 12)
print(find_parameters(x,x_series,length(x_series),A,H,B,D,initial_params))


# Need to look at ACF and PACF, to identify GARCH?

# Finds parameters by maximum likelihood estimation
find_parameters_garch <- function(x, initial_params){
  likelihood_garch <- function(params){
    alpha0 <- exp(params[1])
    alpha1 <- exp(params[2])

    if (alpha0 <= 0 || alpha1 <= 0) {
      return(Inf) #TODO brute force exp params? not needed? ----
    }

    ll <- 0
    for (i in 2:length(x)){
      variance_term <- alpha0 + alpha1 * x[i-1]^2
      ll <- ll + 0.5 * (log(variance_term) + x[i]^2 / variance_term)
    }
    return(ll)
  }

  # Use logs of parameters for optimization to ensure positivity
  initial_params_log <- log(initial_params)
  new_params_log <- optim(initial_params_log, likelihood_garch, method = "BFGS", control = list(maxit = 500))$par

  # Transform parameters back to original scale
  new_params <- exp(new_params_log)
  return(new_params)
}

# Simulation:
# TEST
# Testing using some random rt-N(0,alpha0 + alpha1*r_t-1)...
alpha1 <- 0.3
alpha2 <- 0.7
n <- 200
r_vec <- rep(0, n)
# Draw from normal distribution
for (i in 2:n){
  r_vec[i] <- rnorm(1, sd = sqrt(alpha1 + alpha2*r_vec[i-1]^2))
}

initial_params_garch <- c(0.5, 0.1)
alpha_test <- find_parameters_garch(r_vec, initial_params_garch)
# Estimated parameters - We see that they are close to 0.3 and 0.7 -> Code seems to work well
alpha_test

# The estimated volatility
sigma_garch <- rep(0, n-1)
for (i in 2:n){
  sigma_garch[i] <- sqrt(alpha_test[1] + alpha_test[2]*r_vec[i-1]^2)
}

# Plotting to check if volatility captures variability of r
plot(1:length(r_vec),r_vec, type = "l", col = "blue")
lines(1:length(sigma_garch), sigma_garch, col = "red")


# ARCH(1)
# Save inflation in variable x_series, to test for different rt
x_series <- df$Inflation
# rt inflation data differenced
r <- rep(0, length(x_series)-1)
for (i in 2:length(x_series)){
  r[i-1] <- (x_series[i] - x_series[i-1])
}
# Finding estimated alpha's
initial_params_garch <- c(0.5,0.5)
alpha_garch <- find_parameters_garch(r, initial_params_garch)
alpha_garch

# Simulate the series using the estimated parameters
sigma_garch <- rep(0, length(r)-1)
for (i in 2:length(sigma_garch)){
  sigma_garch[i] <- sqrt(alpha_garch[1] + alpha_garch[2]*r[i-1]^2)
}
# Idk about initialization... Just set it like this for a smoother plot
sigma_garch[1] <- sigma_garch[2]

# Plot r_t and the volatility magnitude (not squared...) to see the trend more clearly
par(mfrow=c(1,2))
plot(1:length(r),r, type = "l", col = "blue", main = " Volatility of the differenced inflation rate", xlab = "Month", ylab = "Volatility [%]",
     xlim = c(0,175), cex.main=0.8, cex.axis = 0.8, cex.lab = 0.8)
lines(1:length(sigma_garch), sigma_garch, col = "red")
legend(0.1, 1.8, legend=c("Volatility σ", "r"),
       col=c("red", "blue"), lty=1:1, cex=0.5)
# QQ-plot for rt
qqnorm(r, main = "Normal Q-Q plot for r(t)")
qqline(r)


# Forecasting volatility
set.seed(1)
# Dataset to predict from
x.pred <- na.omit(df2$newInflation)[1:12]
r.pred <- rep(0, length(x.pred))
sigma.pred <- rep(0,length(x.pred))
# Predicting unknown r_t+1 and volatility
for (i in 2:length(r.pred)){
  r.pred[i] <- rnorm(1,sd = sqrt(alpha_garch[1] + alpha_garch[2]*r.pred[i-1]^2))
  sigma.pred[i] <- sqrt(alpha_garch[1] + alpha_garch[2]*r.pred[i-1]^2)
}

# The actual r_t and volatility in the testing set newInlation
x.actual <- na.omit(df2$newInflation)
r.actual <- rep(0, length(x.actual))
sigma.actual <- rep(0,length(x.actual))
# Initial value equal to last value in the dataset Inflation - for smoother plots
sigma.actual[1] <- sigma_garch[length(sigma_garch)]
r.actual[1] <- r[length(r)]

# We set the initial value of the predicted data equal to the actual data - smoother plots
r.pred[1] <- r[length(r)]
sigma.pred[1] <- sigma_garch[length(sigma_garch)]

# Simulating the actual r_t and volatility in the test set
for (i in 2:length(x.actual)){
  r.actual[i] <- (x.actual[i] - x.actual[i-1])
  sigma.actual[i] <- sqrt(alpha_garch[1] + alpha_garch[2]*r.actual[i-1]^2)
}

# Plotting the prediction with the actual volatility
par(mfrow=c(1,1))
plot(1:length(sigma_garch),sigma_garch^2, type = "l", col = "blue", main = "1 year prediction of volatility", xlab = "Month", ylab = "Volatility [%]",
     xlim = c(0,195), cex.main=0.8, cex.axis = 0.8, cex.lab = 0.8)
# x_values tells 'lines' where to place the test set - at the end
x_values <- seq(length(sigma_garch), length.out = length(r.pred), by = 1)
lines(x_values, sigma.pred^2, col = "red")
lines(x_values, sigma.actual^2, col = "green")
legend(0.1, 2.2, legend=c("actual σ^2", "predicted σ^2", "actual σ^2"),
       col=c("green", "red", "blue"), lty=1:1:1, cex=0.5)


