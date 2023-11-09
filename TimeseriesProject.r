setwd("C:\\Users\\celin\\OneDrive\\Dokumenter\\FYSMAT 7.semester\\Tidsrekker\\Project")
df <- read.csv('projectdata.csv', header = T, sep=";", dec=",", stringsAsFactors=FALSE)
head(df)
source("functions.r")

# Run ordinary regression of Y_t on z_t1,z_t2,.., and estimate the betas.
model <- lm(Inflation~InterestRate+Consumption+Bankrupt, data=df)
summary(model)

# Identify ARMA models using the residuals: x_t = Y_t - sum(beta*z_ti)
x <- df$Inflation - model$fitted.values

# Plot residuals
year_divider <- c(0:15)*12
plot(c(1:length(x)), x, type="l", main="residuals", xlab="month", ylab="error") + abline(v=year_divider, lty=2, lw=1, col="gray")

# Need to look at ACF and PACF, to identify ARMA
plot_acf(100,x)
PACF(100,x)

# Difference
diff_x <- difference(x,lag=12)
plot(c(1:length(diff_x)), diff_x, type="l", main="residuals", xlab="month", ylab="error") + abline(v=year_divider, lty=2, lw=1, col="gray")
plot_acf(20,diff_x)
PACF(20,diff_x)
qqnorm(diff_x)

# Least squares or MLE to find, beta, phi and theta

# Adjust residuals if needed


