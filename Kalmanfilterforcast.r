set.seed(1)

n=50

make_y <- function(y0, y1){
  #n = 100
  w = rnorm(n,0,1)
  y <- c(y0,y1)
  for (t in 3:n) {
    y[t] <- 0.1 * y[t-1] + 0.9 * y[t-2] + 0.5*w[t]
  }
  return(y)
}

par(mfrow=c(1,1))
y <- make_y(0,0)
plot(c(1:n), y, xlab = "timesteps", type="b", main = bquote(y[t] == 0.1*y[t-1] + 0.9*y[t-2] + 0.5*w[t]))

forecast <- function(obs, init){  # init: (alpha1, alpha2, sigma1)

  tt <- length(obs)

  # make obj function
  obj <- function(param){
    alpha1 <- param[1]  # in front of y_{t-1}
    alpha2 <- param[2]  # in front of y_{t-2}
    sigma <- param[3]   # in front of w_t

    # define the A matrix, 3x3:
    A <- matrix(0,3,3)
    diag(A[2:3,1:2]) <- 1
    A[1,1] <- param[1]  # alpha1
    A[1,2] <- param[2]  # alpha2

    # observation matrix
    B <- matrix(0,1,3)
    B[1,2] <- alpha1
    B[1,3] <- alpha2

    # system noise
    H <- matrix(0,3,1)
    H[1,1] <- sigma
    HH <- H %*% t(H)

    # observation noise
    D <- sigma
    DD <- D * D

    # initial condition
    x <- matrix(0,3,1)
    S <- diag(3)

    # log likelihood
    ll <- 0

    for (t in 1:tt) {
      # compute innovations
      I <- y[t] - B%*%x

      # innovation variance
      C <- B %*% S %*% t(B) + DD

      # compute M matrix
      M <- S %*% t(B) %*% solve(C)

      # update x, S
      S <- A %*% (diag(3) - M%*%B) %*% S %*% t(A) + HH
      x <- A %*% (x + M%*%I)

      # add to ll
      ll = ll + log(C) + I^2/C

    } # end of for-loop

    return(ll)

  } # end of obj function


  # minimize obj function using optim
  opt <- optim(par=init, fn=obj, method="BFGS",
               control=list(maxit=300, reltol=1e-8))
  # check if convergence
  stopifnot(opt$convergence==0)  #print(opt$message)
  # parameter estimates
  param <- opt$par


  # filter and forecast function
  filterandforecast <- function(param){
    # define vectors
    res <- c() #matrix(ncol=0, nrow=1)
    sds <- c() #matrix(ncol=0, nrow=1)
    sds2 <- c() #matrix(ncol=0, nrow=1)
    scaled_innov <- c() #matrix(ncol=0, nrow=1)

    # system matrix
    A <- matrix(0,3,3)
    diag(A[2:3,1:2]) <- 1
    A[1,1] <- param[1]
    A[1,2] <- param[2]

    # observation matrix
    B <- matrix(0,1,3)
    B[1,2] <- param[1]
    B[1,3] <- param[2]

    # system noise
    H <- matrix(0,3,1)
    H[1,1] <- param[3]
    HH <- H %*% t(H)

    # observation noise
    D <- param[3]
    DD <- D * D

    # init conditions
    #x <- matrix(obs[3:1],3,1)
    #x <- matrix(obs[1:3],3,1)
    x <- matrix(0,3,1)
    S <- diag(3) * var(obs)
    #S <- diag(3)

    # Filter
    for (t in 1:tt) {
      res <- c(res, B%*%x)
      sds <- c(sds, sqrt(B %*% S %*% t(B) + DD))
      sds2 <- c(sds2, sqrt(B %*% S %*% t(B)))

      # innovations
      I <- obs[t] - B%*%x

      # innovations variance
      C <- B %*% S %*% t(B) + DD

      # scaled innovations
      scaled_innov <- c(scaled_innov, I/sqrt(C))

      # projection matrix, M
      M <- S %*% t(B) %*% solve(C)

      # update S
      S <- A %*% (diag(3)-M%*%B) %*% S %*% t(A) + HH
      print(S)

      # update x
      x <- A %*% (x + M%*%I)
    } # end of filter for-loop


    # Forecast
    for (t in 1:20) {
      x <- A %*% x
      S <- A %*% S %*% t(A) + HH

      res <- c(res, B%*%x)
      sds <- c(sds,  sqrt(B %*% S %*% t(B) + DD))
      sds2 <- c(sds2,  sqrt(B %*% S %*% t(B)))
      scaled_innov <- c(scaled_innov, NA)
    } # end of forecast for-loop

    df <- data.frame("res"=res, "sds"=sds, "sds2"=sds2, "scaled_innov"=scaled_innov)
    return(df)

  } # end of filterandforecast function

  result <- filterandforecast(init)   # 'init' must be changed to 'param' if obj() function must be used

  # make plots

  # return values
  return(result)

} # end of forecast function


test <- forecast(y, c(0.1,0.9,0.5))
upper <- test$res+1.96*test$sds
lower <- test$res-1.96*test$sds

k = n+20

plot(c(1:n), y, type="l", ylim=c(-10,10), xlim=c(1,k))
lines(x=1:k, y=test$res, col="red")
lines(x=1:k, y=upper, col="gray", lty=2)
lines(x=1:k, y=lower, col="gray", lty=2)


B <- matrix(0,1,3)
B[1,2] <- 0.1
B[1,3] <- 0.9
x <- matrix(y[3:1],3,1)
B
x

B%*%x
