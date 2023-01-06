Case_1 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = runif(1,-1/4,1/4) + .3+ (1/1.3)*sin(x*1.3) + x^3 + rnorm(length(x),0,0.4)
    Y[,i+n] = runif(1,-1/4,1/4) + 1.0+ (1/1.2)*sin(x*1.3) + x^3 + rnorm(length(x),0,0.4)
    Y[,i+2*n] = runif(1,-1/4,1/4) + .2+ (1/4)*sin(x*1.3) + x^3 + rnorm(length(x),0,0.4)
  }
  return(t(Y))
}

Case_2 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = runif(1,-1/4,1/4) + (1/1.8)*exp(x*1.1) - x^3+ rnorm(length(x),0,0.3)
    Y[,i+n] = runif(1,-1/4,1/4) + (1/1.7)*exp(x*1.4) - x^3+ rnorm(length(x),0,0.3)
    Y[,i+2*n] = runif(1,-1/4,1/4) + (1/1.5)*exp(x*1.5) - x^3+ rnorm(length(x),0,0.3)
  }
  return(t(Y))
}


Case_3 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4,nbasis=6)
  B <-getbasismatrix(x, basisBspline, nderiv=0)
  
  mu.1 <- c(1.5, 1, 1.8, 2.0, 1, 1.5)
  mu.2 <- c(2.8, 1.4, 1.8, 0.5, 1.5, 2.5)
  mu.3 <- c(0.4, 0.6, 2.4, 2.6, 0.1, 0.4)
  
  for (i in 1:n)
  {
    Y[, i] = B%*%mu.1 + rnorm(length(x), 0, 0.4)
    Y[, i+n] = B%*%mu.2 + rnorm(length(x), 0, 0.4)
    Y[, i+2*n] = B%*%mu.3 + rnorm(length(x), 0, 0.4)
  }
  return(t(Y))
}


Case_4 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4,nbasis=6)
  B <-getbasismatrix(x, basisBspline, nderiv=0)
  # following would be the true mean curves
  coef.1 <- c(1.5, 1, 1.6, 1.8, 1, 1.5)
  coef.2 <- c(1.8, 0.6, 0.4, 2.6, 2.8, 1.6)
  coef.3 <- c(1.2, 1.8, 2.2, 0.8, 0.6, 1.8)
  for (i in 1:n)
  {
    Y[, i] = B%*%coef.1 + rnorm(length(x), 0, 0.4)
    Y[, i+n] = B%*%coef.2 + rnorm(length(x), 0, 0.4)
    Y[, i+2*n] = B%*%coef.3 + rnorm(length(x), 0, 0.4)
  }
  return(t(Y))
}


Case_5 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = 0.1 * (0.4 + exp(-(x - 6)^2 / 3) + 0.2 * exp(-(x - 12)^2 / 25) + 0.5 * exp(-(x - 19)^2 / 4)) + rnorm(length(x), 0, 0.012)
    Y[,i+n] = 0.1 * (0.2 + exp(-(x - 5)^2 / 4) + 0.25 * exp(-(x - 18)^2 / 5)) + rnorm(length(x), 0, 0.012)
    Y[,i+2*n] = 0.1 * (0.2 +  exp(-(x - 3)^2 / 4) + 0.25 * exp(-(x - 16)^2 / 5)) + rnorm(length(x), 0, 0.012)
  }
  return(t(Y))
}


Case_6 <- function(data_params) # 4 clusters
{
  K = 4
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = runif(1,-1/3,1/3) +0.2- sin(pi*x*1.1) +x^3 + rnorm(length(x),0,0.4)
    Y[,i+n] = runif(1,-1/3,1/3) +0.5- sin(pi*x*1.4) +x^3 + rnorm(length(x),0,0.4)
    Y[,i+2*n] = runif(1,-1/3,1/3) +0.7- sin(pi*x*1.6) +x^3 + rnorm(length(x),0,0.4)
    Y[,i+3*n] = runif(1,-1/3,1/3) +1.3- sin(pi*x*1.8) +x^3 + rnorm(length(x),0,0.4)
  }
  return(t(Y))
}

