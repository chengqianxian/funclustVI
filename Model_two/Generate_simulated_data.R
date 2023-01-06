library(fda)

# get the true basis splines coefficients
get_true_m <- function(cluster_true_curve_data, x, nbasis){
  K <- nrow(cluster_true_curve_data)
  true_m <- matrix(0, nrow = K, ncol = nbasis)
  basisobj <- fda::create.bspline.basis(c(min(x), max(x)), nbasis)
  for (i in 1:K) {
    true_m[i, ] <- fda::smooth.basis(argvals = x, y = cluster_true_curve_data[i, ], 
                                     fdParobj = basisobj)$fd$coef
  }
  return(true_m)
}
  


# Case 7, three clusters, Scenario 7 in the paper
Case_7 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  for (i in 1:n)
  {
    Y[,i] = rnorm(1, 0, 0.4) -0.25 + (1/1.3)*sin(x*1.3) + x^3 + rnorm(length(x),0,0.2)
    Y[,i+n] = rnorm(1, 0, 0.4) + 1.25 + (1/1.2)*sin(x*1.3) + x^3 + rnorm(length(x),0,0.2)
    Y[,i+2*n] = rnorm(1, 0, 0.4) + 2.50 + (1/4)*sin(x*1.3) + x^3 + rnorm(length(x),0,0.2)
  }
  cluster_true_curve_data <- matrix(0, nrow = K, ncol = length(x))
  cluster_true_curve_data[1, ] <- -0.25 + (1/1.3)*sin(x*1.3) + x^3
  cluster_true_curve_data[2, ] <- 1.25 + (1/1.2)*sin(x*1.3) + x^3
  cluster_true_curve_data[3, ] <- 2.50 + (1/4)*sin(x*1.3) + x^3
  return(list(Y = t(Y),
              cluster_true_curve_data = cluster_true_curve_data))
}

# Case 8, 3 clusters, with smaller variances, Scenario 8 in the paper
Case_8 <- function(data_params) 
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4,nbasis=6)
  B <- getbasismatrix(x, basisBspline, nderiv=0)
  
  mu.1 <- c(1.5, 1, 1.8, 2.0, 1, 1.5) 
  mu.2 <- c(2.8, 1.4, 1.8, 0.5, 1.5, 2.5)
  mu.3 <- c(0.4, 0.6, 2.4, 2.6, 0.1, 0.4)
  mu.matrix <- matrix(c(mu.1, mu.2, mu.3), nrow = 3, byrow = T)
  random.intercept <- numeric(K * n)
  for (i in 1:n){
    random.intercept[i] = rnorm(1, 0, 0.05)
    random.intercept[i+n] = rnorm(1, 0, 0.05)
    random.intercept[i+2*n] = rnorm(1, 0, 0.05)
    Y[, i] = rep(random.intercept[i], length(x)) + B%*%mu.1 + rnorm(length(x), 0, 0.4)
    Y[, i+n] = rep(random.intercept[i+n], length(x)) + B%*%mu.2 + rnorm(length(x), 0, 0.4)
    Y[, i+2*n] = rep(random.intercept[i+2*n], length(x)) + B%*%mu.3 + rnorm(length(x), 0, 0.4)
  }
  return(list(Y = t(Y),
              cluster_true_coef = mu.matrix))
}


# Case 9, three clusters, but more crossing among clusters, Scenario 9 in the paper
Case_9 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4,nbasis=6)
  B <- getbasismatrix(x, basisBspline, nderiv=0)
  
  mu.1 <- c(1.5, 1, 1.8, 2.0, 1, 1.5) 
  mu.2 <- c(2.8, 1.4, 1.8, 0.5, 1.5, 2.5)
  mu.3 <- c(0.4, 0.6, 2.4, 2.6, 0.1, 0.4)
  mu.matrix <- matrix(c(mu.1, mu.2, mu.3), nrow = 3, byrow = T)
  random.intercept <- numeric(K * n)
  for (i in 1:n){
    random.intercept[i] = rnorm(1, 0, 0.3)
    random.intercept[i+n] = rnorm(1, 0, 0.3)
    random.intercept[i+2*n] = rnorm(1, 0, 0.3)
    Y[, i] = rep(random.intercept[i], length(x)) + B%*%mu.1 + rnorm(length(x), 0, 0.15)
    Y[, i+n] = rep(random.intercept[i+n], length(x)) + B%*%mu.2 + rnorm(length(x), 0, 0.15)
    Y[, i+2*n] = rep(random.intercept[i+2*n], length(x)) + B%*%mu.3 + rnorm(length(x), 0, 0.15)
  }
  return(list(Y = t(Y),
              cluster_true_coef = mu.matrix))
}


# Case 10, three clusters, Scenario 10 in the paper
Case_10 <- function(data_params) # 3 clusters
{
  K = 3
  x = data_params$x
  n = data_params$curves_per_cluster
  Y = matrix(0,length(x), K*n) ## each curve is in a row
  basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4,nbasis=6)
  B <- getbasismatrix(x, basisBspline, nderiv=0)
  
  mu.1 <- c(1.5, 1, 1.8, 2.0, 1, 1.5) 
  mu.2 <- c(2.8, 1.4, 1.8, 0.5, 1.5, 2.5)
  mu.3 <- c(0.4, 0.6, 2.4, 2.6, 0.1, 0.4)
  mu.matrix <- matrix(c(mu.1, mu.2, mu.3), nrow = 3, byrow = T)
  random.intercept <- numeric(K * n)
  for (i in 1:n){
    random.intercept[i] = rnorm(1, 0, 0.6)
    random.intercept[i+n] = rnorm(1, 0, 0.6)
    random.intercept[i+2*n] = rnorm(1, 0, 0.6)
    Y[, i] = rep(random.intercept[i], length(x)) + B%*%mu.1 + rnorm(length(x), 0, 0.4)
    Y[, i+n] = rep(random.intercept[i+n], length(x)) + B%*%mu.2 + rnorm(length(x), 0, 0.4)
    Y[, i+2*n] = rep(random.intercept[i+2*n], length(x)) + B%*%mu.3 + rnorm(length(x), 0, 0.4)
  }
  return(list(Y = t(Y),
              cluster_true_coef = mu.matrix))
}



