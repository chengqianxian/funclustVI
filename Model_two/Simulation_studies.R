### Scenario 7, Case 7 ------------------------------------------------------------------------------------------------
set.seed(1)
x <- seq(from=0, to=1, by=0.01)
curves_per_cluster = 50 
data_params <- list(x = x, curves_per_cluster = curves_per_cluster)

Case.7.data <- Case_7(data_params)
Y <- Case.7.data$Y
cluster_true_curve_data <- Case.7.data$cluster_true_curve_data
K = 3
nbasis = 6
true_coeff <- get_true_m(cluster_true_curve_data, x, nbasis)
basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4,nbasis=6)
B <- getbasismatrix(x, basisBspline, nderiv=0)
Y.1.true <- B %*% true_coeff[1, ]
Y.2.true <- B %*% true_coeff[2, ]
Y.3.true <- B %*% true_coeff[3, ]


# test
Y <- Y
x <- seq(from=0, to=1, by=0.01)
K = 3
curves_per_cluster = 50
true_cluster_assignments <- rep(1:K,each = curves_per_cluster)
init = "km"
nbasis = 6
convergence_threshold <- 0.1
# convergence_threshold_2 <- 0.1
max_iterations <- 1000
b_not = 100 * 20
r_not = 4 * 20
d_not_vector = c(1/3, 1/3, 1/3)
m_not_matrix = true_coeff
v_not_vector = rep(10, nbasis)
alpha_not = 100 * 5
beta_not = 16 * 5
ext.model.test.result <- funcslustVI.ext(x, Y, K, true_cluster_assignments, init, nbasis, 
                                         convergence_threshold, max_iterations, 
                                         b_not, r_not, d_not_vector, m_not_matrix, v_not_vector,
                                         alpha_not, beta_not)

ext.model.test.result$iteration
ext.model.test.result$cluster_assignments
predicited_clusters <- ext.model.test.result$cluster_assignments
true_clusters <- true_cluster_assignments
Mismatch.case.7 <- Mismatch(predicited_clusters, true_clusters, K)
Vmeasure.case.7 <- sabre::vmeasure(predicited_clusters, true_cluster_assignments)$v_measure


## plot

plot(x, Y[1, ], col = "purple", type="l", lty = 2, lwd = 0.5, ylim=c(-1.5, 5), ylab="f(x)", xlab="x")
for (i in 2:50) {
  lines(x, Y[i, ], col = "purple", type="l", lty = 2, lwd = 0.5)
}

for (i in 51:100) {
  lines(x, Y[i, ], col = "blue", type="l", lty = 2, lwd = 0.5)
}

for (i in 101:150) {
  lines(x, Y[i, ], col = "green", type="l", lty = 2, lwd = 0.5)
}
lines(x, Y.1.true, col = "black", lwd = 3)
lines(x, Y.2.true, col = "black", lwd = 3)
lines(x, Y.3.true, col = "black", lwd = 3)

Y.1.est <- B %*% t(ext.model.test.result$m_list[[1]])
Y.2.est <- B %*% t(ext.model.test.result$m_list[[2]])
Y.3.est <- B %*% t(ext.model.test.result$m_list[[3]])

lines(x, Y.1.est, col = "red", lwd = 3)
lines(x, Y.2.est, col = "red", lwd = 3)
lines(x, Y.3.est, col = "red", lwd = 3)
legend(x= "topleft", cex=.7, legend=c("True Mean Curves", "VB Estimated Mean Curves"), col=c("black", "red"), lty=c(1,1))


### Scenario 8, case 8 ------------------------------------------------------------------------------------------------
x <- seq(from=0, to=1, by=0.01)
curves_per_cluster = 50 
data_params <- list(x = x, curves_per_cluster = curves_per_cluster)

set.seed(1)
Case_8_data <- Case_8(data_params)
Y <- Case_8_data$Y


# test
Y <- Y
x <- seq(from=0, to=1, by=0.01)
K = 3
curves_per_cluster = 50
true_cluster_assignments <- rep(1:K,each = curves_per_cluster)
init = "km"
nbasis = 6
convergence_threshold <- 0.01

max_iterations <- 1000
b_not = 78.125 * 30
r_not = 12.5 * 30
d_not_vector = c(1/3, 1/3, 1/3)
m_not_matrix = Case_3_data$cluster_true_coef
v_not_vector = rep(10, nbasis)
alpha_not = 600 * 300
beta_not = 1.5 * 300
ext.model.test.3.result <- funcslustVI.ext(x, Y, K, true_cluster_assignments, init, nbasis, 
                                           convergence_threshold, max_iterations, 
                                           b_not, r_not, d_not_vector, m_not_matrix, v_not_vector,
                                           alpha_not, beta_not)

ext.model.test.3.result$iteration
ext.model.test.3.result$cluster_assignments

predicited_clusters <- ext.model.test.3.result$cluster_assignments
true_clusters <- true_cluster_assignments
Mismatch(predicited_clusters, true_clusters, K)
sabre::vmeasure(predicited_clusters, true_cluster_assignments)$v_measure

## plot
mu.1 <- c(1.5, 1, 1.8, 2.0, 1, 1.5) 
mu.2 <- c(2.8, 1.4, 1.8, 0.5, 1.5, 2.5)
mu.3 <- c(0.4, 0.6, 2.4, 2.6, 0.1, 0.4)
basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4,nbasis=6)
B <- getbasismatrix(x, basisBspline, nderiv=0)
Y.1.true <- B %*% mu.1
Y.2.true <- B %*% mu.2
Y.3.true <- B %*% mu.3

plot(x, Y[1, ], col = "purple", type="l", lty = 2, lwd = 0.5, ylim=c(-0.5, 4), ylab="f(x)", xlab="x")
for (i in 2:50) {
  lines(x, Y[i, ], col = "purple", type="l", lty = 2, lwd = 0.5)
}

for (i in 51:100) {
  lines(x, Y[i, ], col = "blue", type="l", lty = 2, lwd = 0.5)
}

for (i in 101:150) {
  lines(x, Y[i, ], col = "green", type="l", lty = 2, lwd = 0.5)
}

lines(x, Y.1.true, col = "black", lwd = 3)
lines(x, Y.2.true, col = "black", lwd = 3)
lines(x, Y.3.true, col = "black", lwd = 3)



Y.1.est <- B %*% t(ext.model.test.3.result$m_list[[1]])
Y.2.est <- B %*% t(ext.model.test.3.result$m_list[[2]])
Y.3.est <- B %*% t(ext.model.test.3.result$m_list[[3]])

lines(x, Y.1.est, col = "red", lwd = 3)
lines(x, Y.2.est, col = "red", lwd = 3)
lines(x, Y.3.est, col = "red", lwd = 3)
legend(x= "topleft", cex=.7, legend=c("True Mean Curves", "VB Estimated Mean Curves"), col=c("black", "red"), lty=c(1,1))



### Scenario 9, Case 9 ------------------------------------------------------------------------------------------------
x <- seq(from=0, to=1, by=0.01)
curves_per_cluster = 50 
data_params <- list(x = x, curves_per_cluster = curves_per_cluster)

set.seed(1)
Case_9_data <- Case_9(data_params)
Y <- Case_9_data$Y
  

# test
Y <- Y
x <- seq(from=0, to=1, by=0.01)
K = 3
curves_per_cluster = 50
true_cluster_assignments <- rep(1:K,each = curves_per_cluster)
init = "km"
nbasis = 6
convergence_threshold <- 0.01
# convergence_threshold_2 <- 0.1
max_iterations <- 1000
b_not = 176 * 30
r_not = 4 * 30
d_not_vector = c(1/3, 1/3, 1/3)
m_not_matrix = Case_9_data$cluster_true_coef
v_not_vector = rep(10, nbasis)
alpha_not = 110 * 10
beta_not = 10 * 10
ext.model.test.2.result <- funcslustVI.ext(x, Y, K, true_cluster_assignments, init, nbasis, 
                                         convergence_threshold, max_iterations, 
                                         b_not, r_not, d_not_vector, m_not_matrix, v_not_vector,
                                         alpha_not, beta_not)

ext.model.test.2.result$iteration
ext.model.test.2.result$cluster_assignments

predicited_clusters <- ext.model.test.2.result$cluster_assignments
true_clusters <- true_cluster_assignments
Mismatch.case.9 <- Mismatch(predicited_clusters, true_clusters, K)
Vmeasure.case.9 <- sabre::vmeasure(predicited_clusters, true_cluster_assignments)$v_measure

## plot

mu.1 <- c(1.5, 1, 1.8, 2.0, 1, 1.5) 
mu.2 <- c(2.8, 1.4, 1.8, 0.5, 1.5, 2.5)
mu.3 <- c(0.4, 0.6, 2.4, 2.6, 0.1, 0.4)
basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4, nbasis=6)
B <- getbasismatrix(x, basisBspline, nderiv=0)
Y.1.true <- B %*% mu.1
Y.2.true <- B %*% mu.2
Y.3.true <- B %*% mu.3

plot(x, Y[1, ], col = "purple", type="l", lty = 2, lwd = 0.5, ylim=c(-0.5, 4), ylab="f(x)", xlab="x")
for (i in 2:50) {
  lines(x, Y[i, ], col = "purple", type="l", lty = 2, lwd = 0.5)
}

for (i in 51:100) {
  lines(x, Y[i, ], col = "blue", type="l", lty = 2, lwd = 0.5)
}

for (i in 101:150) {
  lines(x, Y[i, ], col = "green", type="l", lty = 2, lwd = 0.5)
}

lines(x, Y.1.true, col = "black", lwd = 3)
lines(x, Y.2.true, col = "black", lwd = 3)
lines(x, Y.3.true, col = "black", lwd = 3)



Y.1.est <- B %*% t(ext.model.test.2.result$m_list[[1]])
Y.2.est <- B %*% t(ext.model.test.2.result$m_list[[2]])
Y.3.est <- B %*% t(ext.model.test.2.result$m_list[[3]])

lines(x, Y.1.est, col = "red", lwd = 3)
lines(x, Y.2.est, col = "red", lwd = 3)
lines(x, Y.3.est, col = "red", lwd = 3)
legend(x= "topleft", cex=.7, legend=c("True Mean Curves", "VB Estimated Mean Curves"), col=c("black", "red"), lty=c(1,1))


### Scenario 10, Case 10 ------------------------------------------------------------------------------------------------
x <- seq(from=0, to=1, by=0.01)
curves_per_cluster = 50 
data_params <- list(x = x, curves_per_cluster = curves_per_cluster)

set.seed(1)
Case_10_data <- Case_10(data_params)
Y <- Case_10_data$Y


# test
Y <- Y
x <- seq(from=0, to=1, by=0.01)
K = 3
curves_per_cluster = 50
true_cluster_assignments <- rep(1:K,each = curves_per_cluster)
init = "km"
nbasis = 6
convergence_threshold <- 0.01
# convergence_threshold_2 <- 0.1
max_iterations <- 1000
b_not = 78.125 * 30
r_not = 12.5 * 30
d_not_vector = c(1/3, 1/3, 1/3)
m_not_matrix = Case_10_data$cluster_true_coef
v_not_vector = rep(10, nbasis)
alpha_not = 28 * 10
beta_not = 10 * 10
ext.model.test.4.result <- funcslustVI.ext(x, Y, K, true_cluster_assignments, init, nbasis, 
                                           convergence_threshold, max_iterations, 
                                           b_not, r_not, d_not_vector, m_not_matrix, v_not_vector,
                                           alpha_not, beta_not)

ext.model.test.4.result$iteration
ext.model.test.4.result$cluster_assignments

predicited_clusters <- ext.model.test.4.result$cluster_assignments
true_clusters <- true_cluster_assignments
Mismatch.case.10 <- Mismatch(predicited_clusters, true_clusters, K)
Vmeasure.case.10 <- sabre::vmeasure(predicited_clusters, true_cluster_assignments)$v_measure

## plot

mu.1 <- c(1.5, 1, 1.8, 2.0, 1, 1.5) 
mu.2 <- c(2.8, 1.4, 1.8, 0.5, 1.5, 2.5)
mu.3 <- c(0.4, 0.6, 2.4, 2.6, 0.1, 0.4)
basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4, nbasis=6)
B <- getbasismatrix(x, basisBspline, nderiv=0)
Y.1.true <- B %*% mu.1
Y.2.true <- B %*% mu.2
Y.3.true <- B %*% mu.3

plot(x, Y[1, ], col = "purple", type="l", lty = 2, lwd = 0.5, ylim=c(-0.5, 4), ylab="f(x)", xlab="x")
for (i in 2:50) {
  lines(x, Y[i, ], col = "purple", type="l", lty = 2, lwd = 0.5)
}

for (i in 51:100) {
  lines(x, Y[i, ], col = "blue", type="l", lty = 2, lwd = 0.5)
}

for (i in 101:150) {
  lines(x, Y[i, ], col = "green", type="l", lty = 2, lwd = 0.5)
}

lines(x, Y.1.true, col = "black", lwd = 3)
lines(x, Y.2.true, col = "black", lwd = 3)
lines(x, Y.3.true, col = "black", lwd = 3)



Y.1.est <- B %*% t(ext.model.test.4.result$m_list[[1]])
Y.2.est <- B %*% t(ext.model.test.4.result$m_list[[2]])
Y.3.est <- B %*% t(ext.model.test.4.result$m_list[[3]])

lines(x, Y.1.est, col = "red", lwd = 3)
lines(x, Y.2.est, col = "red", lwd = 3)
lines(x, Y.3.est, col = "red", lwd = 3)
legend(x= "topleft", cex=.7, legend=c("True Mean Curves", "VB Estimated Mean Curves"), col=c("black", "red"), lty=c(1,1))
