### real data analysis: Growth-------------------------------------------------------------------------------------------
## following is an example of real data analysis using Model 2
library(fda)
boy.ht <- t(as.matrix(growth$hgtm))
girl.ht <- t(as.matrix(growth$hgtf))
ht.data <- rbind(boy.ht, girl.ht)
colnames(ht.data) <- NULL
rownames(ht.data) <- NULL

# Data Parameters
x = seq(1:31)
Y = ht.data 
K = 2
curves_per_cluster = c(39, 54)
true_cluster_assignments = rep(1:K, curves_per_cluster)

# Model Parameters 
init = "km"
nbasis = 10 
convergence_threshold = 0.001
max_iterations = 1000
b_not = 100 *10
r_not = 5 *10
d_not_vector = c(1/3, 2/3)
m_not_matrix = matrix(c(70,   82,   85,  122,  141,   148,  177,  180,  181,   181,
                        63,   78,   83,  118,  135,   140,  150,  158,  158,   158), 
                      nrow = 2, byrow = T)
v_not_vector = rep(10, nbasis)
alpha_not = 100 * 10
beta_not = 10 * 10

# Fit the model 
growth.clust <- funcslustVI.ext(x, Y, K, true_cluster_assignments, init, nbasis, 
                               convergence_threshold, max_iterations, 
                               b_not, r_not, d_not_vector, m_not_matrix, v_not_vector,
                               alpha_not, beta_not)


Mismatch(growth.clust$cluster_assignments, true_cluster_assignments, 2)
sabre::vmeasure(growth.clust$cluster_assignments, true_cluster_assignments)$v_measure


# plot
plot(x, Y[1, ], col = "green", type="l", lty = 2, lwd = 0.5, ylim=c(66, 196), 
     ylab="f(x)", xlab="x", main = "Growth")
for (i in 2:39) {
  lines(x, Y[i, ], col = "green", type="l", lty = 2, lwd = 0.5)
}

for (i in 40:93) {
  lines(x, Y[i, ], col = "blue", type="l", lty = 2, lwd = 0.5)
}

basisBspline <- create.bspline.basis(c(min(x),max(x)),norder=4,nbasis=10)
B <- getbasismatrix(x, basisBspline, nderiv=0)
Y.1.est <- B %*% t(growth.clust$m_list[[1]])
Y.2.est <- B %*% t(growth.clust$m_list[[2]])

lines(x, Y.1.est, col = "red", lwd = 3)
lines(x, Y.2.est, col = "red", lwd = 3)

h.boy <- Y[1:39, ]
h.mean.boy <- apply(h.boy, 2, mean)
h.girl <- Y[40:93, ]
h.mean.girl <- apply(h.girl, 2, mean)
lines(x, h.mean.boy, col = "black", type="l", lty = 1, lwd = 3)
lines(x, h.mean.girl, col = "black", type="l", lty = 1, lwd = 3)
legend(x= "topleft", cex=.7, legend=c("Empirical Mean Curves", "VB Estimated Mean Curves", "Raw Curves"), 
       col=c("black", "red", "grey"), lty=c(1, 1, 2))
