### real data analysis: Growth-------------------------------------------------------------------------------------------
# following we provide an example of real data analysis of Growth data
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
convergence_threshold = 0.01
max_iterations = 100
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(100, 100) * 20
gamma_dist_config_matrix[2, ] = c(5, 5) * 20
d_not_vector = c(1/3, 2/3)
m_not_vector = matrix(c(70,   82,   85,  122,  141,   148,  177,  180,  181,   181,
                        63,   78,   83,  118,  135,   140,  150,  158,  158,   158), 
                      nrow = 2, byrow = T)
v_not_vector = rep(10, 10)
verbose = FALSE
draw = TRUE
plot_params = list()
plot_params$xlim = NULL
plot_params$ylim = c(66, 196)
plot_params$show_curves = TRUE
plot_params$title = "Growth"

# Fit the model 
growth.clust = funcslustVI(x, Y, K, true_cluster_assignments, 
                    init, nbasis, convergence_threshold, max_iterations, 
                    gamma_dist_config_matrix, d_not_vector, m_not_vector, v_not_vector,
                    verbose, draw, plot_params)

Mismatch(growth.clust$cluster_assignments, true_cluster_assignments, 2)
sabre::vmeasure(growth.clust$cluster_assignments, true_cluster_assignments)$v_measure

# repetitions
growth.model <- function(seed){
  set.seed(seed)
  model = funcslustVI(x, Y, K, true_cluster_assignments, 
                      init, nbasis, convergence_threshold, max_iterations, 
                      gamma_dist_config_matrix, d_not_vector, m_not_vector, v_not_vector, 
                      verbose, draw, plot_params)
  
  mis.mat <- Mismatch(model$cluster_assignments, true_cluster_assignments, 2)
  v.measure <- sabre::vmeasure(model$cluster_assignments, true_cluster_assignments)$v_measure
  ELBO_converged <- model$`Converged ELBO`
  return(c(mis.mat, v.measure, ELBO_converged))
}


growth.analysis <- function(repetition.times){
  start.time <- Sys.time()
  growth.metrics <- sapply(1:repetition.times, growth.model)
  max_ELBO_seed <- which.max(growth.metrics[3, ])
  mis.mean <- mean(growth.metrics[1, ])
  v.mean <- mean(growth.metrics[2, ])
  mis.sd <- sd(growth.metrics[1, ])
  v.sd <- sd(growth.metrics[2, ])
  end.time <- Sys.time()
  time.length <- end.time - start.time
  return(list(mismatch.mean = mis.mean,
              vmeasure.mean = v.mean,
              mismatch.sd = mis.sd,
              vmeasure.sd = v.sd,
              runing.time = time.length,
              max.ELBO.seed = max_ELBO_seed))
}
growth.performance <- growth.analysis(50)

# plot of Growth
cluster_number = true_cluster_assignments[1]
col = 2 + cluster_number
plot(x, Y[1, ], col=col, type="l", ylim=c(66, 196), xlim=c(1, 31), main="Growth", ylab="f(x)", xlab="x")
for(i in 2:NROW(Y)) {
  cluster_number = true_cluster_assignments[i]
  col = 2 + cluster_number
  lines(x, Y[i, ], type="l", col=col)
}

