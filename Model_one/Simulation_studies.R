# Scenario 1 in the paper, Case 1 in codes
# Initialization
number_of_simulations = 1
save_path = NULL

# Data Parameters
x = seq(from=0,to=pi/3, length = 100)
K = 3
curves_per_cluster = 50 
true_cluster_assignments = rep(1:K,each = curves_per_cluster)
seeds = 10

# Pack into data parameter list
data_params = list()
data_params$x = x
data_params$K = K
data_params$curves_per_cluster = curves_per_cluster
data_params$true_cluster_assignments = true_cluster_assignments
data_params$seeds = seeds 
data_params$generate_data = Case_1

# Model Parameters
init = "km"
nbasis = 6
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 30
gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5) * 30
d_not_vector = c(1/3, 1/3, 1/3)
m_not_vector = matrix(c(0.30, 0.41, 0.63, 1.13, 1.68, 2.04,
                        1.00, 1.12, 1.36, 1.88, 2.44, 2.80,
                        0.20, 0.24, 0.31, 0.62, 1.10, 1.44), nrow = 3, byrow = T)
v_not_vector = c(50, 50, 50, 50, 50, 50)
convergence_threshold = 0.01
max_iterations = 1000
verbose = FALSE
draw = TRUE

# Pack into model parameter list 
model_params = list()
model_params$model_func = get_funclustVI_cluster_assignments
model_params$init = "km"
model_params$nbasis = 6
model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
model_params$d_not_vector = d_not_vector
model_params$m_not_vector = m_not_vector 
model_params$v_not_vector = v_not_vector   
model_params$convergence_threshold = convergence_threshold
model_params$max_iterations = max_iterations 
model_params$save_path = save_path
model_params$verbose = verbose
model_params$draw = draw
plot_params = list()
plot_params$xlim = NULL
plot_params$ylim = c(-1, 4)
plot_params$show_curves = TRUE
plot_params$title = 'Scenario 1'
model_params$plot_params = plot_params

# Evaluation parameter list 
eval_func_list = list()
eval_func_list$mismatch = get_mismatches
eval_func_list$vmeasure = get_v_measure

scenarios.1 <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)
####---------------------------------------------------------------------------------------------
# Scenario 2 in the paper, Case 2 in codes
# Initialization
number_of_simulations = 1
save_path = NULL

# Data Parameters
x = seq(from=0,to=pi/3, length = 100)
K = 3
curves_per_cluster = 50 
true_cluster_assignments = rep(1:K,each = curves_per_cluster)
seeds = 10

# Pack into data parameter list
data_params = list()
data_params$x = x
data_params$K = K
data_params$curves_per_cluster = curves_per_cluster
data_params$true_cluster_assignments = true_cluster_assignments
data_params$seeds = seeds 
data_params$generate_data = Case_2

# Model Parameters
init = "km"
nbasis = 6
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 30
gamma_dist_config_matrix[2, ] = c(7.031, 7.031, 7.031) * 30
d_not_vector = c(1/3, 1/3, 1/3)
m_not_vector = matrix(c(0.56, 0.63, 0.80, 0.91, 0.77, 0.61,
                        0.59, 0.68, 0.92, 1.25, 1.37, 1.40,
                        0.67, 0.78, 1.08, 1.56, 1.88, 2.06), nrow = 3, byrow = T)
v_not_vector = c(10, 10, 10, 10, 10, 10)
convergence_threshold = 0.01
max_iterations = 1000
verbose = FALSE
draw = FALSE

# Pack into model parameter list 
model_params = list()
model_params$model_func = get_funclustVI_cluster_assignments
model_params$init = "km"
model_params$nbasis = 6
model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
model_params$d_not_vector = d_not_vector
model_params$m_not_vector = m_not_vector 
model_params$v_not_vector = v_not_vector   
model_params$convergence_threshold = convergence_threshold
model_params$max_iterations = max_iterations 
model_params$save_path = save_path
model_params$verbose = verbose
model_params$draw = draw
plot_params = list()
plot_params$xlim = NULL
plot_params$ylim = c(-0.5, 3)
plot_params$show_curves = FALSE
plot_params$title = 'Scenario 2'
model_params$plot_params = plot_params

# Evaluation parameter list 
eval_func_list = list()
eval_func_list$mismatch = get_mismatches
eval_func_list$vmeasure = get_v_measure

scenarios.2 <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)

####---------------------------------------------------------------------------------------------
# Scenario 3 in the paper, Case 3 in codes
# Initialization
number_of_simulations = 1
save_path = NULL

# Data Parameters
x <- seq(from=0, to=1, by=0.01)
K = 3
curves_per_cluster = 50 
true_cluster_assignments = rep(1:K,each = curves_per_cluster)
seeds = 1

# Pack into data parameter list
data_params = list()
data_params$x = x
data_params$K = K
data_params$curves_per_cluster = curves_per_cluster
data_params$true_cluster_assignments = true_cluster_assignments
data_params$seeds = seeds 
data_params$generate_data = Case_3

# Model Parameters
init = "km"
nbasis = 6
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 10
gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5) * 10
d_not_vector = c(1/3, 1/3, 1/3)
m_not_vector = matrix(c(1.5, 1, 1.8, 2.0, 1, 1.5,
                        2.8, 1.4, 1.8, 0.5, 1.5, 2.5,
                        0.4, 0.6, 2.4, 2.6, 0.1, 0.4), nrow = 3, byrow = T)
v_not_vector = c(10, 10, 10, 10, 10, 10)
convergence_threshold = 0.01
max_iterations = 100
verbose = FALSE
draw = FALSE

# Pack into model parameter list 
model_params = list()
model_params$model_func = get_funclustVI_cluster_assignments
model_params$init = "km"
model_params$nbasis = 6
model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
model_params$d_not_vector = d_not_vector
model_params$m_not_vector = m_not_vector 
model_params$v_not_vector = v_not_vector  
model_params$convergence_threshold = convergence_threshold
model_params$max_iterations = max_iterations 
model_params$save_path = save_path
model_params$verbose = verbose
model_params$draw = draw
plot_params = list()
plot_params$xlim = NULL
plot_params$ylim = c(-1, 4)
plot_params$show_curves = FALSE
plot_params$title = 'Scenario 3'
model_params$plot_params = plot_params

# Evaluation parameter list 
eval_func_list = list()
eval_func_list$mismatch = get_mismatches
eval_func_list$vmeasure = get_v_measure

# sensitivity analysis
scenarios.3 <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)

####---------------------------------------------------------------------------------------------
# Scenario 4 in the paper, Case 4 in codes
# Initialization
number_of_simulations = 1
save_path = NULL

# Data Parameters
x <- seq(from=0, to=1, by=0.01)
K = 3
curves_per_cluster = 50 
true_cluster_assignments = rep(1:K,each = curves_per_cluster)
seeds = 1

# Pack into data parameter list
data_params = list()
data_params$x = x
data_params$K = K
data_params$curves_per_cluster = curves_per_cluster
data_params$true_cluster_assignments = true_cluster_assignments
data_params$seeds = seeds 
data_params$generate_data = Case_4

# Model Parameters
init = "km"
nbasis = 6
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125) * 10
gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5) * 10
d_not_vector = c(1/3, 1/3, 1/3)
m_not_vector = matrix(c(1.5, 1, 1.6, 1.8, 1, 1.5,
                        1.8, 0.6, 0.4, 2.6, 2.8, 1.6,
                        1.2, 1.8, 2.2, 0.8, 0.6, 1.8), nrow = 3, byrow = T)
v_not_vector = c(10, 10, 10, 10, 10, 10)
convergence_threshold = 0.01
max_iterations = 100
verbose = FALSE
draw = TRUE

# Pack into model parameter list 
model_params = list()
model_params$model_func = get_funclustVI_cluster_assignments
model_params$init = "km"
model_params$nbasis = 6
model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
model_params$d_not_vector = d_not_vector
model_params$m_not_vector = m_not_vector 
model_params$v_not_vector = v_not_vector  
model_params$convergence_threshold = convergence_threshold
model_params$max_iterations = max_iterations 
model_params$save_path = save_path
model_params$verbose = verbose
model_params$draw = draw
plot_params = list()
plot_params$xlim = NULL
plot_params$ylim = c(-1, 4)
plot_params$show_curves = TRUE
plot_params$title = 'Scenario 4'
model_params$plot_params = plot_params

# Evaluation parameter list 
eval_func_list = list()
eval_func_list$mismatch = get_mismatches
eval_func_list$vmeasure = get_v_measure

scenarios.4 <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path) 

####---------------------------------------------------------------------------------------------
# Scenario 5, energy consumption, Case 5 in codes
# Initialization
number_of_simulations = 1
save_path = NULL

# Data Parameters
x <- seq(from = 0, to = 24, length.out = 96)
K = 3
curves_per_cluster = 50 
true_cluster_assignments = rep(1:K,each = curves_per_cluster)
seeds = 1

# Pack into data parameter list
data_params = list()
data_params$x = x
data_params$K = K
data_params$curves_per_cluster = curves_per_cluster
data_params$true_cluster_assignments = true_cluster_assignments
data_params$seeds = seeds 
data_params$generate_data = Case_5

# Model Parameters
init = "km"
nbasis = 12
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(347.22, 347.22, 347.22) 
gamma_dist_config_matrix[2, ] = c(0.05, 0.05, 0.05) 
d_not_vector = c(1/3, 1/3, 1/3)
m_not_vector = matrix(c( 0.03, 0.07, -0.03, 0.19, 0.07, 0.05, 0.07, 0.03, 0.12, 0.05, 0.04, 0.04,
                         0.02, 0.01, 0.03, 0.17, -0.01, 0.03, 0.01, 0.03, 0.05, 0.01, 0.02, 0.02,
                         0.03, 0.03, 0.18, 0.02, 0.02, 0.02, 0.02, 0.06, 0.02, 0.02, 0.02, 0.02), nrow = 3, byrow = T)
v_not_vector = rep(10, 12)
convergence_threshold = 0.01
max_iterations = 100
verbose = FALSE
draw = TRUE

# Pack into model parameter list 
model_params = list()
model_params$model_func = get_funclustVI_cluster_assignments
model_params$init = "km"
model_params$nbasis = 12
model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
model_params$d_not_vector = d_not_vector
model_params$m_not_vector = m_not_vector 
model_params$v_not_vector = v_not_vector
model_params$convergence_threshold = convergence_threshold
model_params$max_iterations = max_iterations 
model_params$save_path = save_path
model_params$verbose = verbose
model_params$draw = draw
plot_params = list()
plot_params$xlim = NULL
plot_params$ylim =c(-0.05, 0.20)
plot_params$show_curves = TRUE
plot_params$title = 'Scenario 5'
model_params$plot_params = plot_params

# Evaluation parameter list 
eval_func_list = list()
eval_func_list$mismatch = get_mismatches
eval_func_list$vmeasure = get_v_measure

scenarios.5 <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)

####---------------------------------------------------------------------------------------------
# Scenario 6 in the paper, Case 6 in codes
# Initialization
number_of_simulations = 1
save_path = NULL

#Data Parameters
x = seq(from=0,to=pi/3, length = 100)
K = 4
curves_per_cluster = 50 
true_cluster_assignments = rep(1:K,each = curves_per_cluster)
seeds = 50

# Pack into data parameter list
data_params = list()
data_params$x = x
data_params$K = K
data_params$curves_per_cluster = curves_per_cluster
data_params$true_cluster_assignments = true_cluster_assignments
data_params$seeds = seeds 
data_params$generate_data = Case_6

# Model Parameters
init = "km"
nbasis = 6
gamma_dist_config_matrix = matrix(0, 2, K)
gamma_dist_config_matrix[1, ] = c(78.125, 78.125, 78.125, 78.125) * 25
gamma_dist_config_matrix[2, ] = c(12.5, 12.5, 12.5, 12.5) * 25
d_not_vector = c(1/4, 1/4, 1/4, 1/4)
m_not_vector = matrix(c(0.06,   -0.34,   -1.13,   -0.54,    0.93,    1.66,
                        0.69,    0.18,   -0.78,    0.81,    2.44,    2.82,
                        0.64,    0.03,   -0.96,    1.43,    2.65,    2.62,
                        1.57,    0.83,   -0.08,    3.05,     3.4,    3.03), nrow = 4, byrow = T)
v_not_vector = c(10, 10, 10, 10, 10, 10)
convergence_threshold = 0.01
max_iterations = 100
verbose = FALSE
draw = TRUE

# Pack into model parameter list 
model_params = list()
model_params$model_func = get_funclustVI_cluster_assignments
model_params$init = "km"
model_params$nbasis = 6
model_params$gamma_dist_config_matrix = gamma_dist_config_matrix
model_params$d_not_vector = d_not_vector
model_params$m_not_vector = m_not_vector 
model_params$v_not_vector = v_not_vector
model_params$convergence_threshold = convergence_threshold
model_params$max_iterations = max_iterations 
model_params$save_path = save_path
model_params$verbose = verbose
model_params$draw = draw
plot_params = list()
plot_params$xlim = NULL
plot_params$ylim = c(-2, 4)
plot_params$show_curves = TRUE
plot_params$title = 'Scenario 6'
model_params$plot_params = plot_params

#Evaluation parameter list 
eval_func_list = list()
eval_func_list$mismatch = get_mismatches
eval_func_list$vmeasure = get_v_measure

scenarios.6 <- simulate(data_params, model_params, eval_func_list, number_of_simulations, save_path)


