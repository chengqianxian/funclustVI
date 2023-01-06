library(fda)
library(MASS)

#' Generates cluster assignments and related information given functional data. 
#'
#' @param Y A matrix in which the rows represent the curves 
#' @param K The number of clusters in the data
#' @param nbasis The number of basis functions
#' @param x The dependent variable 
#' @param init The initialization method for the algorithim
#' @param true_cluster_assignments The true cluster assignments 
#' @param gamma_dist_config_matrix A matrix where the rows are the a_not and r_not for each cluster
#' @param d_not_vector A vector of d which is the paramter of Dirichlet distribution for pi
#' @param m_not_vector A matrix of mean vector of phi, the coefficient vector,for all clusters
#' @param v_not_vector A vector of the precision of phi
#' @param convergence_threshold The threshold that determines when the model has converged 
#' @param max_iterations The maximum amount of iterations for the algorithim
#' @param verbose A boolean indicating whether or not to print the inner parameters of the model
#' @param plot_params List of parameters corresponding to characteristics of the plot. Must include vectors xlim and ylim corresponding to the x and y limits of the plot. 
#' @param draw Parameter that determines whether or not to plot the real vs estimated curves after fitting the model 
#' 
#' @return A list with entries containing various information about the fitted model

funcslustVI <- function(x, Y, K, true_cluster_assignments, init, nbasis, 
                        convergence_threshold = 0.01, max_iterations = 1000, 
                        gamma_dist_config_matrix = NULL, d_not_vector = NULL, m_not_vector = NULL, v_not_vector = NULL,
                        verbose, draw, plot_params) {
  probability_matrix = NULL
  
  if (init == 'hcl') {
    probability_matrix = get_approx_probability_matrix_hcl(Y, K)
  } else if (init == 'tpm') {
    probability_matrix = get_true_probability_matrix(true_cluster_assignments, K)
  } else if (init == "cust") {
    probability_matrix = init
  } else {
    probability_matrix = get_approx_probability_matrix_km(Y, K, x)
  }
  
  init_cluster_assignments = get_final_cluster_assignments(probability_matrix) 
  
  tau_list = get_tau_list(Y, probability_matrix, K) # this can be moved to line 63, else....
  
  if (is.null(true_cluster_assignments) == FALSE) {
    true_m_not = get_true_m_not(x, Y, K, nbasis, true_cluster_assignments)
  }
  if(is.null(m_not_vector)){
    phi_matrix = get_approx_phi_matrix(Y, K, nbasis, probability_matrix, x)
    m_not_vector = phi_matrix
    warning("No prior information for m_not, approximated values have been used based on your data!")
  }
  A_vector = NULL
  R_vector = NULL
  alpha_vector = NULL
  beta_vector = NULL
  if (is.null(gamma_dist_config_matrix) != TRUE) {
    alpha_vector = gamma_dist_config_matrix[1, ]
    R_vector = c(gamma_dist_config_matrix[2, ])
    beta_vector = c(gamma_dist_config_matrix[2, ])
  } else {
    alpha_vector = get_alpha_vector(tau_list)
    beta_vector = get_beta_vector(tau_list)
    R_vector = beta_vector
    warning("No prior information for a_not and r_not, approximated values have been used based on your data!")
  }
  
  if (verbose == TRUE) {
    print("Probability Matrix Initialization")
    print(probability_matrix)
  }
  
  B = get_B(x, nbasis)
  
  converged = FALSE
  prev_elbo = 0
  iteration = 0
  while (converged == FALSE & iteration <= max_iterations) {
    iteration = iteration + 1 
    A_vector = update_A_k_vector(alpha_vector, Y, K, probability_matrix)
    sigma_list = update_sigma_list(Y, phi_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, v_not_vector)
    m_list = update_m_list(Y, m_not_vector, phi_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, sigma_list, v_not_vector)
    R_vector = update_R_vector(Y, K, probability_matrix, x, sigma_list, m_list, nbasis, beta_vector)
    d_vector = update_d_vector(Y, K, probability_matrix, d_not_vector)
    probability_matrix = update_probability_matrix(Y, K, x, sigma_list, m_list, A_vector, R_vector, d_vector, nbasis)
    phi_matrix = get_approx_phi_matrix(Y, K, nbasis, probability_matrix, x)
    elbo = get_elbo(x, Y, K, phi_matrix, m_not_vector, nbasis, sigma_list, m_list, A_vector, R_vector, d_vector, probability_matrix, alpha_vector, beta_vector)
    curr_elbo = elbo[1] 
    converged = check_convergence(prev_elbo, curr_elbo, convergence_threshold)
    prev_elbo = curr_elbo
    
    if (verbose == TRUE) {
      cat("Iteration: ", toString(iteration))
      print("Sigma List")
      print(sigma_list)
      print("M List")
      print(m_list)
      print("R Vector")
      print(R_vector)
      print("D vector")
      print(d_vector)
      print("Probability Matrix")
      print(probability_matrix)
      cat("ELBO: ", toString(curr_elbo))
    }
  }
  # if(itert)
  cluster_assignments = get_final_cluster_assignments(probability_matrix)
  
  if (draw == TRUE & is.null(true_cluster_assignments) == FALSE) {
    plot_data(x, Y, B, m_list, true_m_not, true_cluster_assignments, plot_params)
  }
  result_list = list("probability_matrix" = probability_matrix, "cluster_assignments" = cluster_assignments, 
                     "init_cluster_assignments" = init_cluster_assignments,
                     "Total iteration" = iteration,
                     "A vector" = A_vector,
                     "R vector" = R_vector,
                     "Sigma List" = sigma_list,
                     "M list" = m_list,
                     "Converged ELBO" = elbo,
                     "True basis coefficients" = true_m_not)
  return(result_list)
}

#' Generates a matrix in which each row represents the variances of the curves in each cluster
#'
#' @param Y The matrix containing rows corresponding to the curves
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param K The number of entries in each curve vector
#'
#' @return A matrix in which each row represents the variances of the curves in each cluster

get_tau_list <- function(Y, probability_matrix, K) {
  cumulative_matrix_list = vector("list", length = K)
  
  for (i in 1:NROW(Y)) {
    max_prob_index = which(probability_matrix[i,] == max(probability_matrix[i,]))
    if(length(max_prob_index) != 1) {
      max_prob_index = max_prob_index[sample(1:length(max_prob_index), 1)]
    }
    
    if(is.null(cumulative_matrix_list[[max_prob_index]]) == TRUE) {
      mat = matrix(0, 1, NCOL(Y))
      mat[1, ] = Y[i, ]
      cumulative_matrix_list[[max_prob_index]] = mat
    } else {
      cumulative_matrix_list[[max_prob_index]] = rbind(cumulative_matrix_list[[max_prob_index]], Y[i, ])
    }
  }
  # the purpose of this part above is to group the curves into K clusters
  tau_list = vector("list", length = K)
  for (i in 1:K) {
    vec <- vector()
    tau_list[[i]] <- vec
  }
  
  for (i in 1:K) {
    mat = cumulative_matrix_list[[i]]
    
    if (is.null(mat) == TRUE) {
      tau_list[[i]] = c(1:NCOL(Y))*0 
    } else {
      for(j in 1:NCOL(Y)) {
        Var = var(mat[,j])
        if(is.na(Var) == TRUE) Var = 10^(-20)
        tau = 1 / Var 
        tau_list[[i]] = c(tau_list[[i]], tau) 
      }
    }
  }
  
  return(tau_list)
  
}

#' Generates a vector A
#'
#' @param alpha_vector A vector that in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#' @param observations_per_curve The number of entries in each curve vector
#'
#' @return A vector A


update_A_k_vector <- function(alpha_vector, Y, K, probability_matrix){
  A_vector <- numeric(K)
  n <- ncol(Y)
  for (k in 1:K) {
    A_vector[k] <- alpha_vector[k] + n / 2 * sum(probability_matrix[, k])   
  }
  return(A_vector)
}

#' Generates a vector in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @param tau_list A list in which each index is the tau vector for that cluster
#'
#' @return A vector that in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'

get_alpha_vector <- function(tau_list) {
  alpha_vector = c(1:length(tau_list))*0
  for (cluster_number in 1:length(tau_list)) {
    tau_k = tau_list[[cluster_number]]
    if(length(unique(tau_k)) == 1) {
      alpha_vector[cluster_number] = 0 
    } else {
      expected_value = mean(tau_k)
      variance = var(tau_k)
      alpha = expected_value ^ 2 / variance # use Gamma property
      alpha_vector[cluster_number] = alpha
    }
  }
  return(alpha_vector)
}

#' Generates a vector in which the entries are the beta parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @param tau_list A list in which each index is the tau vector for that cluster
#'
#' @return A vector that in which the entries are the beta parameters of the gamma distribution (1 / variance) of the curves in each cluster

get_beta_vector <- function(tau_list) {
  beta_vector = c(1:length(tau_list))*0
  for (cluster_number in 1:length(tau_list)) {
    tau_k = tau_list[[cluster_number]]
    if(length(unique(tau_k)) == 1) {
      beta_vector[cluster_number] = 0
    } else {
      expected_value = mean(tau_k)
      variance = var(tau_k)
      beta = expected_value / variance
      beta_vector[cluster_number] = beta
    }
  }
  return(beta_vector)
}

#' Generates initial probability matrix with approximated probability values via K-mean clustering for each cluster
#' 
#' @param Y The matrix containing rows corresponding the curves
#' @param K The number of clusters in cluster data
#' @param x The dependent variable 
#'
#' @return A probability matrix with approximated probability values for each cluster using kmeans 

get_approx_probability_matrix_km <- function(Y, K, x) {  # In this function, x is not necessary currently
  res = kmeans(Y, K) 
  predictions = res$cluster
  probability_matrix = matrix(0, NROW(Y), K)
  for (i in 1:length(predictions)) {
    cluster_prediction = predictions[[i]]
    probability_matrix[i, cluster_prediction] = 1
  }
  
  return(probability_matrix)
}

#' Generates initial probability matrix with approximated probability values via hierarchical clustering for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The number of clusters in cluster data
#' @param x The x used to generate the clusters
#'
#' @return A probability matrix with probability values for each cluster

get_approx_probability_matrix_hcl <- function(Y, K, x) {
  predictions = cutree(hclust(dist(Y), method = "ward"), K) # use distance for clustering
  probability_matrix = matrix(0, NROW(Y), K)
  for (i in 1:length(predictions)) {
    cluster_prediction = predictions[[i]]
    probability_matrix[i, cluster_prediction] = 1
  }
  
  return(probability_matrix)
}

#' Generates B matrix
#'
#' @param x The x used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @return A matrix in which row represent curves of various clusters


get_B <- function(x, nbasis) {
  rangeval = c(0, x[length(x)])
  basisobj = fda::create.bspline.basis(rangeval, nbasis)
  B <- fda::getbasismatrix(x, basisobj=basisobj)
  return(B)
}

#' Gets matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param nbasis The number of basis functions
#' @param probability_matrix The x used to generate the clusters
#' @param x The x used to generate the clusters
#' 
#' @return A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B

get_approx_phi_matrix <- function(Y, K, nbasis, probability_matrix, x) {
  function_data = get_approx_function_data(Y, K, probability_matrix)
  observations_per_curve = length(x)
  min_arg = x[1]
  max_arg = x[observations_per_curve]
  basisobj = fda::create.bspline.basis(c(min_arg, max_arg), nbasis)
  phi_matrix = matrix(0, K, nbasis)
  for (i in 1:K) {
    f = function_data[i, ]
    phi = fda::smooth.basis(argvals = x, y = f, fdParobj = basisobj)$fd$coef
    phi_matrix[i, ] = c(phi)
  }
  return(phi_matrix)
}

#' Gets the vnot parameter
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param x The x used to generate the clusters
#' @param probability_matrix The x used to generate the clusters
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param K The total number of clusters
#' @param nbasis The number of basis functions
#'
#' @return The vnot parameter


get_v_not_vector <- function(Y, x, probability_matrix, phi_matrix, K, nbasis) {
  cumulative_phi_matrix = get_cumulative_phi_matrix(Y, x, nbasis)
  cumulative_phi_matrix_list = vector("list", length = K)
  
  for (i in 1:NROW(Y)) {
    max_prob_index = which(probability_matrix[i,] == max(probability_matrix[i,]))
    if(length(max_prob_index) != 1) {
      max_prob_index = max_prob_index[sample(1:length(max_prob_index), 1)]
    }
    
    if (is.null(cumulative_phi_matrix_list[[max_prob_index]]) == TRUE)   {
      mat = matrix(0, 1, nbasis)
      mat[1, ] = cumulative_phi_matrix[i, ]
      cumulative_phi_matrix_list[[max_prob_index]] = mat
    } else {
      cumulative_phi_matrix_list[[max_prob_index]] = rbind(cumulative_phi_matrix_list[[max_prob_index]], cumulative_phi_matrix[i, ])
    }
  }
  
  v_not_vector = c(1:K)
  for(i in 1:K) {
    var_sum = 0
    mat = cumulative_phi_matrix_list[[i]]
    if (is.null(mat) == TRUE) {
      v_not = 0
    } else {
      for(j in 1:nbasis) {
        var_sum = var_sum + var(mat[, j])
      } 
      if(is.na(var_sum) == TRUE) var_sum = 10^(-20)
      avg_var = var_sum / nbasis
      v_not = 1 / avg_var
    }
    v_not_vector[i] = v_not
  }
  return(v_not_vector)
}

#' Gets the d_not parameter for each cluster
#'
#' @param K The total number of clusters
#'
#' @return The d_not parameter vector for each cluster

get_d_not_vector <- function(K) {
  val = 1 / K
  d_not_vector = c(1:K)
  
  for (i in 1:K) {
    d_not_vector[i] = val
  }
  
  return(d_not_vector)
}

#' Gets functional data approximations, use mean for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#'
#' @return Matrix with approximations of functional data for each cluster in a row

get_approx_function_data <- function(Y, K, probability_matrix) {
  sum_count_vector = c(1:K)*0
  sum_matrix = matrix(0, K, NCOL(Y))
  
  for (i in 1:NROW(Y)) {
    max_prob_index = which(probability_matrix[i,] == max(probability_matrix[i,]))
    if(length(max_prob_index) != 1) {
      max_prob_index = max_prob_index[sample(1:length(max_prob_index), 1)]
    }
    
    sum_count_vector[max_prob_index] = sum_count_vector[max_prob_index] + 1
    sum_matrix[max_prob_index, ] = sum_matrix[max_prob_index, ] + Y[i, ]
  }
  
  function_data = matrix(0, K, NCOL(Y))
  for (i in 1:K) {
    function_data[i, ] = sum_matrix[i, ] / sum_count_vector[i]
  }
  
  return(function_data)
}

#' Gets a matrix of the coffecient vectors, phi, for each curve
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param x The x used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @return a matrix of the coffecient vectors, phi, for each curve

get_cumulative_phi_matrix <- function(Y, x, nbasis) {
  observations_per_curve = length(x)
  min_arg = x[1]
  max_arg = x[observations_per_curve]
  basisobj = fda::create.bspline.basis(c(min_arg, max_arg), nbasis)
  cumulative_phi_matrix = matrix(0, NROW(Y), nbasis)
  for (i in 1:NROW(Y)) {
    f = Y[i, ]
    phi = fda::smooth.basis(argvals = x, y = f, fdParobj = basisobj)$fd$coef
    cumulative_phi_matrix[i, ] = c(phi)
  }
  return(cumulative_phi_matrix)
}

#' Gets a vector of the final cluster assignments based on the probability matrix 
#'
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' 
#' @return A vector of the final cluster assignments

get_final_cluster_assignments <- function(probability_matrix) {
  final_cluster_assignments = c(1:NROW(probability_matrix))
  
  for (i in 1:NROW(probability_matrix)) {
    max_prob_index = which(probability_matrix[i,] == max(probability_matrix[i,]))
    if(length(max_prob_index) != 1) {
      max_prob_index = max_prob_index[sample(1:length(max_prob_index), 1)]
    }
    
    final_cluster_assignments[i] = max_prob_index
  }
  return(final_cluster_assignments)
}

#' Gets the true probability matrix using the true cluster assignments 
#'
#' @param true_cluster_assignments A vector containing the true cluster assignments 
#' @param K The number of clusters in the data 
#' 
#' @return The true probability matrix of the true cluster assignments

get_true_probability_matrix <- function(true_cluster_assignments, K) {
  number_of_curves = length(true_cluster_assignments)
  probability_matrix = matrix(0, number_of_curves, K)
  for (i in 1:number_of_curves) {
    col = true_cluster_assignments[i]
    probability_matrix[i, col] = 1
  }
  return(probability_matrix)
}

#' Gets the vector with the true m_not values
#'
#' @param x The x used to generate the clusters
#' @param Y The matrix containing rows corresponding the curves
#' @param K The number of clusters in the data 
#' @param nbasis The number of basis functions
#' @param true_cluster_assignments A vector containing the true cluster assignments 
#' 
#' 
#' @return The vector with the true m_not values 


get_true_m_not <- function(x, Y, K, nbasis, true_cluster_assignments) {
  true_probability_matrix = get_true_probability_matrix(true_cluster_assignments, K)
  true_m_not = get_approx_phi_matrix(Y, K, nbasis, true_probability_matrix, x)
  return(true_m_not)
}

#' @param m_list A list of the updated m parameters for each cluster
#' @param true_m_not A matrix containing the true m_not vectors for each cluster

get_mse <- function(B, true_m_not, m_list){
  number_of_clusters = NROW(true_m_not)
  mse <- numeric(number_of_clusters)
  for (i in 1:number_of_clusters) {
    mse[i] <- mean(((B %*% true_m_not[i, ]) - (B %*% t(m_list[[i]])))^2) 
  }
  return(mse)
}


