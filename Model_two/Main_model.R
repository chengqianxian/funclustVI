library(fda)
# construct a function for the model
funcslustVI.ext <- function(x, Y, K, true_cluster_assignments, init, nbasis, 
                            convergence_threshold, max_iterations, 
                            b_not = NULL, r_not = NULL, d_not_vector = NULL, m_not_matrix = NULL, v_not_vector = NULL,
                            alpha_not = NULL, beta_not = NULL) {
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
  
  B = get_B(x, nbasis)
  # initialization
  R_vector = rep(r_not, K)
  mu_a_vector = rep(0, nrow(Y))
  beta = beta_not
  # constant
  alpha <- alpha_not + nrow(Y) / 2
  converged = FALSE
  prev_elbo = NULL
  iteration = 0
  #curr_m_list <- NULL
  while (converged == FALSE & iteration <= max_iterations) {
    iteration = iteration + 1
    A_vector <- update_A_k_vector(b_not, Y, K, probability_matrix)
    sigma_list = update_sigma_list(Y, K, A_vector, R_vector, probability_matrix, x, nbasis, v_not_vector)
    m_list = update_m_list(Y, m_not_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, sigma_list, v_not_vector, mu_a_vector)
    sigma_a_vector = update_sigma_a_vector(Y, K, A_vector, R_vector, probability_matrix, alpha, beta)
    mu_a_vector = update_mu_a_vector(Y, x, nbasis, K, m_list, sigma_a_vector, probability_matrix, A_vector, R_vector)
    R_vector = update_R_vector(Y, K, probability_matrix, x, sigma_list, m_list, nbasis, r_not, sigma_a_vector, mu_a_vector)
    beta = update_beta(beta_not, sigma_a_vector, mu_a_vector)
    d_vector = update_d_vector(Y, K, probability_matrix, d_not_vector)
    probability_matrix = update_probability_matrix(Y, K, x, sigma_list, m_list, A_vector, R_vector, d_vector, nbasis, sigma_a_vector, mu_a_vector)
    elbo = get_elbo(x, Y, K, d_not_vector, v_not_vector, b_not, r_not, m_not_matrix, alpha_not, beta_not, 
                    nbasis, sigma_list, m_list, A_vector, R_vector, d_vector, sigma_a_vector, mu_a_vector, 
                    alpha, beta, probability_matrix)
    curr_elbo = elbo[1] 
    converged = check_convergence(prev_elbo, curr_elbo, convergence_threshold)
    prev_elbo = curr_elbo
  }
  cluster_assignments <- get_final_cluster_assignments(probability_matrix)
  return(list(iteration = iteration,
              cluster_assignments = cluster_assignments,
              A_vector = A_vector,
              sigma_list = sigma_list,
              m_list = m_list,
              sigma_a_vector = sigma_a_vector,
              mu_a_vector = mu_a_vector,
              R_vector = R_vector,
              beta = beta,
              d_vector = d_vector,
              probability_matrix = probability_matrix,
              elbo = elbo))
}

check_convergence.2 <- function(m_list, K, curr_m_list, convergence_threshold){
  if(is.null(curr_m_list)) return(FALSE)
  else{
    m.list.sum <- numeric(K)
    for (i in 1:K) {
      m_list.diff <- numeric(K)
      for (j in 1:K) {
        m_list.diff[j] <- mean(abs(m_list[[i]] - curr_m_list[[j]]))
      }
      m.list.sum[i] <- max(m_list.diff)
    }
    if(sum(m.list.sum) <= convergence_threshold) return(TRUE)
    else return(FALSE)
  }
}


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

get_approx_probability_matrix_hcl <- function(Y, K, x) {
  predictions = cutree(hclust(dist(Y), method = "ward"), K) # use distance for clustering
  probability_matrix = matrix(0, NROW(Y), K)
  for (i in 1:length(predictions)) {
    cluster_prediction = predictions[[i]]
    probability_matrix[i, cluster_prediction] = 1
  }
  
  return(probability_matrix)
}

get_true_probability_matrix <- function(true_cluster_assignments, K) {
  number_of_curves = length(true_cluster_assignments)
  probability_matrix = matrix(0, number_of_curves, K)
  for (i in 1:number_of_curves) {
    col = true_cluster_assignments[i]
    probability_matrix[i, col] = 1
  }
  return(probability_matrix)
}

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

get_B <- function(x, nbasis) {
  rangeval = c(0, x[length(x)])
  basisobj = fda::create.bspline.basis(rangeval, nbasis)
  B <- fda::getbasismatrix(x, basisobj=basisobj)
  return(B)
}

