#' Update the sigma parameter for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param K The total number of clusters
#' @param A_vector A vector
#' @param R_vector A vector
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param x The x used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @return A list of the updated sigma parameters for each cluster

update_sigma_list <- function(Y, phi_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, v_not_vector) {
  ev_q_tau = A_vector / R_vector
  I = diag(nbasis)
  if(is.null(v_not_vector)){
    v_not_vector = get_v_not_vector(Y, x, probability_matrix, phi_matrix, K, nbasis)
    warning("No prior information for v_not_vector, approximated values have been used based on your data!")
  }
  B = get_B(x, nbasis)
  sigma_list = list()
  for (i in 1:K) {
    sum_matrix = matrix(0, nbasis, nbasis)
    for (j in 1:NROW(Y)) {
      p = probability_matrix[j, i]
      # v_not = v_not_vector[i]
      temp = p*(t(B) %*% B) ### check, change
      sum_matrix = sum_matrix + temp
    }
    
    mat = ev_q_tau[i] * sum_matrix + v_not_vector[i] * I # here change
    sigma = MASS::ginv(mat)
    sigma_list[[i]] = sigma
  }
  return(sigma_list)
}

#' Update the m parameter for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param m_not_vector The vector containing m_not values for each cluster
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param K The total number of clusters
#' @param A_vector A vector
#' @param R_vector A vector
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param x The x used to generate the clusters
#' @param nbasis The number of basis functions
#' @param sigma_list A list of the updated sigma parameters for each cluster
#'
#' @return A list of the updated m parameters for each cluster

update_m_list <- function(Y, m_not_vector, phi_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, sigma_list, v_not_vector) {
  ev_q_tau = A_vector / R_vector
  if(is.null(v_not_vector)){
    v_not_vector = get_v_not_vector(Y, x, probability_matrix, phi_matrix, K, nbasis)
    warning("No prior information for v_not_vector, approximated values have been used based on your data!")
  }
  B = get_B(x, nbasis)
  m_list = list()
  
  for (i in 1:K) {
    sum_matrix = matrix(0, 1, nbasis)
    for (j in 1:NROW(Y)) {
      p = probability_matrix[j, i]
      temp = p * (Y[j, ] %*% B)
      sum_matrix = sum_matrix + temp
    }
    
    m = (ev_q_tau[i] * sum_matrix + v_not_vector[i] * m_not_vector[i, ]) %*% sigma_list[[i]] 
    m_list[[i]] = m
  }
  return(m_list)
}

#' Update the d parameter for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param x The x used to generate the clusters
#' @param sigma_list A list of the updated sigma_list parameters for each cluster
#' @param m_list A vlist of the updated m parameters for each cluster
#' @param nbasis The number of basis functions]
#' @param beta_vector A vector of the beta parameters for each cluster
#' 
#' @return A vector of the updated R parameters for each cluster

update_R_vector <- function(Y, K, probability_matrix, x, sigma_list, m_list, nbasis, beta_vector) {
  
  B = get_B(x, nbasis)
  
  R_vector = c(1:K)
  for (i in 1:K) {
    r_not = beta_vector[i]
    sum = 0
    for (j in 1:NROW(Y)) {
      p = probability_matrix[j, i]
      ev_phi = sum(diag(B %*% sigma_list[[i]] %*% t(B))) + t((Y[j, ] - B %*% t(m_list[[i]]))) %*% (Y[j, ] - B %*% t(m_list[[i]]))
      sum = sum + (p * ev_phi)
    }
    R = r_not + 1/2 * sum
    R_vector[i] = R
  }
  return(R_vector)
}

#' Update the d parameter for each cluster
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' 
#' @return A vector of the updated d parameters for each cluster

update_d_vector <- function(Y, K, probability_matrix, d_not_vector) {
  if(is.null(d_not_vector)){
    d_not_vector = get_d_not_vector(K)
    warning("No prior information for d_not, 1/K has been used for each cluster.")
  }
  
  d_vector = c(1:K)*0
  for (i in 1:K) {
    sum = 0
    for (j in 1:NROW(Y)) {
      p = probability_matrix[j, i]
      sum = sum + p
    }
    d = d_not_vector[i] + sum
    d_vector[i] = d
  }
  return(d_vector)
}

#' Update the probabilty matrix
#'
#' @param Y The matrix containing rows corresponding the curves
#' @param K The total number of clusters
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param x The x used to generate the clusters
#' @param m_list A list of the m parameters for each cluster
#' @param A_vector A vector
#' @param R_vector A vector
#' @param d_vector A vector of the d parameters for each cluster
#' @param nbasis The number of basis functions
#' 
#' @return The updated probability matrix 

update_probability_matrix <- function(Y, K, x, sigma_list, m_list, A_vector, R_vector, d_vector, nbasis) {
  probability_matrix = matrix(0, NROW(Y), K)
  observations_per_curve = NCOL(Y)
  B = get_B(x, nbasis)
  
  coef = observations_per_curve / 2
  for (k in 1:K) {
    for (i in 1:NROW(Y)) {
      ev_phi = sum(diag(B %*% sigma_list[[k]] %*% t(B))) + t((Y[i, ] - B %*% t(m_list[[k]]))) %*% (Y[i, ] - B %*% t(m_list[[k]]))
      ev_tau_k_log_tau_k = digamma(A_vector[k]) - log(R_vector[k])
      ev_tau_k_tau_k = A_vector[k] / R_vector[k]
      ev_pi_log_pi_k = digamma(d_vector[k]) - digamma(sum((d_vector)))
      
      a_i_k = coef * ev_tau_k_log_tau_k - 1/2 * ev_tau_k_tau_k * ev_phi + ev_pi_log_pi_k
      num = exp(a_i_k)
      
      probability_matrix[i, k] = num
    }
  }
  probability_matrix_row_sum = apply(probability_matrix, 1, sum)
  for (i in 1:nrow(probability_matrix)) {
    if(probability_matrix_row_sum[i] == 0) {
      maxindex <- which.max(probability_matrix[i, ])
      probability_matrix[i, maxindex] <- 1
      probability_matrix[i, -maxindex] <- 0
    }
    else probability_matrix[i, ] <- probability_matrix[i, ] / probability_matrix_row_sum[i]
  }
  return(probability_matrix)
}