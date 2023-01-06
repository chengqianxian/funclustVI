update_sigma_list <- function(Y, K, A_vector, R_vector, probability_matrix, x, nbasis, v_not_vector) {
  ev_q_tau = A_vector / R_vector
  I = diag(nbasis)
  B = get_B(x, nbasis)
  
  sigma_list = list()
  for (i in 1:K) {
    sum_matrix = matrix(0, nbasis, nbasis)
    for (j in 1:nrow(Y)) {
      p = probability_matrix[j, i]
      temp = p*(t(B) %*% B)
      sum_matrix = sum_matrix + temp
    }
    
    mat = ev_q_tau[i] * sum_matrix + v_not_vector[i] * I
    sigma = MASS::ginv(mat)
    sigma_list[[i]] = sigma
  }
  return(sigma_list)
}

update_m_list <- function(Y, m_not_matrix, K, A_vector, R_vector, probability_matrix, x, nbasis, sigma_list, v_not_vector, mu_a_vector) {
  ev_q_tau = A_vector / R_vector
  
  B = get_B(x, nbasis)
  m_list = list()
  Y_star <- Y
  for (i in 1:nrow(Y)) {
    Y_star[i, ] <- Y[i, ] - mu_a_vector[i]
  }
  
  for (i in 1:K) {
    sum_matrix = matrix(0, 1, nbasis)
    for (j in 1:nrow(Y_star)) {
      p = probability_matrix[j, i]
      temp = p * (Y_star[j, ] %*% B)
      sum_matrix = sum_matrix + temp
    }
    
    m = (ev_q_tau[i] * sum_matrix + v_not_vector[i] * m_not_matrix[i, ]) %*% sigma_list[[i]]
    m_list[[i]] = m
  }
  return(m_list)
}

update_sigma_a_vector <- function(Y, K, A_vector, R_vector, probability_matrix, alpha, beta){
  n <- ncol(Y)
  N <- nrow(Y)
  ev_q_tau = A_vector / R_vector
  sigma_a_vector <- numeric(N)
  for (i in 1:N) {
    sigma_a_vector[i] <- 1 / (n * sum(probability_matrix[i, ] * ev_q_tau) + alpha / beta)
  }
  return(sigma_a_vector)
}

update_mu_a_vector <- function(Y, x, nbasis, K, m_list, sigma_a_vector, probability_matrix, A_vector, R_vector){
  B <- get_B(x, nbasis)
  mu_a_vector <- numeric(nrow(Y))
  ev_q_tau = A_vector / R_vector
  N <- nrow(Y)
  for (i in 1:N) {
    mu_a_k <- numeric(K)
    for (k in 1:K) {
      mu_a_k[k] <- ev_q_tau[k] * probability_matrix[i, k] * (sum(Y[i, ] - (B %*% t(m_list[[k]])))) # error here? * to -
    }
    mu_a_vector[i] <- sum(mu_a_k) * sigma_a_vector[i]
  }
  return(mu_a_vector)
}

update_A_k_vector <- function(b_not, Y, K, probability_matrix){
  A_vector <- numeric(K)
  n <- ncol(Y)
  for (k in 1:K) {
    A_vector[k] <- b_not + n / 2 * sum(probability_matrix[, k])   
  }
  return(A_vector)
}
update_R_vector <- function(Y, K, probability_matrix, x, sigma_list, m_list, nbasis, r_not, sigma_a_vector, mu_a_vector) {
  B <- get_B(x, nbasis)
  R_vector <- numeric(K)
  n <- ncol(Y)
  for (k in 1:K) {
    sum_R = 0
    for (i in 1:nrow(Y)) {
      p = probability_matrix[i, k]
      ev_phi_a = sum(diag(B %*% sigma_list[[k]] %*% t(B))) + n * sigma_a_vector[i] + 
        t((Y[i, ] - B %*% t(m_list[[k]]) - mu_a_vector[i])) %*% (Y[i, ] - B %*% t(m_list[[k]]) - mu_a_vector[i])
      sum_R = sum_R + (p * ev_phi_a)
    }
    R = r_not + 1/2 * sum_R
    R_vector[k] = R
  }
  return(R_vector)
}

update_beta <- function(beta_not, sigma_a_vector, mu_a_vector){
  beta <- beta_not + 1/2 * sum(sigma_a_vector + mu_a_vector^2)
  return(beta)
}

update_d_vector <- function(Y, K, probability_matrix, d_not_vector) {
  if(is.null(d_not_vector)){
    d_not_vector = get_d_not_vector(K)
    warning("No prior information for d_not, 1/K has been used for each cluster.")
  }
  
  d_vector = numeric(K)
  for (k in 1:K) {
    sum_d = 0
    for (i in 1:nrow(Y)) {
      p = probability_matrix[i, k]
      sum_d = sum_d + p
    }
    d = d_not_vector[k] + sum_d
    d_vector[k] = d
  }
  return(d_vector)
}

update_probability_matrix <- function(Y, K, x, sigma_list, m_list, A_vector, R_vector, d_vector, nbasis, sigma_a_vector, mu_a_vector) {
  probability_matrix = matrix(0, nrow(Y), K)
  n = ncol(Y)
  B = get_B(x, nbasis)
  
  for (k in 1:K) {
    for (i in 1:nrow(Y)) {
      ev_phi_a = sum(diag(B %*% sigma_list[[k]] %*% t(B))) + n * sigma_a_vector[i] + 
        t((Y[i, ] - B %*% t(m_list[[k]]) - mu_a_vector[i])) %*% (Y[i, ] - B %*% t(m_list[[k]]) - mu_a_vector[i])
      ev_tau_k_log_tau_k = digamma(A_vector[k]) - log(R_vector[k]) 
      ev_tau_k_tau_k = A_vector[k] / R_vector[k]
      ev_pi_log_pi_k = digamma(d_vector[k]) - digamma(sum((d_vector)))
      
      a_i_k = n / 2 * ev_tau_k_log_tau_k - 1/2 * ev_tau_k_tau_k * ev_phi_a + ev_pi_log_pi_k
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
