get_elbo <- function(x, Y, K, d_not_vector, v_not_vector, b_not, r_not, m_not_vector, alpha_not, beta_not, 
                     nbasis, sigma_list, m_list, A_vector, R_vector, d_vector, sigma_a_vector, mu_a_vector, 
                     alpha, beta, probability_matrix) {
  B = get_B(x, nbasis)
  n <- ncol(Y)
  diff_pi <- sum((d_not_vector - d_vector) * (digamma(d_vector) - digamma(sum(d_vector)))) 
  
  diff_z_sum_1 = 0
  diff_z_sum_2 = 0
  for (i in 1:nrow(Y)) {
    for (k in 1:K) {
      ev_pi_log_pi_k = digamma(d_vector[k]) - digamma(sum((d_vector)))
      diff_z_sum_1 = probability_matrix[i, k] * ev_pi_log_pi_k + diff_z_sum_1
      if (probability_matrix[i, k] <= .0000001) {
        diff_z_sum_2 = 0 + diff_z_sum_2
      } else {
        diff_z_sum_2 = probability_matrix[i, k] * log(probability_matrix[i, k]) + diff_z_sum_2
      }
    }
  }
  diff_z = diff_z_sum_1 - diff_z_sum_2
  
  diff_phi_sum <- 0
  for (k in 1:K) {
    diff_phi_sum <-  (- v_not_vector[k] * (sum(diag(sigma_list[[k]])) + (m_list[[k]] - m_not_vector[k, ]) %*% 
                                                t(m_list[[k]] - m_not_vector[k, ])) + log(det(sigma_list[[k]]))) / 2 + diff_phi_sum
  }
  diff_phi <- diff_phi_sum
  
  diff_tau_sum_1 = 0
  diff_tau_sum_2 = 0
  for (k in 1:K) {
    ev_tau_k_log_tau_k = digamma(A_vector[k]) - log(R_vector[k])
    diff_tau_sum_1 = (b_not - 1) * ev_tau_k_log_tau_k - r_not * (A_vector[k] / R_vector[k]) + diff_tau_sum_1
    # diff_tau_sum_2 = log(dgamma(ev_tau_k_log_tau_k, A_vector[k], R_vector[k])) + diff_tau_sum_2
    # diff_tau_sum_2 = A_vector[k] * (log(R_vector[k]) - 1) + (A_vector[k] - 1) * ev_tau_k_log_tau_k + diff_tau_sum_2
    diff_tau_sum_2 = 1 / 2 * log(A_vector[k]) + 1 / (2 * A_vector[k]) + diff_tau_sum_2
  }
  diff_tau = diff_tau_sum_1 - diff_tau_sum_2
  
  diff_a = - 1/ 2 * alpha / beta * sum(sigma_a_vector +  mu_a_vector^2) + sum(log(sqrt(sigma_a_vector)))
  
  diff_tau_a <- (alpha_not - alpha) * (digamma(alpha) - log(beta)) - (beta_not - beta) * alpha / beta - alpha * log(beta)
  
  ev_log_lik_sum = 0
  for (i in 1:nrow(Y)) {
    for (k in 1:K) {
      ev_tau_k_log_tau_k = digamma(A_vector[k]) - log(R_vector[k])
      ev_phi_a = sum(diag(B %*% sigma_list[[k]] %*% t(B))) + n * sigma_a_vector[i] + 
        t((Y[i, ] - B %*% t(m_list[[k]]) - mu_a_vector[i])) %*% (Y[i, ] - B %*% t(m_list[[k]]) - mu_a_vector[i])
      ev_log_lik_sum = probability_matrix[i, k] * (ncol(Y) / 2 * ev_tau_k_log_tau_k - 1/2 * (A_vector[k] / R_vector[k]) * ev_phi_a) + ev_log_lik_sum
    }
  }
  ev_log_lik = as.numeric(ev_log_lik_sum)
  
  elbo = diff_pi + diff_z + diff_phi + diff_tau + diff_a +  diff_tau_a + ev_log_lik
  return(c(elbo, ev_log_lik))
}


check_convergence <- function(prev_elbo, curr_elbo, convergence_threshold) {
  if(is.null(prev_elbo) == TRUE) {
    return(FALSE)
  }
  else{
    dif = abs(curr_elbo - prev_elbo)
    if(dif  <= convergence_threshold) return(TRUE)
    else return(FALSE)
  }
}
