library(gtools)

Mismatch <- function(predicited_clusters, true_clusters, K)
{
  sigma = gtools::permutations(n=K,r=K,v=1:K)
  
  Miss = length(which( true_clusters != predicited_clusters))  ## for permutation 1, 2,... K
  
  mm_aux = predicited_clusters
  for (ind in 2:dim(sigma)[1])
  {
    
    for (j in 1:K)
      mm_aux[which(predicited_clusters == j)] = sigma[ind,j]
    
    Miss[ind] =  length(which( true_clusters != mm_aux))
    
  }
  mis.percent <- min(Miss) / length(predicited_clusters) # mistmatch rate
  return(round(mis.percent, 4))
}


#' Gets the deviance information criterion (DIC) for a clustering result
#'
#' @param gamma_dist_config_matrix The priors of Gamma distribution, first row for alpha, second row for beta
#' @param A_vector The vector A star
#' @param cluster_assignments The resulted cluster assignments from the model 
#' @param Y The observed response data for clustering 
#' @param B The matrix of B-splines smoothing
#' @param m_list The list of converged coefficients of basis functions for each cluster
#' @param D One of component of converged ELBO used to calculating DIC
#' @return dic

get_dic <- function(gamma_dist_config_matrix, A_vector, cluster_assignments, Y, B, m_list, D){
  r_0 <- gamma_dist_config_matrix[2, ]
  r_star <- numeric(length(cluster_assignments))
  for (i in 1:length(cluster_assignments)) {
    j <- cluster_assignments[i]
    r_star[i] <- t((Y[i, ] - B %*% t(m_list[[j]]))) %*% (Y[i, ] - B %*% t(m_list[[j]]))
  }
  R_star <- r_0 + 0.5 * sum(r_star)
  A_over_B <- A_vector / R_star
  
  D_i_bar <- numeric(length(cluster_assignments))
  for (i in 1:length(cluster_assignments)) {
    j <- cluster_assignments[i]
    D_i_bar[i] <- 0.5 * log(A_over_B[j]) - 0.5 * A_over_B[j] *(t((Y[i, ] - B %*% t(m_list[[j]]))) %*% (Y[i, ] - B %*% t(m_list[[j]])))
  }
  D_bar <- sum(D_i_bar)
  
  dic <- -4 * D + 2 * D_bar
  return(dic)
}

#' Gets the number of mismatches for a single iteration of the simulation 
#'
#' @param cluster_assignments A vector where each entry is the cluster assignment for the corresponging curve 
#' @param K The number of clusters in the data
#' @param curves_per_cluster The number of curves per cluster 
#'
#' @return The number of mismatches 
 
get_mismatches <- function(cluster_assignments, data_params) {
  true_cluster_assignments = data_params$true_cluster_assignments
  mismatches = Mismatch(cluster_assignments, true_cluster_assignments, K)
  return(mismatches)
}

#' Gets the number of vmeaure for a single iteration of the simulation 
#'
#' @param cluster_assignments A vector where each entry is the cluster assignment for the corresponging curve 
#' @param curves_per_cluster The number of curves per cluster 
#'
#' @return The vmeasure
get_v_measure <- function(cluster_assignments, data_params) {
  true_cluster_assignments = data_params$true_cluster_assignments
  v_measure = sabre::vmeasure(cluster_assignments, true_cluster_assignments)$v_measure
  return(v_measure)
}



