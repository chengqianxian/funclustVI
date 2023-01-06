library(fda)
library(sabre)

#' Runs and evaluates the model on the specified simulated function data 
#'
#' @param data_params List object containing the parameters required for generating the functional data and its characterics. Must include function named generate_data that accepts the data params list and returns a matrix of curves in the rows. Must also include vector named true_cluster_assignments that contains the actual cluster assignments for each curve. Must also include vector named seeds that contains seeds the seed for each iteration of the simulation. Other then these requirements, anything can be included in the data_params object. 
#' @param model_params List object containing the parameters required for modelling the data and generating the cluster assignments. Must include function named model_func that accepts Y, data_params and model_params and returns a vector containing the cluster assignment of each curve. Must also include list plot_params as with the requirements as referenced above. Other then these requirements, anything can be included in the model_params object. 
#' @param eval_func_list List object containing the functions corresponding to the various evaluations metrics evaluation the performance of the algorithim. Each function must accept a vector cluster_assignments that is generated from the function model_func in the model_params list as well as the data_params object which contains the vector true cluster assignments as referenced above. 
#' @param number_of_simulations The number of simulations 
#' @param save_path The file path to save the results from the simulations 



simulate <- function(data_params, model_params, eval_func_list, number_of_simulations, save_path) {
  start_time = Sys.time()
  num_of_eval_funcs = length(eval_func_list)
  eval_func_name_list = names(eval_func_list)
  eval_metric_sum_vector = c(1:num_of_eval_funcs)*0 
  final_res_mat = matrix(0, number_of_simulations, num_of_eval_funcs)
  
  count = 0
  for (simulation_number in 1:number_of_simulations) {
    seed = data_params$seeds[simulation_number]
    res = compute_function(seed, data_params, model_params, eval_func_list)
    
    cat("seed ", seed, ": ")
    for (i in 1:num_of_eval_funcs) {
      prev_sum = eval_metric_sum_vector[i]
      curr_val = res[[i]]
      final_res_mat[simulation_number, i] = curr_val
      eval_metric_sum_vector[i] = prev_sum + curr_val
      eval_func_name = eval_func_name_list[i]
      cat(eval_func_name, " = ", curr_val, " ")
    }
    cat("\n")
    
    count = count + 1
  }
  
  eval_metric_avg_vector = eval_metric_sum_vector / count
  
  for(i in 1:num_of_eval_funcs) {
    eval_func_name = eval_func_name_list[i]
    avg_val = eval_metric_avg_vector[i]
    cat("Average ", eval_func_name, " = ", avg_val, "\n")
  }
  
  if(is.null(save_path) == FALSE) {
    write(final_res_mat, save_path)
  }
  
  end_time = Sys.time()
  
  simulation_length = end_time - start_time 
  
  res = compute_function(data_params$seeds[1], data_params, model_params, eval_func_list) # return one set of simulation results
  result_list = c(list("result_matrix" = final_res_mat, "simulation_length" = simulation_length, "eval_metric_avg_vector" = eval_metric_avg_vector), res)
  return(result_list)
}

#' Generates simulated with a certain seed, evaluates the data using the model and returns a list with the results from the various evaluation metrics. 
#'
#' @param seed The seed 
#' @param data_params List object containing the parameters required for generating the functional data and its characterics. Must include function named generate_data that accepts the data params list and returns a matrix of curves in the rows. Must also include vector named true_cluster_assignments that contains the actual cluster assignments for each curve. Must also include vector named seeds that contains seeds the seed for each iteration of the simulation. Other then these requirements, anything can be included in the data_params object. 
#' @param model_params List object containing the parameters required for modelling the data and generating the cluster assignments. Must include function named model_func that accepts Y, data_params and model_params and returns a vector containing the cluster assignment of each curve. Must also include list plot_params as with the requirements as referenced above. Other then these requirements, anything can be included in the model_params object. 
#' @param eval_func_list List object containing the functions corresponding to the various evaluations metrics evaluation the performance of the algorithim. Each function must accept a vector cluster_assignments that is generated from the function model_func in the model_params list as well as the data_params object which contains the vector true cluster assignments as referenced above. 

compute_function <- function(seed, data_params, model_params, eval_func_list) {
  set.seed(seed)
  num_of_eval_funcs = length(eval_func_list)
  eval_func_name_list = names(eval_func_list)
  
  generate_data = data_params$generate_data
  
  Y = generate_data(data_params)
  
  eval_metric_res_vector = c(1:num_of_eval_funcs)*0
  
  model_func = model_params$model_func
  
  clf = model_func(Y, data_params, model_params)
  cluster_assignments = clf$cluster_assignments
  
  result_list = list()
  
  for(i in 1:num_of_eval_funcs) {
    eval_func_name = eval_func_name_list[i]
    eval_func = eval_func_list[[i]]
    res = eval_func(cluster_assignments, data_params)
    result_list[[eval_func_name]] = res
  }
  
  return(c(result_list, clf))
}

#' Model function wrapper for the funclustVI 
#'
#' @param Y A matrix in which the rows represent the curves 
#' @param data_params List object containing the parameters required for generating the functional data and its characterics
#' @param model_params List object containing the parameters required for model the data and generating the cluster assignments


get_funclustVI_cluster_assignments <- function(Y, data_params, model_params) {
  x = data_params$x
  K = data_params$K
  init = model_params$init
  nbasis = model_params$nbasis
  convergence_threshold = model_params$convergence_threshold
  max_iterations = model_params$max_iterations
  gamma_dist_config_matrix = model_params$gamma_dist_config_matrix
  d_not_vector = model_params$d_not_vector
  m_not_vector = model_params$m_not_vector
  v_not_vector = model_params$v_not_vector
  plot_params = model_params$plot_params
  true_cluster_assignments = data_params$true_cluster_assignments
  verbose = model_params$verbose 
  draw = model_params$draw
  clf = funcslustVI(x, Y, K, true_cluster_assignments, init, nbasis, convergence_threshold, max_iterations, gamma_dist_config_matrix, 
                    d_not_vector, m_not_vector, v_not_vector,
                    verbose, draw, plot_params)
  return(clf)
}






