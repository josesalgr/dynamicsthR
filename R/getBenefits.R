#' Solve a mathematical optimization problem using specified solver.
#'
#' This function takes an optimization problem represented by a list of matrices and vectors and solves it using either the CBC (Coin-OR branch-and-cut) or Gurobi solver.
#'
#' @param a A list containing information about the optimization problem, including coefficients, constraints, and other parameters.
#' @param solver A character specifying the solver to be used ("gurobi" or "cbc").
#' @param cores An integer specifying the number of processor cores to be used.
#' @param verbose A logical value indicating whether to print solver output.
#' @param gap_limit A numeric scalar specifying the relative MIP (Mixed-Integer Programming) optimality gap as a stopping condition.
#' @param time_limit An integer specifying the time limit for the solver in seconds.
#'
#' @return A list containing the solution to the optimization problem, including variable values and solver statistics.
#' @export
#' @examples
#' solve(a, solver = "gurobi", cores = 2, verbose = TRUE, gap_limit = 0.0, time_limit = 100)
getBenefits <- function(sol) {
  
  if(sol$solver == "gurobi"){
    name_sol = "x"
  }
  else if(sol$solver == "cbc"){
    name_sol = "column_solution"
  }
  else{
    stop("Invalid solution input")
  }
  
  # Variable extraction
  #Variable X
  # size_x <- length(stringr::str_subset(sol$model$n_variable, "x"))
  # indexes_x <- which(!is.na(stringr::str_extract(sol$model$n_variable, "x")))
  # var_x <- data.frame(id = 1:size_x,
  #                     i = rep(0, size_x),
  #                     k = rep(0, size_x),
  #                     t = rep(0, size_x),
  #                     n = rep(0, size_x),
  #                     sol = round(sol$solution[name_sol][[1]][indexes_x[1]:indexes_x[length(indexes_x)]]), 0)
  # 
  # for(i in indexes_x){
  #   id <- which(indexes_x == i)
  #   indexes <- stringr::str_locate_all(sol$model$n_variable[i], "_")
  #   index_i <- indexes[[1]][1]
  #   index_k <- indexes[[1]][2]
  #   index_t <- indexes[[1]][3]
  #   index_n <- indexes[[1]][4]
  #   
  #   var_x$i[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_i+1, end = index_k-1)) + 1
  #   var_x$k[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_k+1, end = index_t-1)) + 1
  #   var_x$t[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_t+1, end = index_n-1)) + 1
  #   var_x$n[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_n+1, end = stringr::str_length(sol$model$n_variable[i])))
  # }
  # 
  
  #Variable V
  # size_v <- length(stringr::str_subset(sol$model$n_variable, "v"))
  # indexes_v <- which(!is.na(stringr::str_extract(sol$model$n_variable, "v")))
  # var_v <- data.frame(id = indexes_v,
  #                     i = rep(0, size_v),
  #                     k = rep(0, size_v),
  #                     t = rep(0, size_v),
  #                     n = rep(0, size_v),
  #                     sol = round(sol$solution[name_sol][[1]][indexes_v[1]:indexes_v[length(indexes_v)]]), 0)
  # 
  # for(i in indexes_v){
  #   id <- which(indexes_v == i)
  #   indexes <- stringr::str_locate_all(sol$model$n_variable[i], "_")
  #   index_i <- indexes[[1]][1]
  #   index_k <- indexes[[1]][2]
  #   index_t <- indexes[[1]][3]
  #   index_n <- indexes[[1]][4]
  #   
  #   var_v$i[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_i+1, end = index_k-1)) + 1
  #   var_v$k[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_k+1, end = index_t-1)) + 1
  #   var_v$t[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_t+1, end = index_n-1)) + 1
  #   var_v$n[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_n+1, end = stringr::str_length(sol$model$n_variable[i])))
  # }
  # 
  
  #Variable w
  # size_w <- length(stringr::str_subset(model$n_variable, "w"))
  # indexes_w <- which(!is.na(stringr::str_extract(model$n_variable, "w")))
  # var_w <- data.frame(id = indexes_w,
  #                     i = rep(0, size_w),
  #                     t = rep(0, size_w),
  #                     sol = round(solution$x[indexes_w[1]:indexes_w[length(indexes_w)]]), 0)
  # 
  # for(i in indexes_w){
  #   id <- which(indexes_w == i)
  #   indexes <- stringr::str_locate_all(model$n_variable[i], "_")
  #   index_i <- indexes[[1]][1]
  #   index_t <- indexes[[1]][2]
  #   
  #   var_w$i[id] <- as.numeric(stringr::str_sub(model$n_variable[i], start = index_i+1, end = index_t-1)) + 1
  #   var_w$t[id] <- as.numeric(stringr::str_sub(model$n_variable[i], start = index_t+1, end = stringr::str_length(model$n_variable[i]))) + 1
  # }
  
  #Variable y
  size_y <- length(stringr::str_subset(sol$model$n_variable, "y"))
  indexes_y <- which(!is.na(stringr::str_extract(sol$model$n_variable, "y")))
  var_y <- data.frame(id = indexes_y,
                      i = rep(0, size_y),
                      t = rep(0, size_y),
                      sol = sol$solution[name_sol][[1]][indexes_y[1]:indexes_y[length(indexes_y)]])

  for(i in indexes_y){
    id <- which(indexes_y == i)
    indexes <- stringr::str_locate_all(sol$model$n_variable[i], "_")
    index_i <- indexes[[1]][1]
    index_s <- indexes[[1]][2]
    index_t <- indexes[[1]][3]

    var_y$i[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_i+1, end = index_s-1))
    var_y$s[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_s+1, end = index_t-1))
    var_y$t[id] <- as.numeric(stringr::str_sub(sol$model$n_variable[i], start = index_t+1, end = stringr::str_length(sol$model$n_variable[i])))
  }

  
  # #Variable bd
  # size_bd <- length(stringr::str_subset(model$n_variable, "bd"))
  # indexes_bd <- which(!is.na(stringr::str_extract(model$n_variable, "bd")))
  # var_bd <- data.frame(id = indexes_bd,
  #                      t = rep(0, size_bd),
  #                      sol = solution$x[indexes_bd[1]:indexes_bd[length(indexes_bd)]])
  # 
  # for(i in indexes_bd){
  #   id <- which(indexes_bd == i)
  #   indexes <- stringr::str_locate_all(model$n_variable[i], "_")
  #   index_t <- indexes[[1]][1]
  #   
  #   var_bd$t[id] <- as.numeric(stringr::str_sub(model$n_variable[i], start = index_t+1, end = stringr::str_length(model$n_variable[i]))) + 1
  # }
  
  return(var_y)
}