#' Do values in a numeric vector fall in specified range?
#'
#' This is a shortcut for `x >= left & x <= right`, implemented
#' efficiently in C++ for local values, and translated to the
#' appropriate SQL for remote tables.
#'
#' @param x A numeric vector of values.
#' @param left,right Boundary values.
#' @export
#' @examples
#' between(1:12, 7, 9)
#'
#' x <- rnorm(1e2)
#' x[between(x, -1, 1)]
solve <- function(a, solver = "gurobi", cores = 2, verbose = TRUE, gap_limit = 0.0, time_limit = 100) {
  
  model = a
  model$A <- Matrix::sparseMatrix(i = model$A_i + 1, j = model$A_j + 1, x = model$A_x, dims = c(length(a$rhs + 1), length(a$obj + 1)))
  
  if(solver == "cbc")
  {
    ## cbc model
    constraints_minus_equal <- which(model$sense == "<=")
    constraints_plus_equal <- which(model$sense == ">=")
    row_ub <- model$rhs
    row_ub[constraints_plus_equal] <- Inf
    row_lb <- model$rhs
    row_lb[constraints_minus_equal] <- -Inf
    
    ## cbc parameters
    cbc_args <- list()
    cbc_args$threads <- cores
    cbc_args$log <- as.integer(verbose)
    cbc_args$verbose <- 15
    # Stop condition: Relative MIP optimality gap
    cbc_args$ratio <- gap_limit
    # Stop condition: Time limit
    cbc_args$sec <- time_limit
    cbc_args$timem <- "elapsed"
    # Activate heuristics methods
    cbc_args$heuristicsOnOff <- "on"
    
    runtime_cbc <- system.time(
      solution <- rcbc::cbc_solve(obj = model$obj,
                                  mat = model$A,
                                  is_integer = ifelse(model$vtype == "B", TRUE, FALSE),
                                  row_ub = row_ub,
                                  row_lb = row_lb,
                                  col_lb = model$lb,
                                  col_ub = model$ub,
                                  max = ifelse(model$modelsense == "min", FALSE, TRUE),
                                  cbc_args = cbc_args)
    )[[1]]
    
    #Gap_limit
    if(isTRUE(solution$status == 0L)){
      solution$gap <- gap_limit
    }
    else{
      solution$gap <- "No reported"
    }
  }
  else if(solver == "gurobi")
  {
    
    ## Gurobi model
    model$sense <- replace(model$sense, model$sense == "==", "=")
    
    ## Gurobi parameters
    params <- list()
    params$Threads <- cores
    params$LogToConsole <- as.integer(verbose)
    params$NodefileStart <- 0.5
    # Stop condition: Relative MIP optimality gap
    params$MIPGap <- gap_limit
    # Stop condition: Time limit
    params$TimeLimit <- time_limit
    
    solution <- gurobi::gurobi(model, params)
  }
  
  #extraction
  #Variable X
  size_x <- length(str_subset(model$n_variable, "x"))
  indexes_x <- which(!is.na(str_extract(model$n_variable, "x")))
  var_x <- data.frame(id = 1:size_x,
                      i = rep(0, size_x),
                      k = rep(0, size_x),
                      t = rep(0, size_x),
                      n = rep(0, size_x),
                      sol = round(solution$column_solution[indexes_x[1]:indexes_x[length(indexes_x)]]), 0)
  
  for(i in indexes_x){
    id <- which(indexes_x == i)
    indexes <- str_locate_all(model$n_variable[i], "_")
    index_i <- indexes[[1]][1]
    index_k <- indexes[[1]][2]
    index_t <- indexes[[1]][3]
    index_n <- indexes[[1]][4]
    
    var_x$i[id] <- as.numeric(str_sub(model$n_variable[i], start = index_i+1, end = index_k-1)) + 1
    var_x$k[id] <- as.numeric(str_sub(model$n_variable[i], start = index_k+1, end = index_t-1)) + 1
    var_x$t[id] <- as.numeric(str_sub(model$n_variable[i], start = index_t+1, end = index_n-1)) + 1
    var_x$n[id] <- as.numeric(str_sub(model$n_variable[i], start = index_n+1, end = str_length(model$n_variable[i])))
  }
  
  
  #Variable V
  size_v <- length(str_subset(model$n_variable, "v"))
  indexes_v <- which(!is.na(str_extract(model$n_variable, "v")))
  var_v <- data.frame(id = indexes_v,
                      i = rep(0, size_v),
                      k = rep(0, size_v),
                      t = rep(0, size_v),
                      n = rep(0, size_v),
                      sol = round(solution$column_solution[indexes_v[1]:indexes_v[length(indexes_v)]]), 0)
  
  for(i in indexes_v){
    id <- which(indexes_v == i)
    indexes <- str_locate_all(model$n_variable[i], "_")
    index_i <- indexes[[1]][1]
    index_k <- indexes[[1]][2]
    index_t <- indexes[[1]][3]
    index_n <- indexes[[1]][4]
    
    var_v$i[id] <- as.numeric(str_sub(model$n_variable[i], start = index_i+1, end = index_k-1)) + 1
    var_v$k[id] <- as.numeric(str_sub(model$n_variable[i], start = index_k+1, end = index_t-1)) + 1
    var_v$t[id] <- as.numeric(str_sub(model$n_variable[i], start = index_t+1, end = index_n-1)) + 1
    var_v$n[id] <- as.numeric(str_sub(model$n_variable[i], start = index_n+1, end = str_length(model$n_variable[i])))
  }
  
  
  #Variable w
  size_w <- length(str_subset(model$n_variable, "w"))
  indexes_w <- which(!is.na(str_extract(model$n_variable, "w")))
  var_w <- data.frame(id = indexes_w,
                      i = rep(0, size_w),
                      t = rep(0, size_w),
                      sol = round(solution$column_solution[indexes_w[1]:indexes_w[length(indexes_w)]]), 0)
  
  for(i in indexes_w){
    id <- which(indexes_w == i)
    indexes <- str_locate_all(model$n_variable[i], "_")
    index_i <- indexes[[1]][1]
    index_t <- indexes[[1]][2]
    
    var_w$i[id] <- as.numeric(str_sub(model$n_variable[i], start = index_i+1, end = index_t-1)) + 1
    var_w$t[id] <- as.numeric(str_sub(model$n_variable[i], start = index_t+1, end = str_length(model$n_variable[i]))) + 1
  }
  
  #Variable y
  # size_y <- length(str_subset(model$n_variable, "y"))
  # indexes_y <- which(!is.na(str_extract(model$n_variable, "y")))
  # var_y <- data.frame(id = indexes_y,
  #                     i = rep(0, size_y),
  #                     t = rep(0, size_y),
  #                     sol = solution$column_solution[indexes_y[1]:indexes_y[length(indexes_y)]])
  # 
  # for(i in indexes_y){
  #   id <- which(indexes_y == i)
  #   indexes <- str_locate_all(model$n_variable[i], "_")
  #   index_i <- indexes[[1]][1]
  #   index_s <- indexes[[1]][2]
  #   index_t <- indexes[[1]][3]
  #   
  #   var_y$i[id] <- as.numeric(str_sub(model$n_variable[i], start = index_i+1, end = index_s-1))
  #   var_y$s[id] <- as.numeric(str_sub(model$n_variable[i], start = index_s+1, end = index_t-1))
  #   var_y$t[id] <- as.numeric(str_sub(model$n_variable[i], start = index_t+1, end = str_length(model$n_variable[i])))
  # }
  # 
  
  #Variable bd
  size_bd <- length(str_subset(model$n_variable, "bd"))
  indexes_bd <- which(!is.na(str_extract(model$n_variable, "bd")))
  var_bd <- data.frame(id = indexes_bd,
                       t = rep(0, size_bd),
                       sol = solution$column_solution[indexes_bd[1]:indexes_bd[length(indexes_bd)]])
  
  for(i in indexes_bd){
    id <- which(indexes_bd == i)
    indexes <- str_locate_all(model$n_variable[i], "_")
    index_t <- indexes[[1]][1]
    
    var_bd$t[id] <- as.numeric(str_sub(model$n_variable[i], start = index_t+1, end = str_length(model$n_variable[i]))) + 1
  }
  var_sol <- list(var_x, var_v, var_w, var_bd)
  
  return(var_sol)
  
}