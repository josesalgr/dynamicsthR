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
solve <- function(a, solver = "gurobi", cores = 2L, verbose = TRUE, gap_limit = 0.0, time_limit = 100) {
  
  # assert that arguments are valid
  assertthat::assert_that(
    is.list(a),
    solver %in% c("gurobi", "cbc"),
    is.integer(as.integer(cores)),
    is.logical(verbose),
    is.numeric(gap_limit),
    is.numeric(time_limit)
  )
  
  # Rename the input list for better readability
  model = a
  
  # Convert sparse matrix representation
  model$A <- Matrix::sparseMatrix(i = model$A_i + 1, j = model$A_j + 1, x = model$A_x, dims = c(length(a$rhs + 1), length(a$obj + 1)))
  
  # Choose solver-specific implementation
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
  
  sol <- list(model = model, 
              solution = solution,
              solver = solver,
              hola = 1)
  
  return(sol)
}