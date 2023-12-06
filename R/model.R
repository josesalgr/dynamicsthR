#' Create a dynamic model based on input parameters.
#'
#' This function takes input parameters and uses them to create a dynamic model. The model is designed for dynamic actions, and the input should be a list containing various components.
#'
#' @param input A list containing the necessary components for creating the dynamic model. Components include Unidades, Especies, Amenazas, Jk, Ij, Ik, ExpansionType, Dlong, Adyacency, Dradial, and Ck.
#' @param objective A string specifying the type of objective to achieve in the mathematical model. There are two available: "min_prop" and "max_relative_benefit".
#' @param levels An integer specifying the number of levels in the model.
#' @param periods An integer specifying the number of periods for the dynamic model.
#' @param budget A numeric scalar specifying the budget for the model.
#'
#' @return A dynamic model based on the input parameters.
#' @export
#' @examples
#' model(input_data, levels = 3, periods = 10, budget = 100)
model <- function(input = NULL, objective = "min_prop", levels = 3L, periods = 10L, budget = 0) {
  
  # assert that arguments are valid
  assertthat::assert_that(
    !is.null(input),
    is.numeric(budget),
    assertthat::is.scalar(periods),
    is.integer(as.integer(periods)),
    assertthat::is.scalar(levels),
    is.integer(as.integer(levels)),
    is.list(input)
  )
  
  # Check if input is a list
  if(inherits(input, "list")){
    
    if(objective %in% c("min_prop", "max_relative_benefit")){
      
      if(objective == "min_prop"){
        objective = 1
      }
      else{
        objective = 2
      }
      # Create the dynamic model using the provided function (createModel_DynamicActions)
      model <- createModel_DynamicActions(input$Unidades,
                                          input$Especies,
                                          input$Amenazas,
                                          input$Jk,
                                          input$Ij,
                                          input$Ik,
                                          input$ExpansionType,
                                          input$Dlong,
                                          input$Adyacency,
                                          input$Dradial,
                                          input$Ck,
                                          levels, periods, budget, objective)
    }
    else{
      stop("Error: Objective not found.")
    }
    return(model)
  }
  
}