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
model <- function(input = NULL, levels = 3L, periods = 10L, budget = 0) {
  
  # assert that arguments are valid
  assertthat::assert_that(
    !is.null(input),
    assertthat::is.scalar(budget),
    assertthat::is.scalar(periods),
    is.integer(periods),
    assertthat::is.scalar(levels),
    is.integer(levels),
    is.list(input)
  )
  
  if(inherits(input, "list")){

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
                                        levels, periods, budget)
    
    return(model)
  }

}