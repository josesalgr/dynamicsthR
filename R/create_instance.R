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
create_instance <- function(type = NULL, file = "", 
                            n = 10L, 
                            species = 4, 
                            threats = 2, 
                            expansionType = NULL,
                            threats_lloc = NULL, 
                            dlong = NULL, 
                            seed = 49) {
  
  # assert that arguments are valid
  assertthat::assert_that(
    assertthat::is.string(type),
    assertthat::is.scalar(n)
  )
  
  instance <- list()
  
  if(type == "square"){

    set.seed(seed)
    phi <- 0.1
    instance$Unidades <- n
    instance$Especies <- species
    instance$Amenazas <- threats
    
    # ExpansionType
    expType <- matrix(0, nrow = threats, ncol = 3)
    expType[, 1] <- 1
    instance$ExpansionType <- tidyr::tibble(expType)
    
    # Define function to draw random samples from a multivariate normal
    # distribution
    rmvn <- function(n, mu = 0, V = matrix(1)) {
      p <- length(mu)
      if (any(is.na(match(dim(V), p)))) 
        stop("Dimension problem!")
      D <- chol(V)
      t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
    }
    
    # Set up a square lattice region
    simgrid <- expand.grid(1:n, 1:n)
    n_total <- nrow(simgrid)
    
    # Set up distance matrix
    distances_matrix <- as.matrix(dist(simgrid))
    #distances_matrix[lower.tri(distances_matrix)] <- NA
    
    # distances_matrix_df <- as.data.frame(distances_matrix)
    # 
    # # Add row names as a separate column
    # distances_matrix_df <- distances_matrix_df %>%
    #   dplyr::mutate(row = rownames(.)) %>%
    #   dplyr::select(row, everything())
    # 
    # distances_matrix_df <- tidyr::pivot_longer(distances_matrix_df, 
    #                                     cols = -row, 
    #                                     names_to = "col", 
    #                                     names_prefix = "V", 
    #                                     values_drop_na = TRUE)
    
    instance$Dradial <- tidyr::tibble(distances_matrix)
    
    # Cost simulation
    X_sim <- rmvn(1, rep(0, n_total), exp(-phi * distances_matrix))
    X_sim <- abs(X_sim)
    Ck <- round(X_sim*10/max(X_sim), 2)
    instance$Ck <- tidyr::tibble(Ck)
    
    # Species distribution model
    df_species <- data.frame(id = 1:n_total)

    for(i in 1:species){
      # Generate random variable
      X_sim <- rmvn(1, rep(0, n_total), exp(-phi * distances_matrix))
      
      X_sim_dic <- ifelse(X_sim < 0.0, 0, 1)
      # Visualize results
      df_species <- cbind(df_species, X_sim_dic)
    }
    instance$Ij <- tidyr::tibble(df_species, .name_repair = "minimal")
    
    # Threat distribution model
    df_threats <- data.frame(id = 1:n_total)
    
    if(is.null(threats_lloc)){
      for(k in 1:threats){
        X_sim <- sample(c(rep(0, n_total-1), 1))
        df_threats <- cbind(df_threats, X_sim)
      }
    }
    else{
      for(k in 1:threats){
        df_threats <- cbind(df_threats, c(threats_lloc))
      }
    }
    instance$Ik <- tidyr::tibble(df_threats, .name_repair = "minimal")
    
    # Sensitivity
    instance$Jk <- tidyr::tibble(matrix(1, nrow = species, ncol = threats))
    
    # Adyacency matrix
    # Function to check if two points are adjacent
    are_adjacent <- function(p1, p2) {
      abs(p1[1] - p2[1]) <= 1 & abs(p1[2] - p2[2]) <= 1
    }
    
    # Initialize an adjacency matrix
    adjacency_matrix <- matrix(NA, nrow = nrow(simgrid), ncol = nrow(simgrid))
    
    # Populate the adjacency matrix
    for (i in 1:nrow(simgrid)) {
      for (j in 1:nrow(simgrid)) {
        if (are_adjacent(simgrid[i, ], simgrid[j, ]) && i!=j) {
          adjacency_matrix[i, j] <- 1
        }
      }
    }
    adjacency_df <- as.data.frame(adjacency_matrix)
    
    # Add row names as a separate column
    adjacency_df <- adjacency_df %>%
      dplyr::mutate(row = rownames(.)) %>%
      dplyr::select(row, everything())
    
    adjacency_df <- tidyr::pivot_longer(adjacency_df, 
                        cols = -row, 
                        names_to = "col", 
                        names_prefix = "V", 
                        values_drop_na = TRUE)
    
    adjacency_df[,1] <- as.numeric(adjacency_df[,1][[1]])
    adjacency_df[,2] <- as.numeric(adjacency_df[,2][[1]])
    instance$Adyacency <- tidyr::tibble(adjacency_df[,-3])
    
    instance$Dlong <- dlong
  }
  else if(type == "specific"){
    
    df <- read.table(file, sep = "\t")
    
    ind_PU <- which(df == "Unidades")
    ind_Species <- which(df == "Especies")
    ind_Threats <- which(df == "Amenazas")
    ind_ExpansionType <- which(df == "ExpansionType")
    ind_Dradial <- which(df == "Dradial")
    ind_Ck <- which(df == "Ck")
    ind_Ij <- which(df == "Ij")
    ind_Ik <- which(df == "Ik")
    ind_Jk <- which(df == "Jk")
    ind_Adyacency <- which(df == "Adyacency")
    ind_Dlong <- which(df == "Dlong")
    
    instance$Unidades <- as.numeric(df[ind_PU + 1, 1])
    instance$Especies <- as.numeric(df[ind_Species + 1, 1])
    instance$Amenazas <- as.numeric(df[ind_Threats + 1, 1])

    # ExpansionType
    expType <- df[(ind_ExpansionType + 1):(ind_Dradial - 1), 1]
    expType_df <- as.data.frame(expType)
    expType_df <- tidyr::separate(expType_df, col = "expType", sep = " ", into = paste0("V", 1:3), convert = TRUE)
    instance$ExpansionType <- tidyr::tibble(expType_df)
    
    # Radial distance
    distances <- df[(ind_Dradial + 1):(ind_Ck - 1), 1]
    distances_df <- as.data.frame(distances)
    distances_df <- tidyr::separate(distances_df, col = "distances", sep = " ", into = paste0("V", 1:instance$Unidades), convert = TRUE)
    instance$Dradial <- tidyr::tibble(distances_df)

    # Cost
    Ck <- df[(ind_Ck + 1):(ind_Ij - 1), 1]
    instance$Ck <- tidyr::tibble(as.numeric(Ck))
    
    # Species distribution model
    Ij <- df[(ind_Ij + 1):(ind_Ik - 1), 1]
    Ij_df <- as.data.frame(Ij)
    Ij_df <- tidyr::separate(Ij_df, col = "Ij", sep = " ", into = paste0("V", 1:instance$Especies), convert = TRUE)
    instance$Ij <- tidyr::tibble(Ij_df)
    
    # Threat distribution model
    Ik <- df[(ind_Ik + 1):(ind_Jk - 1), 1]
    Ik_df <- as.data.frame(Ik)
    Ik_df <- tidyr::separate(Ik_df, col = "Ik", sep = " ", into = paste0("V", 1:instance$Amenazas), convert = TRUE)
    instance$Ik <- tidyr::tibble(Ik_df)
    
    # Sensitivity
    Jk <- df[(ind_Jk + 1):(ind_Adyacency - 1), 1]
    Jk_df <- as.data.frame(Jk)
    Jk_df <- tidyr::separate(Jk_df, col = "Jk", sep = " ", into = paste0("V", 1:instance$Amenazas), convert = TRUE)
    instance$Jk <- tidyr::tibble(Jk_df)
    
    # Adyacency matrix
    distances <- df[(ind_Adyacency + 1):(ind_Dlong - 1), 1]
    distances_df <- as.data.frame(distances)
    distances_df <- tidyr::separate(distances_df, col = "distances", sep = " ", into = paste0("V", 1:2), convert = TRUE)
    instance$Adyacency <- tidyr::tibble(distances_df)
    
    # Upstream/downstream distance
    distances <- df[(ind_Dlong + 1):(nrow(df)), 1]
    distances_df <- as.data.frame(distances)
    distances_df <- tidyr::separate(distances_df, col = "distances", sep = " ", into = paste0("V", 1:3), convert = TRUE)
    instance$Dlong <- tidyr::tibble(distances_df)
    
    
  }
  else{
    stop("Invalid type.")
  }
  
  return(instance)
}