# To keep all functions for prepare running simulation ----------

## For generate Qinf ---------
Q_init_func <- function(farm_dat, numsim ){
  #' Create Q_init Matrix
  #'
  #' This function generates a `Q_init` matrix for the Sellke construction modeling. 
  #' Each row of the matrix corresponds to a set of random thresholds for infection, 
  #' drawn from an exponential distribution with a rate parameter of 1.
  #'
  #' @param farm_dat A data frame or matrix containing the farm data. The number of rows in this data determines the number of columns in the `Q_init` matrix.
  #' @param numsim An integer specifying the number of simulations (rows) to include in the `Q_init` matrix.
  #'
  #' @return A matrix with `numsim` rows and `nrow(farm_dat)` columns, where each entry is a randomly generated threshold for infection.
  #'
  #' @examples
  #' # Example usage:
  #' farm_data <- data.frame(id = 1:10) # Example farm data with 10 farms
  #' num_simulations <- 5
  #' Q_matrix <- Q_init_func(farm_data, num_simulations)
  #' print(Q_matrix)
  #' 
  Q_init_matrix <- matrix(0,nrow=numsim,ncol=nrow(farm_dat))
  for (ii in 1:numsim){
    Q_init_matrix[ii,] <- rexp(nrow(farm_dat), rate = 1) # random threshold to infection
  }
  return(Q_init_matrix)
}

docstring(Q_init_func) 

?Q_init_func


Q <- Q_init_func(farm_dat = farm_2024, numsim = 1000)

