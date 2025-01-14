# To keep all functions for prepare running simulation ----------

## For generate Qinf ---------
numsim <- 1000 # no. of iteration
totpoints <- nrow(farm_2024)
set.seed(001)
Q_init_matrix <- matrix(0,nrow=numsim,ncol=totpoints)
for (ii in 1:numsim){
  Q_init_matrix[ii,] <- rexp(totpoints, rate = 1) # random threshold to infection
}