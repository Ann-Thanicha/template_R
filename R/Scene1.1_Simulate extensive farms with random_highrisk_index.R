### Simulate Scen1.1 randomly outdoor high index ####

## First specify the packages of interest
packages = c("ggplot2", "sf", "sp","raster", "readxl",
             "dplyr", "RandomFields","RColorBrewer","ggridges", "cowplot")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# set home directory
setwd("C:/Users/6350089/OneDrive - Universiteit Utrecht/Food transition/Scripts/R_spatial")

# Import data -------------
## import spatial data that already clean from qgis
crude_NL <- read_sf("crude_NL.shp")

# import farm data
farm_2024 <- readRDS("farm_2024.rds")
# create suscep3 followed Thomas's parameters + 6 times outdoor
farm_2024$suscep3 <- farm_2024$suscep
summary(as.factor(farm_2024$suscep3)) # this susceptibility from Thomas with 6x for outdoor layer and broilers

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-------------------------------------------

# Kernel, tinf and Qinf -------
## Random Q_inf-----------
numsim <-1000 # no. of iteration
totpoints <- nrow(farm_2024)
set.seed(001)
Q_init_matrix <- matrix(0,nrow=numsim,ncol=totpoints)
for (ii in 1:numsim){
  Q_init_matrix[ii,] <- rexp(totpoints, rate = 1) # random threshold to infection
}

# saveRDS(Q_init_matrix,"Q_init_matrix.rds")
Q_init_matrix <- readRDS("Q_init_matrix.rds")

## Create T_inf from fitted model of within farm model --------------
# function for T_inf of chicken 
t_inf_chick <- function(log_size) {
  mean <-1.0516315 + 0.3599556*log_size + 0.0602508*log_size^2
  var <-  0.1954552 # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# function for T_inf of turkey
t_inf_turkey <- function(log_size) {
  mean <-1.1439855 + 0.8948005*log_size + 0.0365499*log_size^2
  var <-  0.2061337 # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# function for T_inf of duck
t_inf_duck <- function(log_size) {
  mean <-2.93114591 + 0.46725553*log_size 
  var <-  0.01038639*log_size # constant variance
  # Calculate shape and rate based on x
  shape <- mean^2/var
  rate <- mean/var
  rgamma(1, shape = shape, rate = rate )}

# Indicate types
chick_id <- which(farm_2024$type%in% c("LAYER","LAYER_BREEDER" ,"BROILER","BROILER_BREEDER"))
turkey_id <- which(farm_2024$type%in% c("MEAT_TURKEY"))
duck_id <- which(farm_2024$type%in% c("MEAT_DUCK", "DUCK_BREEDER"))
# Create matrix of T_inf
T_inf_matrix <- matrix(0,nrow=numsim,ncol=totpoints)
for (ii in 1:numsim){
    T_inf_matrix[ii,chick_id] <- sapply(log(farm_2024$size[chick_id]), FUN= t_inf_chick) 
    T_inf_matrix[ii,turkey_id] <- sapply(log(farm_2024$size[turkey_id]), FUN= t_inf_turkey) 
    T_inf_matrix[ii,duck_id] <- sapply(log(farm_2024$size[duck_id]), FUN= t_inf_duck) 
}

# Check T_inf with plots
#plot(farm_2024$size[chick_id],T_inf_matrix[1,chick_id] ) 
#plot(farm_2024$size[turkey_id],T_inf_matrix[1,turkey_id] ) 
#plot(farm_2024$size[duck_id],T_inf_matrix[1,duck_id] ) 
# It looks okay.
# saveRDS(T_inf_matrix,"T_inf_matrix.rds") save for back up
T_inf_matrix <- readRDS("T_inf_matrix.rds")

## Function transmission kernel --------
# Calibrated parameters
h0 <- 0.0005741044 # calibrated H0
r0  <- 2.5
alpha  <- 2.2
h_kernel <- function(r){h0/(1 + (r/r0)^alpha)} ; # transmission kernel as a function of r
beta<-1

#~~~~~~~~~~~~~~~~~~~~--------------------------------------
# Define index case ----------
# From Egil report: Index farms were selected randomly in a high density poultry area (HDPA) or
# a low density poultry area (LDPA). DPPA farms were random selected from those farms with the 5% percentile highest point density 
# within 1 km from a farm. SPPA were randomly selected among the 75% percentile lowest point density within 1 km.

# So first, we need to start by finding how many farm with radius 1 km of a farm
farm_df <- farm_2024

# Express the coordinates as complex numbers for fast calculation of the euclidean distance
matrix_points <- farm_df[,c("X","Y")]
totpoints <- nrow(matrix_points) # total number of points
colnames(matrix_points) <- c("xcoord","ycoord")
# add a column for the index
Index_points <- c(1:totpoints)
Coord <- (complex(length.out=2,real=matrix_points$xcoord,imaginary=matrix_points$ycoord));
distancematrix <- as.matrix(abs(outer(Coord,Coord,"-")))
distancematrix <- distancematrix/1000 # change to km

# From the distance matrix count the distance between farm that below 1 km for each column
farm_df$point_1km <- colSums(distancematrix  <= 1)
summary(farm_df$point_1km )

# Define if it is in high or low density area
farm_df$index <- NA
farm_df$index[farm_df$point_1km >= quantile(farm_df$point_1km, 0.95)] <- "high" # the farms above 5% are in high density areas
farm_df$index[farm_df$point_1km <= quantile(farm_df$point_1km, 0.75)] <- "low" # the farms below 75% are in low  density areas

# add this to farm 2024
farm_2024$index <-farm_df$index

# check by plot in the map
points_sf <- st_as_sf(farm_df, coords = c("X", "Y"), crs = 28992)

ggplot(data = points_sf) +
  geom_sf(data = crude_NL, fill = "white", color = "black") +
  geom_sf(aes(color = index), size = 1, ) +
 # scale_color_manual(values = c("low" = "green", "high" = "red"))+
    theme_minimal() +
  labs(title = "Farms in high and low density areas",
       x = "Longitude",
       y = "Latitude",
       color = "Density")
#~~~~~~~~~~~~~~~~~~~~--------------------------------------
# Run simulation for baseline with high density index case ----------------
farm_df <- farm_2024

het_matrix <-  as.matrix(abs(outer(farm_df$suscep,farm_df$infect,"*")))

# Express the coordinates as complex numbers for fast calculation of the euclidean distance
matrix_points <- farm_df[,c("X","Y")]
totpoints <- nrow(matrix_points) # total number of points
colnames(matrix_points) <- c("xcoord","ycoord")
# add a column for the index
Index_points <- c(1:totpoints)
Coord <- (complex(length.out=2,real=matrix_points$xcoord,imaginary=matrix_points$ycoord));
distancematrix <- as.matrix(abs(outer(Coord,Coord,"-")))
distancematrix <- distancematrix/1000 # change to km

# create an hazard matrix evaluating for each host j the chance to be infected by host i as a function of distance 
hazardmatrix_before  <- as.matrix(apply(distancematrix,MARGIN=c(1,2),FUN=h_kernel))

# Create matrix for farm size
theta = 7490
# Create farm size matrix
farm_df$mod_size <- 1-exp(-1*farm_df$size/theta)
size_matrix <- as.matrix(abs(outer(farm_df$mod_size,farm_df$mod_size,"*")))
diag(size_matrix) <- 0

numsim <- 1000

# Now 47.5% of layer farms is outdoor and 6.5% of broiler is outdoor, this will be baseline
# if we calculate from the (outdoor broilers+layers)/(total broilers+layers) = (46+421)/1598 = 0.29

result <- list()
source("event.R")


for(i in 1:numsim) {
  
  # Multiplied
  hazardmatrix<- hazardmatrix_before * het_matrix * size_matrix
  diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0
  
  # take T_inf and Q_inf
  T_inf <- T_inf_matrix[i,]   # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
  Q_init <- Q_init_matrix[i,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1
  # Define index case 
  K <- sample(which(farm_df$index=="high"), size = 1) # random index case 
  
  source("InitSim.R") # initialization and setting of the first infected
  while(nrow(Queue)!=0){
    source("Simloop.R")
    result[[i]] <- History
    
  }
  print(i)
  
}


results <- bind_rows(result, .id = "iter")
results$per_od_layer <- 47.5
results$per_od_broiler <- 6.5
results$scen <- "Baseline"
results$cluster <- "Baseline"
results$index <- "high"
saveRDS(results, "result_baseline_high_index.rds") 

# number of infected farms
n_infected <- results  %>% filter(Type_event ==2) %>% group_by(iter) %>% summarise(n=n())
hist(n_infected$n)
summary(n_infected$n)
# outbreak duration
duration <- results %>% group_by(iter) %>% summarise(duration=max(Event_time))
hist(duration$duration)
summary(duration$duration)

# read results from all outdoor
results <- readRDS("C:/Users/6350089/OneDrive - Universiteit Utrecht/Food transition/Scripts/R_spatial/result_2024_data/all_outdoor_high_index.rds")
# number of infected farms
n_infected <- results  %>% filter(Type_event ==2) %>% group_by(iter) %>% summarise(n=n())
hist(n_infected$n)
summary(n_infected$n)
# outbreak duration
duration <- results %>% group_by(iter) %>% summarise(duration=max(Event_time))
hist(duration$duration)
summary(duration$duration)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~--------------------------
# Run simulation for % of transition with high density index case ----------------
farm_df <- farm_2024
# Instead of changing all farms to indoor before randomly sampling, 
# we will randomly decrease or increase from baseline
table(farm_df$type2, farm_df$system) 
46/(666+46) # 6.46% outdoor broiler
421/(465+421) # 47.5% outdoor layer



# Now 47.5% of layer farms is outdoor and 6.5% of broiler is outdoor, this will be baseline
# We reset the farm system that all layer and broiler farms are indoor system
# then 0%,20%,40%,60%,80% and 100% of broilers and layer farms are outdoor (randomly chosen) 

# Get the id of layer and broiler farms beforehand to save time for sampling
layer_id <- farm_df$host_id[farm_df$type2 =="LAYER"]
layer_id_outdoor <- farm_df$host_id[farm_df$type2 =="LAYER" & farm_df$system =="Outdoor"]
broiler_id <- farm_df$host_id[farm_df$type2 =="BROILER"]
broiler_id_outdoor <- farm_df$host_id[farm_df$type2 =="BROILER" & farm_df$system =="Outdoor"]
numsim <- 1000

# Function to random transition of outdoor farms, and change their susceptibility-------
# we can also use this to track the transition farm, by compare them with baseline susceptibility
# unequal susceptibility are the transition farms
per_od_layer <- c(0,25,50,75,100) # to avoid floating number issue, we use integer for %
per_od_broiler <-c(0,50,100)
per_transit <- expand.grid(per_od_layer,per_od_broiler) # combination of these two % transition
#per_transit <- per_transit[2:14,] # remove 0-0 and 1-1 combination since it is all outdoor and all indoor

for(j in 1:nrow(per_transit)){

transit_suscep <- list()
per_transit_layer <- per_transit[j,1]
per_transit_broiler <- per_transit[j,2]

for(i in 1:numsim){
  suscep <- farm_df$suscep
   # For layer farms 
   # in case 0% outdoor layer farms, all layer farms change to indoor
   if (per_transit_layer==0){
     id_change <- farm_df$host_id[farm_df$type2 =="LAYER" & farm_df$system =="Outdoor" ]
     suscep[id_change] <- 1
   }
  
  # in case 25% outdoor layer farms, random half of outdoor layer farms to be indoor
  if (per_transit_layer==25){
    id_change <- sample(farm_df$host_id[farm_df$type2 =="LAYER" & farm_df$system =="Outdoor" ], size = length(layer_id_outdoor)- round(0.25 *length(layer_id)), replace =FALSE)
    suscep[id_change] <- 1
  }
  
  # in case 50%,75% outdoor layer farms, random indoor layer farms to be outdoor
  if (per_transit_layer %in% c(50,75)){
    id_change <- sample(farm_df$host_id[farm_df$type2 =="LAYER" & farm_df$system =="Indoor" ], size = round(per_transit_layer *length(layer_id)/100)- length(layer_id_outdoor), replace =FALSE)
    suscep[id_change] <- 6.3
  }
  
  # in case 100% outdoor layer farms, random indoor layer farms to be outdoor
  if (per_transit_layer == 100){
    id_change <- farm_df$host_id[farm_df$type2 =="LAYER" & farm_df$system =="Indoor" ]
    suscep[id_change] <- 6.3
  }
  
  #For broiler farms
  # in case 0%, all broiler farms change to indoor
  if (per_transit_broiler==0){
    id_change <- farm_df$host_id[farm_df$type2 =="BROILER" & farm_df$system =="Outdoor" ]
    suscep[id_change] <- 0.134
  }
  
  # in case 25% , random half of outdoor layer farms to be indoor
  if (per_transit_broiler > 0 & per_transit_broiler < 100){
    id_change <- sample(farm_df$host_id[farm_df$type2 =="BROILER" & farm_df$system =="Indoor" ], size = round(per_transit_broiler *length(broiler_id)/100)- length(broiler_id_outdoor), replace =FALSE)
    suscep[id_change] <- 0.8442
  }
  
 
  # in case 100% layer farms, random indoor layer farms to be outdoor
  if (per_transit_broiler == 100){
    id_change <- farm_df$host_id[farm_df$type2 =="BROILER" & farm_df$system =="Indoor" ]
    suscep[id_change] <- 0.8442
  }
  
  transit_suscep[[i]] <- suscep
}

saveRDS(transit_suscep, paste0('scen1','_layer',per_transit[j,1],'_broiler',per_transit[j,2],'.rds'))

}

#Check
# change to character before using == to avoid floating number problem
per_transit_layer =50
per_transit_broiler = 75

farm_df$suscep[farm_df$type2 =="LAYER" & farm_df$system =="Outdoor" ]
farm_df$suscep[farm_df$type2 =="LAYER" & farm_df$system =="Indoor" ]


sum(as.character(suscep)=="0.134") 
sum(as.character(suscep)=="0.8442")/712 # total broiler = 712

sum(as.character(suscep)=="1") 
sum(as.character(suscep)=="6.3")/886 # total layer = 886


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-----------------------
# Prepare matrix before running model --------------
farm_df <- farm_2024

# Calculate hazard matrix and size matrix------------------
# Express the coordinates as complex numbers for fast calculation of the euclidean distance
matrix_points <- farm_df[,c("X","Y")]
totpoints <- nrow(matrix_points) # total number of points
colnames(matrix_points) <- c("xcoord","ycoord")
# add a column for the index
Index_points <- c(1:totpoints)
Coord <- (complex(length.out=2,real=matrix_points$xcoord,imaginary=matrix_points$ycoord));
distancematrix <- as.matrix(abs(outer(Coord,Coord,"-")))
distancematrix <- distancematrix/1000 # change to km

# create an hazard matrix evaluating for each host j the chance to be infected by host i as a function of distance 
hazardmatrix_before  <- as.matrix(apply(distancematrix,MARGIN=c(1,2),FUN=h_kernel))

# Create matrix for farm size
theta = 7490
# Create farm size matrix
farm_df$mod_size <- 1-exp(-1*farm_df$size/theta)
size_matrix <- as.matrix(abs(outer(farm_df$mod_size,farm_df$mod_size,"*")))
diag(size_matrix) <- 0


# Run model --------------
files <- list.files(path ="C:/Users/6350089/OneDrive - Universiteit Utrecht/Food transition/Scripts/R_spatial/result_2024_data/Scen1.1 random_extensive_high")
files <- files[1:15]

for(j in 1:15){

result <- list()
source("event.R")
infect <- farm_df$infect
suscep <- readRDS(paste0("C:/Users/6350089/OneDrive - Universiteit Utrecht/Food transition/Scripts/R_spatial/result_2024_data/Scen1.1 random_extensive_high/",files[j]))


for(i in 1:numsim) {
  # read suscep from randomly transit matrix
  het_matrix <-  as.matrix(abs(outer(suscep[[i]],farm_df$infect,"*")))
  
  # Multiplied
  hazardmatrix<- hazardmatrix_before * het_matrix * size_matrix
  diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0
  
  # take T_inf and Q_inf
  T_inf <- T_inf_matrix[i,]   # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
  Q_init <- Q_init_matrix[i,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1
  # Define index case 
  K <- sample(which(farm_df$index=="high"), size = 1) # random index case 
  
  source("InitSim.R") # initialization and setting of the first infected
  while(nrow(Queue)!=0){
    source("Simloop.R")
    result[[i]] <- History
    
  }
  print(i)
}
# Summarise outcomes -------
# bind list
results <- bind_rows(result, .id = "iter")
per_od_layer <- as.numeric(sub(".*_layer([0-9]+)_.*", "\\1", files[j]))
per_od_broiler  <- as.numeric(sub(".*_broiler([0-9]+)\\.rds", "\\1", files[j]))
results$per_od_layer <- per_od_layer
results$per_od_broiler <- per_od_broiler 
results$scen <- "1.1"
results$cluster <- "random"
results$index <- "high"
saveRDS(results, paste0("Sim1.1_layer",per_od_layer,"_broiler",per_od_broiler,"_high",".rds" )) 
gc()
}





#saveRDS(results, "result_extensive_rand_1.rds")
# number of infected farms
n_infected <- results %>% filter(Type_event ==2) %>% group_by(iter) %>% summarise(n=n())
hist(n_infected$n)
summary(n_infected$n)
# outbreak duration
duration <- results %>% group_by(iter) %>% summarise(duration=max(Event_time))
hist(duration$duration)
summary(duration$duration)

# Plot results in R file plot_all_results
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~--------------------------------

# Calculate Ri -----------
# Function to calculate the reproduction number of farms ---------
# Follow the Ri calculation from Beninca's paper in page 7 Eq6 
# Infectious period from gamma distibution
# Indicate types
chick_id <- which(farm_2024$type%in% c("LAYER","LAYER_BREEDER" ,"BROILER","BROILER_BREEDER"))
turkey_id <- which(farm_2024$type%in% c("MEAT_TURKEY"))
duck_id <- which(farm_2024$type%in% c("MEAT_DUCK", "DUCK_BREEDER"))

# First we calculate mean and variance for infectious duration of each farms
# for mean
farm_2024$Tinf_mean <- NA
farm_2024$Tinf_mean[chick_id] <- 1.0516315 + 0.3599556*log(farm_2024$size[chick_id]) + 0.0602508*log(farm_2024$size[chick_id])^2
farm_2024$Tinf_mean[turkey_id] <- 1.1439855 + 0.8948005*log(farm_2024$size[turkey_id]) + 0.0365499*log(farm_2024$size[turkey_id])^2
farm_2024$Tinf_mean[duck_id] <- 2.93114591 + 0.46725553*log(farm_2024$size[duck_id])
summary(farm_2024$Tinf_mean)

# for variance
farm_2024$Tinf_var <- NA
farm_2024$Tinf_var[chick_id] <- 0.1954552
farm_2024$Tinf_var[turkey_id] <- 0.2061337
farm_2024$Tinf_var[duck_id] <- 0.01038639*log(farm_2024$size[duck_id])
summary(farm_2024$Tinf_var)

# From Boendar https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030071
# if the infectious periods arise from a common gamma distribution with 
# shape parameter = c and scale parameter = T/c, where T = the mean and 1/sqrt(c) = variance of the infectious period probability distribution
shape <- farm_2024$Tinf_mean^2/farm_2024$Tinf_var # calculate shape of gamma distribution = mean^2/var
scale <- farm_2024$Tinf_var/farm_2024$Tinf_mean # scale = 1/rate; rate  = mean/var

farm_df <- farm_2024
matrix_points <- data.frame(farm_df$X, farm_df$Y)
totpoints <- nrow(matrix_points) # total number of points
colnames(matrix_points) <- c("xcoord","ycoord")


# Express the coordinates as complex numbers for fast calculation of the euclidean distance
Coord <- (complex(length.out=2,real=matrix_points$xcoord,imaginary=matrix_points$ycoord));
distancematrix <- as.matrix(abs(outer(Coord,Coord,"-")))
distancematrix <- distancematrix/1000 # change to km
# calculate hazard matrix
hazardmatrix_before  <- as.matrix(apply(distancematrix,MARGIN=c(1,2),FUN=h_kernel))

# Ri for baseline
# het matrix
suscep <- farm_df$suscep
infect <- farm_df$infect

het_matrix <-  as.matrix(abs(outer(suscep,infect,"*")))

# Multiplied
hazardmatrix<- hazardmatrix_before * het_matrix
diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0

# Calculate kernel specific to shape and rate of individual farms
# Ri of each farm is the sum of column of kernel
Ri_matrix <- matrix(NA, ncol = ncol(hazardmatrix), nrow = nrow(hazardmatrix))

for(i in 1:ncol(hazardmatrix)){
  Ri_matrix[,i] <- 1-((1/(1+shape[i]*hazardmatrix[,i]))^scale[i])
}
#rgamma(1, shape = shape[1], rate = rate[1] )

diag(Ri_matrix) <- 0
Ri <- apply(Ri_matrix, MARGIN=2,FUN=sum)  
summary(Ri)


# map plot Ri -----
farm_df$Ri <- Ri
mid_Ri <-(1-min(farm_df$Ri))/(max(farm_df$Ri)-min(farm_df$Ri)) # standardize to see value of Ri = 1 to set the mid range for coulour bar : (value-min)/(max-min) 
# Convert to an sf object
points_sf <- st_as_sf(farm_df, coords = c("X", "Y"), crs = 28992)

ggplot(data = points_sf) +
  geom_sf(data = crude_NL, fill = "white", color = "black") +
  geom_sf(aes(color = Ri, shape = type), size = 1, ) +
  scale_color_gradientn(colours = c("green", "red","red"), 
                        values = c(0, mid_Ri, 1))+
  theme_minimal() +
  labs(title = "100% outdoor layer farms",
       x = "Longitude",
       y = "Latitude",
       color = "Ri")


# 0% outdoor
per_outdoor <- 0
farm_df_new <- farm_df
# het matrix
a <- sample(which(farm_df$type =="LAYER"), size = per_outdoor *length(which(farm_df$type =="LAYER")), replace =FALSE)
farm_df_new$system[a] <- "OUTDOOR-LAYER" # change to outdoor layer
suscep <- farm_df_new$suscep
infect <- farm_df_new$infect
suscep[a] <- 6.3 # change suscep of outdoor layer to 6.3

het_matrix <-  as.matrix(abs(outer(suscep,infect,"*")))
# Multiplied
hazardmatrix<- hazardmatrix_before * het_matrix
diag(hazardmatrix) <- 0 # because the chance of infecting itself is 0

Ri_matrix <- as.matrix(apply(hazardmatrix, MARGIN=c(1,2),FUN= function(x) 1-((1/(1+w*x))^c) ))
diag(Ri_matrix) <- 0
Ri <- apply(Ri_matrix, MARGIN=2,FUN=sum)  
summary(Ri)

summary(as.factor(farm_df_new$system))

# Another Ri map for other %per_outdoor
farm_df_new$Ri <- Ri
mid_Ri <-(1-min(farm_df_new$Ri))/(max(farm_df_new$Ri)-min(farm_df_new$Ri)) # standatize to see value of Ri = 1 to set the mid range for coulour bar : (value-min)/(max-min) 
summary(farm_df_new$Ri)
# Convert to an sf object
points_sf <- st_as_sf(farm_df_new, coords = c("X", "Y"), crs = 28992)

R1 <- ggplot(data = points_sf) +
  geom_sf(data = crude_NL, fill = "white", color = "black") +
  geom_sf(aes(color = Ri, shape = type), size = 1, ) +
  scale_color_gradientn(colours = c("green", "red","red"), 
                        values = c(0, mid_Ri, 1))+
  theme_minimal() +
  labs(title = "100% indoor layer farms",
       x = "Longitude",
       y = "Latitude",
       color = "Ri")

# Combine plot
plot_grid(R1, R2, labels = c('A', 'B'), label_size = 12, ncol = 2)

# Reproduction number in square km


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~--------------------------------
# Probability of HPAI introduction ----------
# environment introduction indoor layer from Gonzales et al., 2013 DOI:10.1111/j.1750-2659.2012.00348.x
# beta1 is RR from Bouwstra,2017 https://doi.org/10.3201/eid2309.170276
beta1 <- data.frame(system =c( "INDOOR-LAYER","OUTDOOR-LAYER", "LAYER-BREEDER", "INDOOR-BROILER", "BROILER-BREEDER",
                               "MEAT-TURKEY","TYRKEY-BREEDER", "MEAT-DUCK",  "DUCK-BREEDER"),
                    beta1 = c(1, 6.3, 0.5, 0.2, 0.4, 12,11.3, 39.5, 25.5))

beta2 = 0.8 # relative risk from Baustrow, 2017

# For now we will use medium water ways until we got the other environmental factors from Jose
# Import water ways-----
medium_water <- read_sf("water3-6m.shp")
medium_water <- sf::st_transform (medium_water, crs = 28992)

# Calculate the shortest distance between farms and medium eater ways-----
farms <- st_as_sf(x = as.data.frame(farm_df[,c("X","Y")]), coords = c("X", "Y"), crs = 28992 )
dist_medwater <- apply(st_distance(farms, medium_water,by_element = FALSE),1,min)
# for other environmental factor check other R script to extract environmental factors
farm_df$dist_medwater <- dist_medwater
# take log
farm_df$ln_dist_medwater <- log(farm_df$dist_medwater)
# calculate beta0 to calibrate it to get the rate of indoor layer farms at average ln medium water distance = 3.5*10^-4 per month
# from Gonzales 2012,The rate of introduction of LPAI to indoor-layer farms was 3.5*10^-4 per month at average log_distance to water way.
# 3.5*10^-4 = exp(log(beta0)+log(beta1)*6.579121)
# log(3.5*10^4) = log(beta0)+log(beta1)*6.579121)
# log(3.5*10^4)-log(beta1)*6.579121) = log(beta0)
# beta0 = exp(log(3.5*10^4)-log(beta1)*6.579121))

beta0 <- exp(log(3.5*10^-4)-(log(0.8)* mean(farm_df$ln_dist_medwater )))/30.44 # /30.44 to change to per day

# Calculate risk of environmental introduction -----
# 0% outdoor
per_outdoor <- 1
farm_df_new <- farm_df
# het matrix
a <- sample(which(farm_df$type =="LAYER"), size = per_outdoor *length(which(farm_df$type =="LAYER")), replace =FALSE)
farm_df_new$system[a] <- "OUTDOOR-LAYER" # change to outdoor layer
# join beta 1 
farm_df_new <-left_join(farm_df_new, beta1, by = "system")


#  probability = 1-exp(-rate*t)
farm_df_new$rate_day <- exp((log(beta0)+ log(farm_df_new$beta1)+ (log(beta2)* log(farm_df_new$dist_medwater)) ))
farm_df_new$prob_year <- 1-exp(-(farm_df_new$rate_day*365))
hist(farm_df_new$prob_year, main="All layer farms become outdoor system",
     xlab="Probability of HPAI introduction per year")
# standardize to get red color at prob = 0.5
# mid_prob <-(0.5-min(farm_df_new$prob_year))/(max(farm_df_new$prob_year)-min(farm_df_new$prob_year))
# Convert to an sf object
points_sf <- st_as_sf(farm_df_new, coords = c("X", "Y"), crs = 28992)

P1 <- ggplot(data = points_sf) +
  geom_sf(data = crude_NL, fill = "white", color = "black") +
  geom_sf(aes(color = prob_year, shape = type ), size = 1) +
  scale_color_gradient(low = "green", high = "red") +
  theme_minimal() +
  labs(title = "100% indoor layer farms",
       x = "Longitude",
       y = "Latitude",
       color = "Yearly probability of introduction")

# prob for 100% outdoor layers 
# Calculate risk of environmental introduction -----
# 100% outdoor
per_outdoor <- 0
farm_df_new <- farm_df
# het matrix
a <- sample(which(farm_df$type =="LAYER"), size = per_outdoor *length(which(farm_df$type =="LAYER")), replace =FALSE)
farm_df_new$system[a] <- "OUTDOOR-LAYER" # change to outdoor layer
# join beta 1 
farm_df_new <-left_join(farm_df_new, beta1, by = "system")

#  probability = 1-exp(-rate*t)
farm_df_new$rate_day <- exp((log(beta0)+ log(farm_df_new$beta1)+ (log(beta2)* log(farm_df_new$dist_medwater)) ))
farm_df_new$prob_year <- 1-exp(-(farm_df_new$rate_day*365))

# standardize to get red color at prob = 0.5
# mid_prob <-(0.5-min(farm_df_new$prob_year))/(max(farm_df_new$prob_year)-min(farm_df_new$prob_year))
# Convert to an sf object
points_sf <- st_as_sf(farm_df_new, coords = c("X", "Y"), crs = 28992)

P2 <- ggplot(data = points_sf) +
  geom_sf(data = crude_NL, fill = "white", color = "black") +
  geom_sf(aes(color = prob_year, shape = type ), size = 1) +
  scale_color_gradient(low = "green", high = "red") +
  theme_minimal() +
  labs(title = "100% outdoor layer farms",
       x = "Longitude",
       y = "Latitude",
       color = "Yearly probability of introduction")

# Combine plot
cowplot::plot_grid(P1, P2, labels = c('A', 'B'), label_size = 12, ncol = 2)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~--------------------------------
# Combine results together------
files <- list.files(path = ("C:/Users/6350089/OneDrive - Universiteit Utrecht/Food transition/Scripts/R_spatial/result_2024_data/Scen1.1 random_extensive_high/sim_result"))
setwd("C:/Users/6350089/OneDrive - Universiteit Utrecht/Food transition/Scripts/R_spatial/result_2024_data/Scen1.1 random_extensive_high/sim_result")
sum_scen1.1 <- lapply(files, readRDS)

# combine list
sum_scen1.1 <- dplyr::bind_rows(sum_scen1.1)


# summary n_infected by iter and per_outdoor
n_infected <- sum_scen1.1 %>% filter(Type_event ==2) %>% group_by(per_od_layer, per_od_broiler, cluster ,iter ) %>% summarise(n=n())

# number of outbreaks that n infected farm >= 10
sum_n_infected <- n_infected %>% group_by(per_od_layer, per_od_broiler, cluster) %>% summarise(major_outbreak = sum(n>=10),
                                                                      max_outbreak_size = max(n),
                                                                      mean_outbreak_size = mean(n))

# add per_outdoor ==0 and per_outdoor ==1 to cluster 1
a <- rbind(sum_n_infected[sum_n_infected$per_outdoor ==0,], sum_n_infected[sum_n_infected$per_outdoor ==1,])
a$cluster <- 1
sum_n_infected <- rbind(sum_n_infected,a)
sum_n_infected$per_outdoor <- sum_n_infected$per_outdoor*100
sum_n_infected$per_major_outbreak <- sum_n_infected$major_outbreak/2000*100 # change to percentage
sum_n_infected$cluster[sum_n_infected$cluster==0] <- "Random"
sum_n_infected$cluster[sum_n_infected$cluster==1] <- "Cluster"

# line plot for major outbreaks
A <- sum_n_infected %>%
  ggplot( ) +
  geom_line( aes(x=per_outdoor, y=per_major_outbreak, color=as.factor(cluster), size =1)) +
  geom_point(aes(x=per_outdoor, y=per_major_outbreak, color=as.factor(cluster), size =5, alpha = 0.5)) + 
  scale_color_manual(values = c('Red',  'Blue'))+
  ylim(0,100)+
  theme_minimal_grid()+ xlab("Percentage of outdoor transition") + ylab("Percentage of major outbreaks")+
  guides(alpha=FALSE, size = FALSE,color = FALSE) 
# + ggtitle("Percentage of major outbreaks (final size >= 10 infected farms)")

# line plot for max outbreak
B <- sum_n_infected %>%
  ggplot( ) +
  geom_line( aes(x=per_outdoor, y=max_outbreak_size, color=as.factor(cluster), size =1)) +
  geom_point(aes(x=per_outdoor, y=max_outbreak_size, color=as.factor(cluster), size =5, alpha = 0.5)) + 
  scale_color_manual(values = c('Red',  'Blue'))+
  ylim(0,1784)+
  theme_minimal_grid()+ xlab("Percentage of outdoor transition") + ylab("Maximum outbreak size")+
  guides(alpha=FALSE, size = FALSE,color = FALSE)
# + ggtitle("Maximum outbreak size")+ 

C <- get_legend(sum_n_infected %>%
                   ggplot( ) +
                   geom_line( aes(x=per_outdoor, y=per_major_outbreak, color=as.factor(cluster), size =1)) +
                   geom_point(aes(x=per_outdoor, y=per_major_outbreak, color=as.factor(cluster), size =5, alpha = 0.5)) + 
                   scale_color_manual(values = c('Red',  'Blue'))+
                   ylim(0,100)+
                   theme_bw()+ xlab("Percentage of outdoor transition") + ylab("Percentage of major outbreaks")+
                   ggtitle("Percentage of major outbreaks (final size >= 10 infected farms)")+ 
                   guides(alpha=FALSE, size = FALSE) + guides(color=guide_legend(title="Transition pattern")))
# Combine plot
plot_grid(A, B, C ,labels = c('A', 'B', ''), label_size = 12, ncol = 3,rel_widths = c(1,1,0.5))
?plot_grid
