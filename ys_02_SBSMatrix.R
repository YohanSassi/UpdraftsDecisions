#___________________________________________________________________________________________________________________________________#
###
###       
###       Project: Explore the balance between personal and social information in thermal updraft selection
###
###        Data come from captive vulture from the Rocher des Aigles, Rocamadour, France experiencing free-flights
###        in a 120m deep canyon.
### 
###       !!! Script used to create Social Bond Strength matrix between individuals from in-aviary camera trap !!!
###
###       Author: Yohan Sassi 
###               yohan.sassi@cefe.cnrs.fr
###
#___________________________________________________________________________________________________________________________________#

rm(list = ls())




#____________________________________________________#
###
#### Required functions & packages ####
###
#____________________________________________________#

# Package needed
library(stringr)
library(ggplot2)
library(cowplot)
library(asnipe)
library(raster)




### function to convert distance in pixel to meters
convert_dist <- function (cam, dist_px){
  if (cam == "CT01"){
    dist_m = (3.10 * dist_px)/dist_CT01
  }
  
  if (cam == "CT02"){
    dist_m = (3.10 * dist_px)/dist_CT02
  }
  
  if (cam == "CT03"){
    dist_m = (3.10 * dist_px)/dist_CT03
  }
  return(dist_m)
}

### function to find the match between pair numbers and individual names
Numextract <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}


### Vulture names depending of the year
vultureNames <- vector("list")

vultureNames[[1]] <- list(year = "2021", 
                          names = c("kazimir","henry","gregoire","leon","kirikou","bulma"),
                          sexes = c("M","F","M","M","F","F"),
                          ageYear = c(2014,1997,2011,2016,2014,2014)) # For year1 = 2020
vultureNames[[2]] <- list(year = "2022", 
                          names = c("kazimir","mathilda","hercule","leon","kirikou","bulma"),
                          sexes = c("M","F","F","M","F","F"),
                          ageYear = c(2014,2018,2012,2016,2014,2014)) # For year2 = 2021




### wing length, allow to order by within or out the distance off the wing span
## Could it be 1.55m ? - half of the perch? The mode of the distribution of distance
#dist_wing <- read.csv("./Data/wing_length.csv",h=T,sep=",")
dist_cam <- read.csv("./Data/distPerches.csv",h=T,sep=";")






#____________________________________________________#
###
####  Social bonds among the group members ####
###
#____________________________________________________#

# Create a unique file containing the picture data for a full years
# Create in same time the association data, that ill be used in asnipe after for social network analysis

year <- c("2021", "2022")
allFiles <- list.files("./Data/Rocamadour_SBSData", ".csv" , all.files = F, full.names = T) # All files containing inter-individual distances on perches
distThreshold <- 1.55

for (y in 1:length(year)){
  
  # Adapt camera distances px to meters
  dist_camYear <- dist_cam[str_detect(unique(dist_cam$pic_name),year[y]),]
  
  dist_CT01 <- as.integer(dist_camYear[1,11])
  dist_CT02 <- as.integer(dist_camYear[2,11])
  dist_CT03 <- as.integer(dist_camYear[3,11])
  
  
  # Files selected
  wantFiles <- allFiles[which(grepl(paste(year[y],".csv", sep=""), allFiles, fixed = TRUE))] 
  
  ## -- Create the full dataset for the specific year
  for (i in 1:length(wantFiles)){
    
    # load a file
    wantdata_sbs <- read.csv(wantFiles[i],h=T,sep=",")
    
    # Add a column corresponding to the camera number
    hyfenindices <- str_locate_all(pattern ='_', wantFiles[i]) # locate the _ in the name
    camName <- substr(wantFiles[i], 
                      as.integer(hyfenindices[[1]][3,1])+1, 
                      as.integer(hyfenindices[[1]][4,1])-1) # Extract the camera name
    
    wantdata_sbs$cam <- camName # Add a column with the camera name
    
    # Concatenate fill one after the other
    if (i==1){
      data_sbs <- wantdata_sbs
    } else {
      data_sbs <- rbind(data_sbs, wantdata_sbs)
    }
  }
  
  ## remove potential mistakes where vult1 == vult2
  data_sbs <- data_sbs[-which(data_sbs$name_vult1 == data_sbs$name_vult2),] 
  
  
  ## convert distances in meters & associate number to vulture name to identify unique pairs --> kazimir_bulma = bulma_kazimir
  match <- data.frame(nb = 1:6, name_vult = vultureNames[[y]]$names) # year specific as group are different
  
  for (r in 1:nrow(data_sbs)){
    
    # distance conversion
    data_sbs$dist_m[r] <- convert_dist(as.character(data_sbs$cam[r]),data_sbs$dist[r])
    
    # associate number to vulture name to identify unique pairs
    data_sbs$num_vult1[r] <- match[which(match$name_vult == data_sbs$name_vult1[r]),1]
    data_sbs$num_vult2[r] <- match[which(match$name_vult == data_sbs$name_vult2[r]),1]
    
  }
  
  # Associate unique pairs numbers
  data_sbs$pair_id <- paste(apply(cbind(data_sbs$num_vult1,data_sbs$num_vult2),1,min),
                            apply(cbind(data_sbs$num_vult1,data_sbs$num_vult2),1,max),sep="_")
  
  # find the match between pairs numbers and pair names
  for (l in 1:nrow(data_sbs)){
    data_sbs$pair_name[l] <- paste(match[which(match$nb == as.integer(Numextract(data_sbs$pair_id[l]))[1]),2],
                                match[which(match$nb == as.integer(Numextract(data_sbs$pair_id[l]))[2]),2],sep="_")
  }
  

  ## --- PREPARE THE ASSOCIATION DATA
  
  ## Create data format needed for asnipe - association between individuals (groups) in pictures (days)
  asso_data <- data.frame(ID = character(0), GROUP = integer(0), DAY = integer(0)) 
  c = 0 # counter for days
  
  for (p in 1:length(unique(data_sbs$pic_name))) { 
    
    # Subsample the data on a specific picture
    sub_pic <- data_sbs[which(data_sbs$pic_name == unique(data_sbs$pic_name)[p]),]
    
    # Is there associations?
    sub_pic_asso <- sub_pic[sub_pic$dist_m <= distThreshold,]
    
    # If there is association between some individuals
    if (nrow(sub_pic_asso) != 0) {
      c = c+1
      
      for (n in 1:nrow(sub_pic_asso)){ # Treat row by row
        
        # Extract identity of the dyad
        ID = c(sub_pic_asso$name_vult1[n],sub_pic_asso$name_vult2[n])
        
        # Create a group number, which should be next number of previous identifies associations
        if (nrow(asso_data) != 0){
          GROUP = rep(max(asso_data$GROUP)+1,2)
        } else {
          GROUP = rep(1,2)
        }
        
        # Give a day number, which is actually here a number by pictures
        DAY = rep(c,2)
        
        # Create a dataframe out of that
        dat <- data.frame(ID,GROUP,DAY) 
        
        # Paste it into the full one
        asso_data <- rbind(asso_data,dat)
      }
    }
  }
  
  
  ## Define the sampling period
  SPs_asso <- get_sampling_periods(asso_data[,c(1,2)],
                                   asso_data[,3],1,data_format="individuals")
  
  ## Define the occurences of each individuals 
  occurs_asso <- get_sampling_periods(asso_data[,c(1,2)],
                                      asso_data[,3],1,data_format="individuals", return="occ")
  
  # Save association data, sampling periods and occurrences for building the social network
  # save(data_sbs,
  #      asso_data,
  #      SPs_asso,
  #      occurs_asso,
  #      file=paste("./Output/Files/AssociationData_",year[y],".RData",sep=""))

}



## Check the distribution of distances between individuals - with threshold of 1.55m shown with red line
load("~/Documents/PhD/Manip Rocamadour/Project_RedAggr/Output/Files/AssociationData_2021.RData")
distData_2021 <- data_sbs

load("~/Documents/PhD/Manip Rocamadour/Project_RedAggr/Output/Files/AssociationData_2022.RData")
distData_2022 <- data_sbs

allDistData <- rbind(distData_2021[,2:18], distData_2022)


ggplot(data = allDistData, aes(x = dist_m)) + #data_sbs
  geom_histogram(binwidth = 0.1, col = "grey75", fill = "grey75") +
  geom_vline(xintercept = 1.50, linetype="dashed", color = "black", linewidth=0.6) + 
  theme_light() +
  background_grid(major = "xy", minor = "xy", color.major = "grey92", color.minor = "grey95") +
  xlab("Inter-individual distances on perch [m]") +
  ylab("Occurences")
  





#____________________________________________________#
###
####  Association matrix  ####
###
#____________________________________________________#

# Parameters to run this part independently
year <- c("2021", "2022")
y = 1

# Load required file
asso_files <- list.files("./Output/Files", "AssociationData_" , all.files = F, full.names = T) # All association files
load(asso_files[which(grepl(paste(year[y],".RData", sep=""), asso_files, fixed = TRUE))] ) # Load the appropriate file corresponding to the year


## Estimate network from asnipe package
SBSMatrix_estimated <- get_network(association_data = SPs_asso,
                                   data_format = "SP",
                                   association_index = "SRI")


## Reorganize matrix to match with sex & species differences matrix after
SBSMatrix_estimatedreord <- matrix(nrow = 6, ncol = 6, dimnames = list(vultureNames[[y]]$names,vultureNames[[y]]$names))
for (n in 1:6){
  # For rows
  ind_1 <- vultureNames[[y]]$names[n]

  for (j in 1:6){
    # for columns
    ind_2 <- vultureNames[[y]]$names[j]

    if (ind_1 == ind_2){
      SBSMatrix_estimatedreord[n,j] <- 0

    } else {
      SBSMatrix_estimatedreord[n,j] <- SBSMatrix_estimated[ind_1,ind_2]
    }
  }
}


# Save
save(SBSMatrix_estimatedreord,
     file=paste("./Output/Files/SBSMatrix_",year[y],".RData",sep=""))



#____________________________________________________#
###
####  Is there biases in its estimate  ####
###
#____________________________________________________#

# Try to see if some bias in the estimation of social relationships due to species / sex differences
# Otherwise generalized association index are needed

## Create a t x N x N matrix where t is the number of permutations
SBSMatrix_random <- network_permutation(association_data = SPs_asso,
                                        data_format = "SP",
                                        association_index = "SRI",
                                        permutations = 1000,
                                        trialSwap = TRUE
)


# Create matrix of species differences where 1 if same species in the int, 0 otherwise
# Create matrix of sex association where 1 if same sex in the int, 0 otherwise
sexMatrix <- matrix(nrow = 6, ncol = 6, dimnames = list(vultureNames[[y]]$names,vultureNames[[y]]$names))
speciesMatrix <- matrix(nrow = 6, ncol = 6, dimnames = list(vultureNames[[y]]$names,vultureNames[[y]]$names))
ageMatrix <- matrix(nrow = 6, ncol = 6, dimnames = list(vultureNames[[y]]$names,vultureNames[[y]]$names))

for (n in 1:6){
  # For rows
  ind_1 <- vultureNames[[y]]$names[n]
  
  for (j in 1:6){
    # for columns
    ind_2 <- vultureNames[[y]]$names[j]
    
    if (ind_1 == ind_2){
      sexMatrix[n,j] <- 0
      speciesMatrix[n,j] <- 1
      ageMatrix[n,j] <- 0 
      
    } else {
      sexMatrix[n,j] <- ifelse(vultureNames[[y]]$sexes[n] == vultureNames[[y]]$sexes[j],0,1)
      ageMatrix[n,j] <- abs(vultureNames[[y]]$ageYear[n] - vultureNames[[y]]$ageYear[j])
        
      if (ind_1 == "kirikou" | ind_2 == "kirikou") {
        speciesMatrix[n,j] <- 0
        
      } else {
        speciesMatrix[n,j] <- 1
        
      }
    }
  }
} 

# save(sexMatrix,
#      speciesMatrix,
#      ageMatrix,
#      file=paste("./Output/Files/ControlMatrix_",year[y],".RData",sep=""))


# Estimate if species / sex differences are important using the MRQAP with data permutation
mrqap.custom.null(SBSMatrix_estimatedreord ~ speciesMatrix + sexMatrix,
                  random.y = SBSMatrix_random,
                  directed = "undirected",
                  diagonal = FALSE,
                  test.statistic = "beta")



#____________________________________________________#
###
####  Is there preferred/Avoided associations  ####
###
#____________________________________________________#
## Is there more preferred/avoided association than by random in the network?
# To do that we will compare coefficient of variation of the association index in the observed network
# against the one based on randomization of data

## Threshold to know if strong affiliates - based on twice the mean of the SRI association ratio (Whitehead 2008)
mean(SBSMatrix_estimatedreord[upper.tri(SBSMatrix_estimatedreord)]) * 2


## CV of observed network
cvObs = cv(SBSMatrix_estimatedreord)


# Compare observed data CV and randomize (x1000) CV to see if preferred/avoided association more than random
nBoot = 1000
cvs_rand <- rep(NA,nBoot)
network_perm = list(SBSMatrix_estimated,SPs_asso)
set.seed(12345)

for (b in 1:nBoot){
  
  network_perm <- network_swap(association_data = network_perm[[2]],
                      data_format = "SP",
                      association_index = "SRI",
                      association_matrix = network_perm[[1]],
                      swap = 1
  )
  
  cvs_rand[b] <- cv(network_perm[[1]])
  
}

# P - value which is the number of time the observed cv is lower than randomized one / number of bootstrap
p = length(cvs_rand[cvs_rand >= cvObs]) / nBoot
p

# Mean + 95% CI around randomised CV
c(mean(cvs_rand),quantile(cvs_rand,c(0.025,0.975)))

# Plot
plot(cvs_rand,pch=20,cex=0.5) 
points(0,cv(SBSMatrix_estimatedreord),cex=1,pch=20,col="red")




#---
## for each dyads, is a significantly preferred/avoided association?
# Compare dyad's observed SRI against randomized ones (x1000) to see if preferred/avoided association
n.perms <- 1000
P.upper <- P.lower <- matrix(NA, nrow=nrow(SBSMatrix_estimated),ncol=ncol(SBSMatrix_estimated),dimnames = dimnames(SBSMatrix_estimated))

for (i in 1:nrow(SBSMatrix_estimated)) {
  for (j in 1:ncol(SBSMatrix_estimated)) {
    P.upper[i,j] <- sum(SBSMatrix_estimated[i,j] <= SBSMatrix_random[,i,j])/n.perms
    P.lower[i,j] <- sum(SBSMatrix_estimated[i,j] >= SBSMatrix_random[,i,j])/n.perms
  }
}


## Reorganize matrix to match with sex & species differences matrix after
P.upper_reord <- matrix(nrow = 6, ncol = 6, dimnames = list(vultureNames[[y]]$names,vultureNames[[y]]$names))
P.lower_reord <- matrix(nrow = 6, ncol = 6, dimnames = list(vultureNames[[y]]$names,vultureNames[[y]]$names))

for (n in 1:6){
  # For rows
  ind_1 <- vultureNames[[y]]$names[n]
  
  for (j in 1:6){
    # for columns
    ind_2 <- vultureNames[[y]]$names[j]
    
    if (ind_1 == ind_2){
      P.upper_reord[n,j] <- 0
      P.lower_reord[n,j] <- 0
      
    } else {
      P.upper_reord[n,j] <- P.upper[ind_1,ind_2]
      P.lower_reord[n,j] <- P.lower[ind_1,ind_2]
      
    }
  }
}

P.upper_reord # if significant -> preferred association
P.lower_reord # if significant -> avoided association















