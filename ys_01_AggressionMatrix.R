#________________________________________________________________________________________________________________________________________________
###
###       
###       Project: Explore the balance between personal and social information in thermal updraft selection
###
###        Data come from captive vulture from the Rocher des Aigles, Rocamadour, France experiencing free-flights
###        in a 120m deep canyon.
### 
###       !!! Script used to create aggression matrix between individuals during recorded feeding events !!!
###
###       Author: Yohan Sassi 
###               yohan.sassi@cefe.cnrs.fr
###
#________________________________________________________________________________________________________________________________________________


rm(list = ls())


#____________________________________________________#
###
##### Packages & functions needed  ####
###
#____________________________________________________#

library(aniDom)
library(stringr)
library(dplyr)





#______________________________________________________________#
##
###
#### Hierarchy based on randomised Elo ranks ####
###
##
#______________________________________________________________#

### TEST HIERARCHY ESTIMATION WITH RANDOMISED ELO RATING ###
# - From: A practical guide for inferring reliable dominance hierarchies and estimating their uncertainty - Sanchez-Tojar et al. 2017


curee.fls <- list.files(path = "./Data/Rocamadour_FeedingData", pattern = "curee_", all.files = FALSE, full.names = TRUE)[c(1,3:10,2)] # ordered files names
cureeDates <- c("2020-12-14", "2020-12-21", "2020-12-28", "2021-01-06", "2021-01-13",
                "2021-11-15", "2021-11-22", "2021-11-29", "2021-12-06", "2021-12-13")

curee_data_elo<- as.data.frame(do.call(rbind,
                                       lapply(1:length(curee.fls), function(f){
                                         
                                         #load data
                                         wantCureeInteractions <- read.csv(curee.fls[f])
                                         
                                         # Reshape data
                                         wantCureeInteractions <- wantCureeInteractions[-c(1:8,10),-c(5:7,10)] ## Remove useless rows
                                         colnames(wantCureeInteractions) <- wantCureeInteractions[c(1),] ## Take the first row as the columns names
                                         wantCureeInteractions <- wantCureeInteractions[-c(1),] # Remove first row as it became the colnames
                                         wantCureeInteractions$date <- cureeDates[f] # associate a date
                                         
                                         wantCureeInteractions
                                       })))

# Clean for missing names, when winner = looser and add a year column
curee_data_elo <- curee_data_elo %>% 
  dplyr::filter(.$`Total number of occurences` >= 1) %>% 
  dplyr::filter(.$Subject != .$Modifiers) %>% 
  dplyr::mutate(year = str_sub(.$date, end = 4))


## transform all loose interaction in win interaction and swap the individuals
all_int_elo <- data.frame(date = character(0), winner = character(0), looser = character(0), draw = character(0))

for (i in 1:nrow(curee_data_elo)){
  
  r1 <- curee_data_elo[i,]
  
  if (r1$Behavior == "win vs"){
    
    date <-  rep(r1$date,as.integer(r1$`Total number of occurences`))
    winner <- rep(r1$Subject,as.integer(r1$`Total number of occurences`))
    looser <- rep(r1$Modifiers,as.integer(r1$`Total number of occurences`))
    draw <-  rep("TRUE",as.integer(r1$`Total number of occurences`))
    
    x <- data.frame(date, winner, looser, draw)
    
    all_int_elo <- rbind(all_int_elo,x)
    
  }else{
    
    date <-  rep(r1$date,as.integer(r1$`Total number of occurences`))
    winner <- rep(r1$Modifiers,as.integer(r1$`Total number of occurences`))
    looser <- rep(r1$Subject,as.integer(r1$`Total number of occurences`))
    draw <-  rep("TRUE",as.integer(r1$`Total number of occurences`))
    
    x <- data.frame(date, winner, looser, draw)
    
    all_int_elo <- rbind(all_int_elo,x)
  }
}

# Create a dataset / year (different group)
all_int_elo2021 <- all_int_elo %>% 
  dplyr::filter(.$date %in% cureeDates[1:5])

all_int_elo2022 <- all_int_elo %>% 
  dplyr::filter(.$date %in% cureeDates[6:10]) %>% 
  dplyr::filter(.$winner != "Gregoire") # Remove one row containing this individual (error)




# List of elo ranks for both years
EloRanks.ls <- vector("list")



## 2021 Data
df <- all_int_elo2021
RandomEloRes<-  elo_scores(winners = df$winner, losers = df$looser, 
                                 identities = unique(df$winner), K = 200, init.score = 1000, 
                                 randomise = TRUE, n.rands = 1000, return.as.ranks = TRUE)

RandomEloRes

# Extract ranks from randomised time
rank <- rowMeans(RandomEloRes)
rank <- rank[order(rank)]
ids <- names(rank)
ranks <- 1:length(rank)

EloRanks.ls[["2021"]] <- as.data.frame(rank)
EloRanks.ls[["2021"]]$ranks <- 1:6 

# Shape of the hierarchy
shape1 <- plot_hierarchy_shape(fitted=FALSE, ids,ranks,
                              df$winner,
                              df$looser)

# Estimation of the uncertainty by repeatability
rept <- estimate_uncertainty_by_repeatability(df$winner, df$looser, 
                                              identities=ids,
                                              init.score=1000,
                                              n.rands = 1000)
round(rept,3) ## [1] 0.822


## 2021 Data
df <- all_int_elo2022
RandomEloRes<-  elo_scores(winners = df$winner, losers = df$looser, 
                           identities = unique(df$winner), K = 200, init.score = 1000, 
                           randomise = TRUE, n.rands = 1000, return.as.ranks = TRUE)

RandomEloRes

# Extract ranks from randomised time
rank <- rowMeans(RandomEloRes)
rank <- rank[order(rank)]
ids <- names(rank)
ranks <- 1:length(rank)

EloRanks.ls[["2022"]] <- as.data.frame(rank)
EloRanks.ls[["2022"]]$ranks <- 1:6 

# Shape of the hierarchy
shape2 <- plot_hierarchy_shape(fitted=FALSE, ids,ranks,
                              df$winner,
                              df$looser)

# Estimation of the uncertainty by repeatability
rept <- estimate_uncertainty_by_repeatability(df$winner, df$looser, 
                                              identities=ids,
                                              init.score=1000,
                                              n.rands = 1000)
round(rept,3) ## [1] 0.825


# Save hierarchy files
#save(EloRanks.ls, file="./Output/Files/EloRanks_ls.RData")



