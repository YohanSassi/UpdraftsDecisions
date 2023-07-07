#________________________________________________________________________________________________________________________________________________
###
###       
###       Project: Explore the balance between personal and social information in thermal updraft selection
###
###        Data come from captive vulture from the Rocher des Aigles, Rocamadour, France experiencing free-flights
###        in a 120m deep canyon.
### 
###       !!! Script used to segment flight data - Soaring / gliding !!!
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

library(stringr)
library(lubridate)
library(move)
library(tidyr)
library(dplyr)
library(plot3D)
library(plyr)
library(doParallel)
library(plotly)



#______________________________________________________________#
###
##### Differentiate soaring from gliding in vulture flights ####
###
#______________________________________________________________#


####
# Credit: Martina Scacco, adapted by Yohan Sassi
# The script imports a csv containing movement data with at least some GPS bursts collected at 1 Hz frequency (the segmentation is only applied to 1 Hz data)
# the dataframe is then converted into a move object to calculate movement metrics needed for the behavioural segmentation 
# into soaring/gliding flight, circling/linear flight and their combination
#
####

#____________________________________________________________________________
### Import data, transform in move object and calculate track variables ####
#____________________________________________________________________________

#Parameters
# If you have multiple cores you can run this in parallel
detectCores()
doParallel::registerDoParallel(detectCores()-2)

## Parameters
minResol <- 1 # 1 to max 2 sec timelag
minBurstDuration <- 30 # we want bursts of at least 30 secs

swV <- 2 #smoothing window of 5 seconds (< min burst duration, 2 before 2 after each loc) for vertical speed for later classification
swT <- 14 #smoothing window of 29 seconds for thermalling behaviour (according to Rolf/Bart's paper)



# ---
# load files names
dailyFilesNames <- list.files(path = "./Output/Files/Fullday_Flights", full.names = TRUE, recursive = FALSE)

underscoreindices <- str_locate_all(pattern ='_', dailyFilesNames[1])
endindices <- str_locate_all(pattern ='.RData', dailyFilesNames[1])
dailyFilesNames <- dailyFilesNames[order(dmy(substr(dailyFilesNames, 
                                                    as.integer(underscoreindices[[1]][3,1])+1, 
                                                    as.integer(endindices[[1]][1,1])-1)))] # Ordered files names


# To transform your dataframe in a move object we need the columns with the geographic coordinates, a column with the time and with the names of the individuals
# Here we assume that the columns are called as follows, change accordingly:
# lat, long (in latlong wgs84 projection)
# time (in the format year-month-day hour:minute:second)
# individual

for (d in 1:length(dailyFilesNames)) { # LOOP ON EACH FULLDAY FLIGHT FILES
  
  # load the considered flight data file
  load(dailyFilesNames[d]) # d # extract all flights -> flight_data.df
  
  # Arrange time to match with move requirement
  df <- flight_data.df %>% 
    mutate(Dates = as.character(as.Date(Date, "%d/%m/%Y"),"%Y-%m-%d")) %>% 
    mutate(Times = as.character(hms::as_hms(Time))) %>% 
    mutate(Timestamp = as.POSIXct(paste(Dates, Times), format="%Y-%m-%d %H:%M:%S")) %>% 
    dplyr::select(c(Timestamp,Latitude,Longitude,Height.above.msl,Ground.speed.km.h,Dates,Behaviour,Indname))
  
  # Create a movestack
  mvstk <- move(x=df$Longitude, y=df$Latitude, proj=CRS("+proj=longlat +ellps=WGS84"),
                time=as.POSIXct(df$Timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
                animal=df$Indname,
                data=df)
  
  # Split the movestack in a list of move objects
  mv_ls <- split(mvstk) 
  
  # Keep only move objects (individuals) with more than 30 locations (minimum needed for thermal segmentation)
  mv_ls <- mv_ls[sapply(mv_ls, n.locs)>=30]
  
  # Calculate variables related to the trajectory -> for segmentation after
  allmv.ls <- lapply(names(mv_ls), function(animalName){
    mv <- mv_ls[[animalName]]
    mv$timelag.sec <- c(NA, timeLag(mv, units="secs"))
    mv$altitude.diff <- c(NA, (mv$Height.above.msl[-1] - mv$Height.above.msl[-nrow(mv)]))
    mv$vert.speed <- mv$altitude.diff/mv$timelag.sec
    mv$turn.angle <- c(NA, turnAngleGc(mv), NA)
    mv$step.length <- c(NA, distance(mv)) # should this give n error try the following line instead:
    #mv$step.length <- c(NA, raster::pointDistance(mv[-sum(n.locs(mv)),], mv[-1,], longlat=isLonLat(mv)))
    #mv$gr.speed <- c(NA, speed(mv)) # should this give n error try the following line instead:
    #mv$gr.speed <- c(NA, mv$step.length/mv$timelag.sec)
    mv
    #animalId <- mv@idData$Indname
    #date <- unique(mv@idData$Dates)
    #save(mv, file = paste0("./Data/Rocamadour_FlightData/ForSeg/",animalId,"_",date,"_mv.RData"))
  })
  
  # Save the list of move object for the day, that will be used in next step
  # date <- unique(mvstk@idData$Dates)
  # save(allmv.ls, file = paste0("./Data/Rocamadour_FlightData/ForSeg/allmv.ls_",date,".RData"))
  
  #_________________________________________________________
  ### Select bursts of continuous 1 sec resolution data ####
  #_________________________________________________________
  
  
  flight_data.annot.ls <- llply(1:length(allmv.ls), function(i){ #fls lapply
    print(i)
    #load(fls[f]) #object mv
    
    # Isolate each individual move object
    mv <- allmv.ls[[i]]
    
    animalId <- mv@idData$Indname
    
    # with cumsum R assigns the same value to all consecutive locations for which the condition is false (timelag <= 1 sec)
    mv$burstID <- c(0, cumsum(mv$timelag.sec[2:n.locs(mv)] > minResol))  #from row 2 (since the first is NA)
    
    # with table we can count all locations with the same ID (each one is a separate burst) and keep those high resolution bursts that have at least a certain number of locations (= minBurstDuration)
    burstDuration <- as.data.frame(table(mv$burstID))
    burstsToKeep <- burstDuration$Var1[which(burstDuration$Freq >= minBurstDuration)]
    
    # use those to subset the move obj and keep only the high resolution bursts
    HRmv <- mv[which(mv$burstID %in% burstsToKeep),]
    if(nrow(HRmv)>0){
      HRdf_bursts <- as.data.frame(HRmv)
      
      # Remove unnecessary columns
      HRdf_bursts <- HRdf_bursts[,-grep("mag|orientation|coords|timestamps|start.timestamp|optional|import|visible|algorithm|battery|decoding|accuracy|manually|activity|checksum|acceleration",
                                        colnames(HRdf_bursts))]
      
      # Split each individual dataframe by burst ID
      burst_ls_corr <- split(HRdf_bursts, HRdf_bursts$burstID)
      
      # Keep only bursts with at least 40 s (30 of smoothing window will be NA) 
      burst_ls_corr_sub <- burst_ls_corr[which(sapply(burst_ls_corr, nrow) >= 40)]
      
      # Compute smoothed turning angle separately for each burst
      HRdf <- do.call(rbind, lapply(burst_ls_corr_sub, function(b){
        b$vertSpeed_smooth <- NA
        b$turnAngle_smooth <- NA
        for(i in (swV+1):(nrow(b)-swV)){ 
          b$vertSpeed_smooth[i] <- mean(b$vert.speed[(i-swV):(i+swV)], na.rm=T)}
        for(i in (swT+1):(nrow(b)-swT)){
          b$turnAngle_smooth[i] <- max(abs(cumsum(b$turn.angle[(i-swT):(i+swT)])))} # !!!! Test with abs after cumsum !!!!
        return(b) # return df with smoothed variables
      }))
      
      # Classify soaring only based on vertical speed
      HRdf <- HRdf[complete.cases(HRdf$vertSpeed_smooth),]
      kmeanV <- kmeans(HRdf$vertSpeed_smooth, 2)   #Get the two clusters
      soarId <- which.max(aggregate(HRdf$vertSpeed_smooth~kmeanV$cluster, FUN=mean)[,2]) # which one is the soaring one?
      soarClust <- rep("glide", length(kmeanV$cluster))
      soarClust[which(kmeanV$cluster==soarId)] <- "soar"
      HRdf$soarClust <- factor(soarClust, levels=c("soar","glide"))  
      
      # Now classify thermalling only based on turning angle (cumulated to a 30 s time window in previous step)
      HRdf$thermalClust <- "other"
      HRdf$thermalClust[which(HRdf$gr.speed >= 2 & HRdf$turnAngle_smooth >= 300)] <- "circular"
      HRdf$thermalClust[which(HRdf$gr.speed >= 2 & HRdf$soarClust=="soar" & HRdf$thermalClust != "circular")] <- "linear"
      HRdf$thermalClust <- factor(HRdf$thermalClust, levels=c("circular","linear","other"))
      
      # Convert to move object
      flight_data.annot <- move(x=HRdf$Longitude, y=HRdf$Latitude, proj=CRS("+proj=longlat +ellps=WGS84"),
                                time=as.POSIXct(HRdf$Timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
                                animal=HRdf$Indname,
                                data=HRdf)
      
      # annotate in binary circular soaring / gliding to plot
      flight_data.annot@data <- flight_data.annot@data %>% 
        mutate(BinarySG = ifelse(as.character(soarClust) == "soar" & as.character(thermalClust) == "circular", 1, 0))
      
      # Save classified dataframe per individual
      #save(flight_data.annot, file = paste0(animalId,"_classifiedBursts_df.rdata"))
    }
    
    flight_data.annot
  }, .parallel=T)
  #}, .parallel=T) #for parallel
  
  # Combine in movestack and save
  flight_data.annot.mvstk <- moveStack(flight_data.annot.ls, forceTz= "UTC")
  date <- unique(flight_data.annot.mvstk@idData$Dates)
  save(flight_data.annot.mvstk, file = paste0("./Output/Files/FlightAnnot/flight_data.annot_",date,".RData"))
}




#_________________________________________________________
### 3D plot for visual check ####
#_________________________________________________________

# Files list of segmented flights
FlightAnnot_fls <- list.files(path = "./Output/Files/FlightAnnot", full.names = TRUE, recursive = FALSE)
load(FlightAnnot_fls[21])

# Choose specific individual
flight_data.annot.mvstk@idData$Indname
#flight_data.annot.mvstk@idData$Dates
ind = 3

# Choose specific flight
unique(flight_data.annot.mvstk[[ind]]$Behaviour)
flightnb <- 4

# Select specific individual & flight that day
want2 <- flight_data.annot.mvstk[[ind]][flight_data.annot.mvstk[[ind]]$Behaviour == unique(flight_data.annot.mvstk[[ind]]$Behaviour)[flightnb]]


# ---- Interactive plot of segmented flight in 3D -----
plot_ly(want3@data,
               x = ~Longitude, 
               y = ~Latitude, 
               z = ~Height.above.msl, 
               type = 'scatter3d', 
               mode = 'lines', 
               color = ~BinarySG,
               colors = c("darkgreen", "orange")) %>% 
  layout(title = paste0(want@idData$Indname, " - ",unique(want@idData$Dates), " - ",unique(want@data$Behaviour)),
         showlegend = FALSE)




# ----- Example of segmented flight in 3D - Fig1  -----

par(mar = c(1,1,1,1))
scatter3D(x = want2$Longitude, 
          y = want2$Latitude, 
          z = want2$Height.above.msl,
          phi = 20, # Angle from which to see the graph
          theta = -60,
          colkey = FALSE,
          colvar = as.integer(want2$BinarySG),
          col = c("#8DA0CB","darkorange3"),
          bty = "b",
          type = "l", # l for lines only
          #ticktype = "detailed", # Detailed axis
          pch = 20, # size of points
          #cex = c(2, 2, 2),
          lwd = 1,
          xlim = c(min(want1$Longitude),max(want1$Longitude)),
          ylim = c(min(want1$Latitude),max(want1$Latitude)),
          zlim = c(min(want1$Height.above.msl),max(want1$Height.above.msl)),
          xlab = "Longitude",
          ylab = "Latitude",
          zlab = 'Altitude',
          main="Segmented flight")

