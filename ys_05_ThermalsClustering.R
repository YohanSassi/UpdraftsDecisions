#________________________________________________________________________________________________________________________________________________
###
###       
###       Project: Explore the balance between personal and social information in thermal updraft selection
###
###        Data come from captive vulture from the Rocher des Aigles, Rocamadour, France experiencing free-flights
###        in a 120m deep canyon.
### 
###       !!! Script used to clusterise thermals / flight !!!
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
library(dplyr)
library(plot3D)
library(cowplot)
library(raster)
library(plotly)
library(rgl) 
library(scatterplot3d)
library(dbscan)

#______________________________________________________________#
###
##### Extract thermals for each individuals / flights ####
###
#______________________________________________________________#

# Parameters
inBetweenThermalThreshold_inSec = 5 # with a mean ground speed in flight = 40km/h, 10s = 100m in straight line
timeInThermal_inSec = 30


# Files list of segmented flights
FlightAnnot_fls <- list.files(path = "./Output/Files/FlightAnnot", full.names = TRUE, recursive = FALSE)

# -- LOOP ON FILES
for (f in 1:length(FlightAnnot_fls)) {
  
  # Load file
  load(FlightAnnot_fls[f])

  # Split movestack into move object
  flight_data.annot.splt <- split(flight_data.annot.mvstk)
  
  # -- LOOP ON FLIGHT
  for (v in 1:length(unique(flight_data.annot.mvstk@data$Behaviour))) {
    
    # Select thermals for each individuals in a specific flight #moveStack(
    allInd_thermal.mvstk <- 
      lapply(names(flight_data.annot.splt), function(x){
        
        # print(unique(x$Behaviour))
        #subdata <- split(flight_data.annot.mvstk)[x]
        subInd <- flight_data.annot.splt[[x]]
        #print(unique(subInd$Behaviour))
        
        if (paste0("Flight",v) %in% unique(subInd$Behaviour)) {
          subInd <- subInd[subInd$Behaviour == paste0("Flight",v) & subInd$BinarySG == 1]
          
          # Estimate time difference between consecutive point (to see if breaks in thermals)
          # add a ThermalID (+1 at the end to start with thermalID at 1)
          subInd$ThermalID <- c(1, cumsum(abs(difftime(subInd$Timestamp[1:length(subInd$Timestamp)-1],
                                                       subInd$Timestamp[2:length(subInd$Timestamp)],
                                                       units="sec")) > inBetweenThermalThreshold_inSec)+1)
          
          # Look for thermals that last for at least timeInThermal_inSec
          thermalToKeep <- subInd@data %>%
            group_by(ThermalID) %>%
            dplyr::summarise(duration = abs(difftime(first(Timestamp),last(Timestamp)))) %>%
            filter(duration >= timeInThermal_inSec)
          
          # Filter out the thermal that are too short (duration < timeInThermal_inSec)
          if (length(thermalToKeep$ThermalID) > 0) {
            subInd <- subInd[subInd$ThermalID %in% thermalToKeep$ThermalID]
          } else {
            subInd <- NULL
          }
          subInd
        }
      }) 
    if (length(allInd_thermal.mvstk %>% 
        discard(is.null)) == 0) {
      next
    } 
    
    allInd_thermal.mvstk <- moveStack(allInd_thermal.mvstk %>% 
                                        discard(is.null) # Remove object that are NULL (2 possible reasons: 
                                      # 1) ind didn't fly in this flight; 
                                      # 2) any of the detected thermal are long enough
    , forceTz= "UTC")
    
    # Save for each day
    date <- unique(allInd_thermal.mvstk@idData$Dates)
    save(allInd_thermal.mvstk , file = paste0("./Output/Files/ThermPerFlights/allInd_thermal_",date,"_Flight",v,".RData"))
    
  }
}



# VISUAL CHECK -- for 3D exploration
load(FlightAnnot_fls[21])
flight_data.annot.mvstk@idData$Indname
unique(flight_data.annot.mvstk$Behaviour)

ind = 3
flightnb <- 3

# Full annotated flight
want1 <- flight_data.annot.mvstk[[ind]][flight_data.annot.mvstk[[ind]]$Behaviour == unique(flight_data.annot.mvstk[[ind]]$Behaviour)[flightnb]]
unique(want1@idData$Indname)


date <- as.character(unique(want1@idData$Dates))

# Keep thermals
thermals_fls <- list.files(path = "./Output/Files/ThermPerFlights", full.names = TRUE, recursive = FALSE)
load(thermals_fls[which(grepl(paste(date,"_Flight",flightnb, sep=""), thermals_fls, fixed = TRUE))])


want2 <- allInd_thermal.mvstk[[want1@idData$Indname]]


# FULL FLIGHT for specific individual
plot_ly(want1@data,
        x = ~Longitude, 
        y = ~Latitude, 
        z = ~Height.above.msl, 
        type = 'scatter3d', 
        mode = 'lines', 
        color = ~BinarySG,
        colors = c("darkgreen", "orange")
        ) %>% 
  layout(title = paste0(want1@idData$Indname, " - ",unique(want1@idData$Dates), " - ",unique(want1@data$Behaviour)))


# ONLY THERMALS for specific individual
#test <- correctedThermalCoords.df[which(correctedThermalCoords.df$Indname == "Henry"),]
plot_ly(want2@data,
        x = ~Longitude, 
        y = ~Latitude, 
        z = ~Height.above.msl, 
        type = 'scatter3d', 
        mode = 'lines', 
        color = ~as.factor(ThermalID)
) %>% 
  layout(title = paste0(want2@idData$Indname, " - ",unique(want2@idData$Dates), " - ",unique(want2@data$Behaviour)),
         showlegend = FALSE)


# ONLY THERMALS BUT FOR ALL INDIVIDUALS IN THIS FLIGHT
test <- split(allInd_thermal.mvstk)

for (i in 1:length(unique(test))){
  ind1 <- test[[i]]@data %>% 
    mutate(Indname = rep(unique(test[[i]]@idData$Indname)),
           Dates = rep(unique(test[[i]]@idData$Dates)))
  
  if (i == 1) {
    allInd_thermal.df <-ind1
  } else {
    allInd_thermal.df <- rbind(allInd_thermal.df,ind1)
  }
}

# Plot
plot_ly(allInd_thermal.df,
        x = ~Longitude, 
        y = ~Latitude, 
        z = ~Height.above.msl, 
        type = 'scatter3d', 
        mode = 'lines', 
        color = ~Indname,
        colors = c("red", "green", "blue", "gold")
) %>% 
  layout(title = paste0( " All individual thermals - ",unique(allInd_thermal.df$Dates), " - ",unique(allInd_thermal.df$Behaviour)))

  


## -- Plots of variables through time
head(want2@data)
plot(want1@data$Timestamp, want1@data$vertSpeed_smooth, type = "l")
plot(want1@data$Timestamp, want1@data$Height.above.msl, type = "l")

want2@idData$Indname







#______________________________________________________________________________________#
###
##### Information on thermals for each individuals / flights & Spatial clustering ####
###
#______________________________________________________________________________________#


thermals_fls <- list.files(path = "./Output/Files/ThermPerFlights", full.names = TRUE, recursive = FALSE)

thermalInfo.df <- lapply(thermals_fls, function(x){
  
  # load data
  load(x)
  #load(thermals_fls[32])
  
  # Split to work on each individual
  allInd_thermal.splt <- split(allInd_thermal.mvstk)
  
  thermalInfo <- lapply(names(allInd_thermal.splt), function(y){
    
    # Select each individuals 
    subInd <- allInd_thermal.splt[[y]]
    
    # Create a new data frame with name of individuals and with a name_thermalnumber variable
    wantData <- as.data.frame(do.call(rbind,
                                      lapply(unique(subInd@data$ThermalID), function(z){
                                        # Select thermals
                                        wantThermals <- subInd@data[which(subInd@data$ThermalID == z),]
                                        
                                        # Add name again in dataset
                                        wantThermals$Indname <- rep(subInd@idData$Indname, nrow(wantThermals))
                                        
                                        # Add variable combining Thermalnb and identity
                                        wantThermals$IndTherm = paste0(wantThermals$Indname,"_",wantThermals$ThermalID)
                                        wantThermals
                                      })
    ))
    
    # Collect information on thermals 
    toAdd <- subInd@data %>% 
      group_by(ThermalID) %>% 
      dplyr::summarise(Dates = unique(subInd@idData$Dates),
                       Flight = unique(Behaviour),
                       IndName =  unique(subInd@idData$Indname),
                       StartTime = first(Timestamp),
                       ThermalID = unique(ThermalID),
                       Duration = abs(difftime(first(Timestamp),last(Timestamp), units = "secs")),
                       IndTherm = paste0(IndName,"_",ThermalID))
    
    return(list(wantData,toAdd))
  })
  
  # Convert back into data frame with all individuals from this flight
  wantData.flight <- as.data.frame(do.call(rbind,lapply(thermalInfo, `[[`, 1)))
  thermalInfo.flight <- as.data.frame(do.call(rbind,lapply(thermalInfo, `[[`, 2)))
  
  ## -- 3D spatial clustering of thermals using DBSCAN -- ##
  # Spatialize coordinates to change CRS to UTM
  wantData.coords <- SpatialPoints(cbind(wantData.flight$Longitude,wantData.flight$Latitude), 
                                   proj4string = CRS("+proj=longlat +ellps=WGS84"))
  wantData.coords <- spTransform(wantData.coords, 
                                 CRS(paste("+proj=utm","+zone=31","+ellps=WGS84", "+datum=WGS84", "+units=m", "+towgs84:0,0,0", sep=" ")))
  
  # Reformat data
  wantData.flight <- as.data.frame(cbind(wantData.flight,wantData.coords)) %>% 
    dplyr::rename(LatitudeUTM = coords.x1, LongitudeUTM = coords.x2) %>% 
    mutate(LatitudeUTM = as.numeric(LatitudeUTM),
           LongitudeUTM = as.numeric(LongitudeUTM),
           Height.above.msl = as.numeric(Height.above.msl))
  
  # Process the clustering algorithm
  #kNNdistplot(wantData.flight[,c("LatitudeUTM", "LongitudeUTM", "Height.above.msl")], k = 4) #esp should be derived from Y axis when the curve become vertical
  wantData.flight_db <- dbscan::dbscan(wantData.flight[,c("LatitudeUTM", "LongitudeUTM", "Height.above.msl")],
                                       eps = 40, 
                                       minPts = 5) # Spatial clustering process
  
  # Visual check 2D & 3D
    # factoextra::fviz_cluster(wantData.flight_db, wantData.flight[,c("LongitudeUTM","LatitudeUTM")], stand = FALSE, ellipse = TRUE, geom = "point")
    # plot3d(wantData.flight$LongitudeUTM, wantData.flight$LatitudeUTM, wantData.flight$Height.above.msl,
    #        col = as.factor(wantData.flight_db$cluster+1), pch = 20, cex = 2, type = "p", ize = 5)

  # Associate cluster to thermals
  wantData.flight$ClusterNb <- wantData.flight_db$cluster
  
  # Check that all thermals is associated to only one cluster & correct if not
  # checkTherm <- wantData.flight %>%
  #   group_by(IndTherm) %>%
  #   dplyr::summarise(IndTherm = unique(IndTherm),
  #                    ClustAssociated = length(unique(ClusterNb))) %>%
  #   filter(ClustAssociated > 1)
  
  allInd_thermalClust.df <- wantData.flight %>%
    add_count(IndTherm, ClusterNb) %>% 
    group_by(IndTherm) %>%
    dplyr::mutate(ClusterNb = ClusterNb[which.max(n)]) %>% 
    ungroup() %>% 
    dplyr::select(-n)
  
  # plot3d(allInd_thermalClust.df$LongitudeUTM, allInd_thermalClust.df$LatitudeUTM, allInd_thermalClust.df$Height.above.msl,
  #        col = rainbow(length(unique(allInd_thermalClust.df$ClusterNb)))[as.factor(allInd_thermalClust.df$ClusterNb)], pch = 20
  #        , cex = 2, type = "p", ize = 5) # Check plot

  
  # Extract name and save file with annotated clusters
  hyfenindices <- str_locate_all(pattern ='_', x)
  EndName <- substr(x, 
                    as.integer(hyfenindices[[1]][2,1]), 
                    str_length(x)) # Extract date and Flight of the considered file to save
  
  # save(allInd_thermalClust.df, 
  #      file = paste0("./Output/Files/ClusteredThermals/allInd_thermalClust",EndName)) # Save
   
  return(list(thermalInfo.flight)) #allInd_thermalClust.df,
})

# Convert back into data frame
#allFlightData.df <- as.data.frame(do.call(rbind,lapply(thermalInfo.df, `[[`, 1)))
thermalInfo.df <- as.data.frame(do.call(rbind,lapply(thermalInfo.df, `[[`, 1)))

# save(thermalInfo.df,
#      file = "./Output/Files/thermalInfo.RData")







# ---- Example of spatial clustering of thermals - full flight in 3D - Fig1  -----
par(mar = c(1,1,1,1))
scatter3D(x = allInd_thermalClust.df$Longitude, 
          y = allInd_thermalClust.df$Latitude, 
          z = allInd_thermalClust.df$Height.above.msl,
          phi = 20, # Angle from which to see the graph
          theta = -60,
          colkey = FALSE,
          colvar = as.integer(allInd_thermalClust.df$ClusterNb),
          col = c("#66C2A5","darkolivegreen2","cornflowerblue"), #brewer.pal(3, "Set2"),
          bty = "b",
          type = "p", # l for lines only
          #ticktype = "detailed", # Detailed axis
          pch = 20,
          cex = 0.5,
          lwd = 1,
          xlim = c(min(want1$Longitude),max(want1$Longitude)),
          ylim = c(min(want1$Latitude),max(want1$Latitude)),
          zlim = c(min(want1$Height.above.msl),max(want1$Height.above.msl)),
          xlab = "Longitude",
          ylab = "Latitude",
          zlab = 'Altitude',
          main="Thermal clustering")

