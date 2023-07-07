#________________________________________________________________________________________________________________________________________________
###
###       
###       Project: Explore the balance between personal and social information in thermal updraft selection
###
###        Data come from captive vulture from the Rocher des Aigles, Rocamadour, France experiencing free-flights
###        in a 120m deep canyon.
### 
###       !!! Script used to associated variables of interest to each thermals !!!
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
library(dplyr)
library(lubridate)
library(purrr)
library(raster)
library(move)



source("./Functions/ys_recoverMissingFiness.R")
source("./Functions/df_to_spatial_df_utm.R") #To convert long/lat to UTM
source("./Functions/getDistance3D.R") #To estimate 3D distances between projected locations




#_________________________________________________________________________________________________#
###
##### Extract inter-thermal glide for each individuals / flights to estimate Finesse ####
###
#_________________________________________________________________________________________________#

# Parameters
timeInGlide_inSec = 5 # with a mean ground speed in flight = 40km/h, 10s = 100m in straight line
AltThreshold = 300

# Files list of segmented flights
FlightAnnot_fls <- list.files(path = "./Output/Files/FlightAnnot", full.names = TRUE, recursive = FALSE)

# Create a file to store all finesse values per flight / individuals
glideToKeep_Allvariables.df <- data.frame(Date = character(0), 
                                          Flight = integer(0), 
                                          IndName = character(0), 
                                          GlideID = integer(0),
                                          Duration = character(0),
                                          CumDist = numeric(0),
                                          LinearDist = numeric(0),
                                          AltLoss = numeric(0),
                                          Straightness = numeric(0),
                                          Finesse = numeric(0)
)


# -- LOOP ON FILES
for (f in 1:length(FlightAnnot_fls)) {
  
  # Load file
  load(FlightAnnot_fls[f])
  
  # Split movestack into move object
  flight_data.annot.splt <- split(flight_data.annot.mvstk)
  
  # -- LOOP ON FLIGHT
  for (v in 1:length(unique(flight_data.annot.mvstk@data$Behaviour))) {
    
    # Select glides for each individuals in a specific flight
    glideToKeep_variables.ls <- 
      lapply(names(flight_data.annot.splt), function(x){
        
        subInd <- flight_data.annot.splt[[x]]
        
        if (paste0("Flight",v) %in% unique(subInd$Behaviour)) {
          subInd <- subInd[subInd$Behaviour == paste0("Flight",v) & subInd$soarClust == "glide" & subInd$thermalClust == "other"]
          
          # Estimate time difference between consecutive point (to see if breaks in gliding phase)
          # add a glideID (+1 at the end to start with glideID at 1)
          subInd$GlideID <- c(1, cumsum(abs(difftime(subInd$Timestamp[1:length(subInd$Timestamp)-1],
                                                     subInd$Timestamp[2:length(subInd$Timestamp)],
                                                     units="sec")) > 1)+1)
          
          # Look for glide that last for at least timeInGlide_inSec and which all glide happen above the AltThreshold
          # Filter out the Glide that are too short (duration < timeInGlide_inSec)
          glideToKeep <- subInd@data %>%
            group_by(GlideID) %>%
            dplyr::summarise(Duration = abs(difftime(first(Timestamp),last(Timestamp))),
                             n = n(),
                             nSupThreshold = sum(Height.above.msl >= AltThreshold)) %>%
            filter(Duration >= timeInGlide_inSec & n == nSupThreshold) %>%
            dplyr::select(.,-c(n, nSupThreshold))
          
          # If there are glides with the previous characteristics, then extract the needed variables
          if (length(glideToKeep$GlideID) > 0) {
            
            glideToKeep_variables <- subInd@data %>%
              mutate(IndName = subInd@idData$Indname) %>% 
              group_by(GlideID) %>%
              dplyr::summarise(Date = unique(date(.$Timestamp[1])),
                               Flight = unique(.$Behaviour),
                               IndName = unique(.$IndName),
                               Duration = abs(difftime(first(Timestamp),last(Timestamp))),
                               n = n(),
                               nSupThreshold = sum(Height.above.msl >= AltThreshold)) %>%
              filter(Duration >= timeInGlide_inSec & n == nSupThreshold) %>%
              dplyr::select(.,-c(n, nSupThreshold)) %>% 
              bind_cols(do.call("rbind",lapply(1:nrow(.), function(g){
                
                # Extract the considered glide
                wantglide <- subInd@data[which(subInd@data$GlideID == .$GlideID[g]),]
                
                # Estimate cumulated distance
                CumDist <- sum(pointDistance(wantglide[1:nrow(wantglide)-1,c("Longitude","Latitude")],
                                             wantglide[2:nrow(wantglide),c("Longitude","Latitude")],
                                             lonlat = TRUE))
                
                # Estimate the linear distance between first and last fix of the glide
                LinearDist <- pointDistance(wantglide[1,c("Longitude","Latitude")],
                                            wantglide[nrow(wantglide),c("Longitude","Latitude")],
                                            lonlat = TRUE)
                
                # Estimate the altitude loss between first and last fix of the glide
                AltLoss <- wantglide[1,"Height.above.msl"] - wantglide[nrow(wantglide),"Height.above.msl"]
                
                # Estimate Straightness
                Straightness <- LinearDist / CumDist 
                
                # Estimate Finesse
                Finesse <- LinearDist / AltLoss

                return(c(CumDist, LinearDist, AltLoss, Straightness, Finesse))
              }))) %>% 
              dplyr::rename(CumDist = ...6,
                            LinearDist = ...7,
                            AltLoss = ...8,
                            Straightness = ...9,
                            Finesse = ...10) %>% 
              relocate(GlideID, .after = IndName)
            
            # Change Duration format otherwise make problem when aggregated in the dataframe
            glideToKeep_variables$Duration <- as.character(glideToKeep_variables$Duration)
            
          } else {
            glideToKeep_variables <- NULL
          }
          glideToKeep_variables
        }
      }) 
    
    if (length(glideToKeep_variables.ls %>% 
               discard(is.null)) == 0) {
      next
    } 
    
    glideToKeep_Flightvariables.df <- as.data.frame(do.call(rbind,glideToKeep_variables.ls %>% 
                                                              discard(is.null))) # Remove object that are NULL (2 possible reasons: 
                                                                                 # 1) ind didn't fly in this flight; 
                                                                                 # 2) any of the detected glide are long enough
    # Save for each Flight into a larger object
    glideToKeep_Allvariables.df <- rbind(glideToKeep_Allvariables.df, 
                                         glideToKeep_Flightvariables.df)
  }
}

# Save full dataframe as a separate file
# save(glideToKeep_Allvariables.df, file = paste0("./Output/Files/GlideToKeep_df.RData"))
#write.csv(glideToKeep_Allvariables.df, file = paste0("./Output/Files/GlideToKeep_df.csv"))


glideFiltered <- glideToKeep_Allvariables.df %>% 
  #mutate(Date_Flight = paste0(Date,"_",Flight)) %>% 
  filter(Straightness > 0.95) %>% # Consider only straight lines
  group_by(Date, Flight, IndName) %>%  #Date, Flight, 
  summarise(n = n(),
            Finesse = max(Finesse), # Needed for the function 
            Min_Finesse = min(Finesse),
            Mean_Finesse = mean(Finesse))


glideFiltered %>% 
  group_by(.$IndName) %>% 
  summarise(n())


# Variability of the Finenesse
ggplot(glideFiltered, aes(glideFiltered$Date_Flight, glideFiltered$`max(Finesse)`, group = glideFiltered$IndName)) + 
  geom_point(aes(col = glideFiltered$IndName)) +
  geom_line(aes(col = glideFiltered$IndName)) +
  theme(axis.text.x = element_text(angle = 45))

# See it into boxplot
ggplot(glideFiltered, aes(glideFiltered$IndName, glideFiltered$`max(Finesse)`)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45))


  



#_________________________________________________________________________________#
###
##### Associated variables of interest to thermals / flights & individuals ####
###
#_________________________________________________________________________________#

# Load all thermals info, Flight info, birds data, list of social bonds matrices and hierarchy levels
load("./Output/Files/thermalInfo.RData") #(open thermalInfo.df)
load("./Data/Rocamadour_FlightData/AllFlightData.info.RData") #(open list allFlightData)
allFlightData[["27:05:2021"]]$Day_flights[3,9:10] <- NA # Problem of GPS - removed otherwise induce problems in further loops


vultInfo.df <- as.data.frame(cbind(c("Kazimir","Henry","Gregoire","Leon","Kirikou","Bulma","Mathilda","Hercule"),
                     c(2014,1997,2011,2016,2014,2014,2018,2012))) %>% 
                dplyr::rename(., Name = V1, yearOfBirth = V2) %>% 
                mutate(Age = 2022 - as.integer(yearOfBirth))

# SBSMatrix list - both years -
SBSMatrix.ls <- vector("list")
load("./Output/Files/SBSMatrix_2021.RData") # 2021 social bond matrices & add capital letters to individual names (to match with data)
SBSMatrix.ls[["2021"]] <- SBSMatrix_estimatedreord
colnames(SBSMatrix.ls[["2021"]]) <- str_to_title(colnames(SBSMatrix.ls[["2021"]]))
rownames(SBSMatrix.ls[["2021"]]) <- str_to_title(rownames(SBSMatrix.ls[["2021"]]))

load("./Output/Files/SBSMatrix_2022.RData") # 2022 social bond matrices
SBSMatrix.ls[["2022"]] <- SBSMatrix_estimatedreord
colnames(SBSMatrix.ls[["2022"]]) <- str_to_title(colnames(SBSMatrix.ls[["2022"]]))
rownames(SBSMatrix.ls[["2022"]]) <- str_to_title(rownames(SBSMatrix.ls[["2022"]]))


# EloRanks - Hierarchy
load("./Output/Files/EloRanks_ls.RData") #(Open list EloRanks.ls)

# Individuals Finesse in flight
load("./Output/Files/GlideToKeep_df.RData") #(Open glideToKeept_Allvariables.df)

glideFiltered <- glideToKeep_Allvariables.df %>%
  #mutate(Date_Flight = paste0(Date,"_",Flight)) %>%
  filter(Straightness > 0.95) %>% # Consider only straight lines
  group_by(Date, Flight, IndName) %>%  #Date, Flight,
  summarise(n = n(),
            Finesse = max(Finesse), # Needed for the function
            Min_Finesse = min(Finesse),
            Mean_Finesse = mean(Finesse))

# Extract the list of names of GPS files per flight / day
thermalClust_fls <- list.files(path = "./Output/Files/ClusteredThermals", full.names = TRUE, recursive = FALSE)


# Extract variables of interest
AllThermalsInfos <- as.data.frame(do.call(rbind,
                                          lapply(thermalClust_fls, function(f){
                                            
                                            # load GPS data file of the flight
                                            load(f)
                                            #load(thermalClust_fls[25])
                                            
                                            # Reorganize GPS data.frame to free memory
                                            allInd_thermalClust.df <- allInd_thermalClust.df %>% 
                                              dplyr::select(-c(gr.speed,altitude.diff,vert.speed,turn.angle,burstID, soarClust, thermalClust, BinarySG))
                                            
                                            # Extract info on thermal that day/flight 
                                            flightDate <- date(allInd_thermalClust.df$Timestamp[1])
                                            flightNb <- unique(allInd_thermalClust.df$Behaviour)
                                            
                                            wantThermalInfo <- thermalInfo.df[which(thermalInfo.df$Dates == flightDate & thermalInfo.df$Flight == flightNb),]
                                            
                                            # Add cluster id, 
                                            # reorder to have order of appearance of thermals, 
                                            # assign if individual discovered or not the considered thermal (Discovery)
                                            # assign to each thermal/individual the delay between entrance and discovery of the thermal (Dt_Discovery)
                                            orderedThermals <- wantThermalInfo %>%
                                              mutate(ClusterNb = unlist(lapply(1:nrow(wantThermalInfo), function(a){
                                                wantCluster <- unique(allInd_thermalClust.df[which(allInd_thermalClust.df$IndTherm == wantThermalInfo[a,"IndTherm"]), "ClusterNb"])
                                              }))) %>% 
                                              arrange(StartTime) %>%
                                              group_by(ClusterNb) %>%
                                              dplyr::mutate(Discovery = ifelse(row_number() == 1,1,0),
                                                            Dt_Discovery = ifelse(row_number() == 1,0,difftime(StartTime[1], StartTime))) %>%
                                              ungroup() %>% 
                                              bind_cols(do.call("rbind",lapply(1:nrow(.), function(b){
                                                # add altitude at the entrance of the considered thermal (AltEntrance)
                                                AltEntrance <- allInd_thermalClust.df[which(allInd_thermalClust.df$IndTherm == as.character(.$IndTherm[b]) & 
                                                                                              allInd_thermalClust.df$Timestamp == .$StartTime[b]), "Height.above.msl"][[1]]
                                                # add age of individual (Age)
                                                Age <- vultInfo.df[which(vultInfo.df$Name == .$IndName[b]),"Age"]
                                                
                                                # Extract date and change format (-> needed to extract information from allFlightData list as slots are dates)
                                                wantDate <- format(as_date(.$Dates[b]), "%d:%m:%Y")
                                                
                                                # add if this individual was release with first or second group (ReleaseOrder)
                                                ReleaseOrder <- allFlightData[[wantDate]]$Day_flights[which(allFlightData[[wantDate]]$Day_flights$Individual == .$IndName[b] & 
                                                                                                              allFlightData[[wantDate]]$Day_flights$Flight_nb == str_sub(flightNb, -1)),"Group"]
                                                # add temperature (Temp)
                                                Temp <- unique(allFlightData[[wantDate]]$Day_flights$Temp)
                                                
                                                # add nebulosity - cloud cover on scale from 0 to 8 (Nebulosity)
                                                Nebulosity <- unique(allFlightData[[wantDate]]$Day_flights$Neb)
                                                
                                                # add windspeed - categorial variable (see later if should be replaced by numeric)
                                                WindSpeed <- ifelse(length(unique(allFlightData[[wantDate]]$Day_flights$Windspeeed)) == 1, 
                                                                    unique(allFlightData[[wantDate]]$Day_flights$Windspeeed),
                                                                    names(sort(table(allFlightData[[wantDate]]$Day_flights$Windspeeed), decreasing = TRUE)[1]))
                                                
                                                # add time first individual take off (start of the flight)
                                                TimeFirstTakeOff <- as.character(hms::as_hms(min(hms::as_hms(str_subset(allFlightData[[wantDate]]$Day_flights[which(
                                                  allFlightData[[wantDate]]$Day_flights$Flight_nb == str_sub(flightNb, -1)),"TakeOff_UTC"], ":")))))
                                                
                                                return(c(Age, Temp, Nebulosity, WindSpeed, ReleaseOrder, TimeFirstTakeOff, AltEntrance))
                                              }))) %>% 
                                              dplyr::rename(Age = ...11,
                                                            Temp = ...12,
                                                            Nebulosity = ...13,
                                                            WindSpeed = ...14,
                                                            ReleaseOrder = ...15,
                                                            TimeFirstTakeOff = ...16,
                                                            AltEntrance = ...17)
                                            

                                            # Variables needing a loop to be added (need to come back to previous thermals etc...)
                                            Dt_LastUser <- rep(NA, nrow(orderedThermals)) # Time delay between its entrance in the thermal and when the last user left
                                            VertSpeedmax <- rep(NA, nrow(orderedThermals)) # Maximal vertical speed reached by individuals before it arrived in the thermal
                                            SBmax <- rep(NA, nrow(orderedThermals)) # Is my preferred affiliates present in the thermal I join?
                                            NbIndUsingTherm <- rep(0, nrow(orderedThermals)) # Number of individual using the thermal before he arrived
                                            NbIndUsedTherm <- rep(0, nrow(orderedThermals)) # Number of individual that used the thermal before focal individual arrived
                                            IndAlreadyUsedThermCluster <- rep(NA, nrow(orderedThermals)) # Either the focal individual had already used this cluster (1/0)
                                            DistThermBefore <- rep(NA, nrow(orderedThermals)) # Distance with the previous thermal this individual used
                                            AltExitBefore <- rep(NA, nrow(orderedThermals)) # Altitude at which the individual left the previous therm
                                            EloInd <- rep(NA, nrow(orderedThermals)) # Elo rank as a proxy of individual position in the hierarchy
                                            Dt_Elo <- rep(NA, nrow(orderedThermals)) # Max difference of Elo rank between in focal ind and those in considered thermal
                                            SBmeanPast <- rep(NA, nrow(orderedThermals)) # Maximal bond strength share with individuals in the thermal before he arrived
                                            FullAbove <- rep(NA, nrow(orderedThermals)) # Is the selected thermal fully above the previous one?
                                            OlderPresent <- rep(0, nrow(orderedThermals)) # Older individual present in the thermal?
                                            
                                            for (i in 1:nrow(orderedThermals)) {
                                              
                                              # Select a specific thermal
                                              wantThermal <- orderedThermals[i,]
                                              
                                              # Other thermals available for this individual?
                                              indThermals <- orderedThermals[which(orderedThermals$IndName == wantThermal$IndName),]
                                              
                                              # If yes, what is the index of wantThermal within indThermals 
                                              indexwantThermal <- which(indThermals$IndTherm == wantThermal$IndTherm)
                                              
                                              # If it's the first row of thermals used by ind it must come from RdA & also true if the individual took only one thermal in the flight 
                                              # Estimate the distance between the center of the Esplanade at RdA & the beginning of the first thermal if indexwantThermal == 1, 
                                              # otherwise with previous thermal
                                              wantThermalFirstloc <- allInd_thermalClust.df[which(allInd_thermalClust.df$IndTherm == wantThermal$IndTherm),] %>% head(1)
                                              if(indexwantThermal == 1){
                                                wantPreviousThermalLastloc <- data.frame(Longitude=1.611650, Latitude=44.801389, Height.above.msl = 264,
                                                                                         Timestamp = as.POSIXct(paste(wantThermal$Dates , as.character(hms::as_hms(hms::as_hms(str_subset(allFlightData[[format(as_date(wantThermal$Dates), "%d:%m:%Y")]]$Day_flights[which(
                                                                                           allFlightData[[format(as_date(wantThermal$Dates), "%d:%m:%Y")]]$Day_flights$Flight_nb == str_sub(wantThermal$Flight, -1) & 
                                                                                             allFlightData[[format(as_date(wantThermal$Dates), "%d:%m:%Y")]]$Day_flights$Individual == wantThermal$IndName),"Release_UTC"], ":")))))))
                                              } else {
                                                wantPreviousThermalLastloc <- as.data.frame(allInd_thermalClust.df[which(allInd_thermalClust.df$IndTherm == indThermals[indexwantThermal-1,]$IndTherm),c("Longitude","Latitude", "Height.above.msl","Timestamp")] %>% 
                                                                                              tail(1))
                                              }
                                              
                          
                                           
                                              
                                              
                                              
                                              # DistThermBefore[i] <- pointDistance(cbind(wantPreviousThermalLastloc$Longitude, wantPreviousThermalLastloc$Latitude),
                                              #                                     cbind(wantThermalFirstloc$Longitude,wantThermalFirstloc$Latitude),
                                              #                                     lonlat = TRUE)
                                              # 
                                              # -- FOR 3D DISTANCE
                                              # Transform previous thermal exit points into UTM
                                              p1 <- data.frame(x = wantPreviousThermalLastloc$Longitude, y = wantPreviousThermalLastloc$Latitude, z = wantPreviousThermalLastloc$Height.above.msl) # Last point in the previous thermal
                                              df_to_spatial_df_utm(dataframe=p1,
                                                                   column_longitude="x",
                                                                   column_latitude="y",
                                                                   utm_zone=31, hemisphere="N")
                                              p1_utm <- as.data.frame(cbind(dataframe.sp_utm@coords, dataframe.sp_utm@data))
                                              colnames(p1_utm)[1:3] <- c("x", "y", "z")

                                              # Transform entrance point in focal thermal into UTM 
                                              p2 <- data.frame(x = wantThermalFirstloc$Longitude, y = wantThermalFirstloc$Latitude, z = wantThermalFirstloc$Height.above.msl)
                                              df_to_spatial_df_utm(dataframe=p2,
                                                                   column_longitude="x",
                                                                   column_latitude="y",
                                                                   utm_zone=31, hemisphere="N")
                                              p2_utm <- as.data.frame(cbind(dataframe.sp_utm@coords, dataframe.sp_utm@data))
                                              colnames(p2_utm)[1:3] <- c("x", "y", "z")

                                              # Estimate 3D distance between thermals
                                              DistThermBefore[i] <- getDistance(p1_utm, p2_utm, DDD = TRUE)[[1]]

                                              # --
                                              
                                              # Was the selected thermal fully above the previous one?
                                              FullAbove[i] <- ifelse(wantPreviousThermalLastloc$Height.above.msl < allInd_thermalClust.df %>% 
                                                                                                                              dplyr::filter(ClusterNb == wantThermal$ClusterNb & Timestamp <= wantThermal$StartTime) %>% 
                                                                                                                              dplyr::filter(Height.above.msl == min(Height.above.msl)) %>% 
                                                                                                                              head(1) %>% pull(Height.above.msl), 1, 0)

                                              # If indexwantThermal == 1 before the individual take-off from RdA, height.above.msl = 264m, otherwise need altitude of last loc of previous thermal
                                              AltExitBefore[i] <- ifelse(indexwantThermal == 1, 264, allInd_thermalClust.df[which(allInd_thermalClust.df$IndTherm == indThermals[indexwantThermal-1,]$IndTherm),] %>% 
                                                                           tail(1) %>% 
                                                                           pull(Height.above.msl))
                                              
                                              # Individual Elo rank that year
                                              EloInd[i] <- EloRanks.ls[[as.character(year(wantThermal$Dates))]][wantThermal$IndName,"ranks"]
                                              
                                              # If the individual did not discover the thermal by itself:
                                              if (wantThermal$Discovery == 0) {
                                                
                                                # Time interval between cluster discovery and focal individual entrance in thermal
                                                discoveryTime <- orderedThermals[which(orderedThermals$ClusterNb == wantThermal$ClusterNb),c("StartTime","ClusterNb","Discovery")] %>% head(1)
                                                decisionTimeInterval <- interval(start = discoveryTime$StartTime, end = wantThermalFirstloc$Timestamp) # Time already in UTC
                                                
                                                # Need to remove 1s to the interval, otherwise at the time focal individual enter the fixe is in the decisionTimeInterval so IndAlreadyUsedThermCluster always == 1
                                                decisionTimeInterval_cor <- interval(start = discoveryTime$StartTime, end = wantThermalFirstloc$Timestamp - seconds(1)) # Time already in UTC
                                                IndAlreadyUsedThermCluster[i] <- ifelse(wantThermal$IndName %in% unique(allInd_thermalClust.df$Indname[which(allInd_thermalClust.df$Timestamp %within% decisionTimeInterval_cor & 
                                                                                                                                                               allInd_thermalClust.df$ClusterNb == wantThermal$ClusterNb)]),1,0)
                                                # Extract all locations in the thermal between its discovery and the time focal individual join it
                                                # it own locations are filtered out
                                                locWithinInterval <- allInd_thermalClust.df[which(allInd_thermalClust.df$Timestamp %within% decisionTimeInterval & 
                                                                                                    allInd_thermalClust.df$ClusterNb == wantThermal$ClusterNb &
                                                                                                    allInd_thermalClust.df$Indname != wantThermal$IndName),]
                                                if (nrow(locWithinInterval) != 0) {
                                                  
                                                  # Number of individual present in the thermal when the focal individual joined it
                                                  NbIndUsingTherm[i] <- nrow(locWithinInterval[which(locWithinInterval$Timestamp == wantThermalFirstloc$Timestamp),])
                                                  
                                                  year <- as.character(year(wantThermal$Dates))
                                                  wantIndName <- unique(wantThermal$IndName)
                                                  
                                                  # Was my preferred affiliates in th thermal (SBmax) & Max Elo difference (with individual in the thermal, or with individual that was present before)
                                                  if(NbIndUsingTherm[i] != 0){ #If individual present in the thermal with focal individual
                                                    
                                                    indInThermal <- unique(locWithinInterval[which(locWithinInterval$Timestamp == wantThermalFirstloc$Timestamp),"Indname"]) %>% pull()
                                                    
                                                    SBmax[i] <- ifelse(max(SBSMatrix.ls[[year]][wantIndName,indInThermal]) == max(SBSMatrix.ls[[year]][wantIndName,]), 1, 0) # Is my preferred affiliates in?
                                                    
                                                    # In the individual present, is there older ones?
                                                    OlderPresent[i] <- ifelse(!any(vultInfo.df[which(vultInfo.df$Name == wantThermal$IndName),"Age"] < 
                                                                                     vultInfo.df[which(vultInfo.df$Name %in% indInThermal),"Age"]),0,1)
                                                    
                                                    EloEstimate <- as_tibble(cbind(rep(EloInd[i],length(indInThermal)),EloRanks.ls[[year]][indInThermal,"ranks"])) %>% 
                                                                   dplyr::rename(., EloInd = V1, EloIndInThermal = V2) %>% 
                                                                   mutate(cubedDiff = (.$EloInd - .$EloIndInThermal)^3)
                                                    
                                                    Dt_Elo[i] <- sum(EloEstimate$cubedDiff) / nrow(EloEstimate) # Mean of the cubed difference of Elo ranks
                                                    #Dt_Elo[i] <- EloInd[i] - max(EloRanks.ls[[year]][indInThermal,"ranks"]) # difference of elo score max
                                                    
                                                  } else { # If no individual present at the time focal individual enter, but was present in the thermal before
                                                    
                                                    SBmax[i] <- 0 #max(SBSMatrix.ls[[year]][wantIndName,indInThermalPreviously]) #If no individual present not interesting to know if my preferred affiliate used it
                                                    Dt_Elo[i] <- 0 #EloInd[i] - max(EloRanks.ls[[year]][indInThermalPreviously,"ranks"]) #If no individual present a 0 will not influence the shape of the cube function
                                                  }
                                                  
                                                  # Weighted (by the number of visit) mean of the social bond with individual that used the thermal
                                                  nbIndVisits <- locWithinInterval %>% 
                                                                 group_by(Indname) %>% 
                                                                 dplyr::summarise(n = length(unique(IndTherm))) %>% 
                                                                 mutate(SBWeighted = SBSMatrix.ls[[year]][wantIndName,.$Indname] * .$n)
                                                  
                                                  SBmeanPast[i] <- sum(nbIndVisits$SBWeighted) / sum(nbIndVisits$n)
                                                  
                                                  # Number of individuals that used the thermal before the focal individual entered it
                                                  NbIndUsedTherm[i] <- length(unique(locWithinInterval$Indname))
                                                  
                                                  # Delay between arrival of the focal individual and last individual present in thermal (0 if an individual is still in the thermal)
                                                  lastTimeInThermal <- locWithinInterval %>% 
                                                                       arrange(Timestamp) %>% 
                                                                       tail(1)
                                                  
                                                  Dt_LastUser[i] <- abs(difftime(wantThermal$StartTime, lastTimeInThermal$Timestamp ,units = "secs"))

                                                  
                                                  # Max vertical speed reached by birds in the thermal between time it left its previous thermal and time at which it join the selected thermal
                                                  releaseInd_time <- as.POSIXct(paste(wantThermal$Dates , as.character(hms::as_hms(hms::as_hms(str_subset(allFlightData[[format(as_date(wantThermal$Dates), "%d:%m:%Y")]]$Day_flights[which(
                                                    allFlightData[[format(as_date(wantThermal$Dates), "%d:%m:%Y")]]$Day_flights$Flight_nb == str_sub(wantThermal$Flight, -1) & 
                                                      allFlightData[[format(as_date(wantThermal$Dates), "%d:%m:%Y")]]$Day_flights$Individual == wantThermal$IndName),"Release_UTC"], ":"))))))
                                                  
                                                  VertSpeedmax[i] <- max(locWithinInterval[which(locWithinInterval$Timestamp >= releaseInd_time),"vertSpeed_smooth"])
                                              
                                                  
                                                }
                                              }
                                            }
                                            
                                            # add variables to data.frame
                                            orderedThermals$DistThermBefore <- as.numeric(DistThermBefore)
                                            orderedThermals$AltExitBefore <- as.numeric(AltExitBefore)
                                            orderedThermals$FullAbove <- as.integer(FullAbove)
                                            orderedThermals$VertSpeedMax <- as.numeric(VertSpeedmax)
                                            orderedThermals$NbIndUsingTherm <- as.integer(NbIndUsingTherm)
                                            orderedThermals$NbIndUsedTherm <- as.integer(NbIndUsedTherm)
                                            orderedThermals$Dt_LastUser <- Dt_LastUser
                                            orderedThermals$SBmax <- as.numeric(SBmax)
                                            orderedThermals$EloInd <- as.integer(EloInd)
                                            orderedThermals$Dt_Elo <- as.integer(Dt_Elo)
                                            orderedThermals$IndAlreadyUsedThermCluster <- as.integer(IndAlreadyUsedThermCluster)
                                            orderedThermals$SBmeanPast <- as.numeric(SBmeanPast)
                                            orderedThermals$OlderPresent <- OlderPresent
                                            
                                            
                                            
                                            orderedThermals
                                          })
))


# Add Finesse value to each individual / thermals
# Remove thermals that happened when individual was flying alone
# Remove the discovery of the first thermal (cannot be something else....)
AllThermalsInfos_Filtered  <- AllThermalsInfos %>% 
  bind_cols(do.call("rbind",lapply(1:nrow(.), function(n){
    
    Finesse = ifelse(nrow(glideFiltered[glideFiltered$Date == .$Dates[n] & glideFiltered$Flight == .$Flight[n] & glideFiltered$IndName == .$IndName[n],]) == 1,
                     glideFiltered[glideFiltered$Date == .$Dates[n] & glideFiltered$Flight == .$Flight[n] & glideFiltered$IndName == .$IndName[n],"Finesse"][[1]],
                     recoverMissingFiness(glideFiltered,.$IndName[n],.$Date[n]))
    
    # Before last individual landing time
    wantFlight <- allFlightData[[format(as_date(.[n,"Dates"]), "%d:%m:%Y")]][["Day_flights"]]
    beforeLast_landingTime <- wantFlight[wantFlight$Flight_nb == str_sub(.[n,"Flight"], -1),] %>% 
      filter(.$Landing_UTC != "NA") %>% 
      arrange(.$Landing_UTC) %>%
      filter(row_number() == nrow(.)-1)
    
    # Add new column to know if yes or not this thermal was when the individual flew alone
    FlyAlone = ifelse(hms::as_hms(.[n,"StartTime"]) > hms::as_hms(beforeLast_landingTime$Landing_UTC), 1, 0)
    
    # Add new column to know if an individual is present or not in the thermal (1 if yes, 0 otherwise)
    IndPresent = ifelse(.$Dt_LastUser[n] != 0 | is.na(.$Dt_LastUser[n]) == TRUE, 0, 1)
    
    return(c(Finesse, FlyAlone, IndPresent))
  }))) %>% 
  dplyr::rename(Finesse = ...31,
                FlyAlone = ...32,
                IndPresent = ...33) %>% 
  filter(.$FlyAlone != 1) %>% # Remove all thermals that happen when individual was flying alone
  dplyr::select(-FlyAlone) %>% 
  group_by(Dates, Flight) %>% 
  filter(row_number()!=1) %>% # Remove the discovery of the first thermal (must be a discovery)
  mutate(Dt_sincefstTakeOff = abs(difftime(ymd_hms(paste0(Dates,TimeFirstTakeOff)), ymd_hms(StartTime), units = "secs"))) %>% # Add a delta time since take off, could be used as a thermal lifetime
  ungroup()


#save(AllThermalsInfos_Filtered, file="./Output/Files/AllThermalsInfos_Filtered.RData")



