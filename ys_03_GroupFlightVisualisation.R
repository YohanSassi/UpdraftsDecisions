#________________________________________________________________________________________________________________________________________________
###
###       
###       Project: Explore the balance between personal and social information in thermal updraft selection
###
###        Data come from captive vulture from the Rocher des Aigles, Rocamadour, France experiencing free-flights
###        in a 120m deep canyon.
### 
###       !!! Script used to visualize group flight data !!!
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

library(moveVis)
library(stringr)
library(lubridate)
library(move)
library(tidyr)
library(dplyr)
library(plot3D)
library(plotly)
library(RColorBrewer)





#____________________________________________________#
###
#####  Create a daily file containing all flights  ####
###
#____________________________________________________#

# Order the data collection days
daysFolderNames <- list.dirs(path = "./Data/Rocamadour_FlightData", full.names = TRUE, recursive = FALSE)
daysFolderNames <- daysFolderNames[order(dmy(daysFolderNames))] # Ordered files names

flightDates <- list.dirs(path = "./Data/Rocamadour_FlightData", full.names = FALSE, recursive = FALSE)
flightDates <- flightDates[order(dmy(flightDates))] #Ordered days of flight
flightDates <- flightDates[-length(flightDates)] #Remove the last element which will be to store cleaned datasets: "Fullday_Flights"

# Load field sheet
fieldData <- read.csv("./Data/Rocamadour_FlightData/Field2021-2022_CleanedFieldSheet.csv", header=T,
                      sep = ";", na.strings = "")
# fieldData$TakeOff_UTC <- ifelse(fieldData$TakeOff_UTC != "NA", str_sub(fieldData$TakeOff_UT,end=-1), "NA")
# fieldData$Landing_UTC <- ifelse(fieldData$Landing_UTC != "NA", str_sub(fieldData$Landing_UTC,end=-4), "NA")


# Create an empty list that will contain important information for further operations
allFlightData <- vector("list")

# Go to work iteratively in each folder corresponding to each days where flying data has been collected (2021 & 2022 data)
for (d in 1:length(flightDates)){
  
  # Collect the list of individuals files contained in this specific day
  daysIndFlightFiles <- list.files(path = daysFolderNames[d], pattern = "format.csv", all.files = FALSE,full.names = TRUE)
  
  
  # Create an empty data frames each day to concatenate individuals flights / days with reduced nb of variables
  flight_data.df <- data.frame(TagID = character(0), Date = character(0), Time = character(0), Latitude = character(0), Longitude = character(0),
                               Height.above.msl = integer(0), Ground.speed.km.h = numeric(0), Behaviour = character(0), Indname = character(0))
  
  
  for (i in 1:length(daysIndFlightFiles)){ # For each individual within a specific day
    
    # Load individual flying data & rename all columns to avoid names problems
    data.df <- read.csv(daysIndFlightFiles[i], header=T, sep = ifelse(str_detect(daysIndFlightFiles[i],"2022"),",",";") , na.strings = "")
    colnames(data.df) <- c("TagID", "Date", "Time", "CENT", "Acc.X", "Acc.Y", "Acc.Z", "Activity", "Depth", "Temperature",
                           "Latitude", "Longitude", "Height.above.msl", "Ground.speed.km.h", "Satellite.count", "Hdop",
                           "Maximum.signal.strength", "Sensor.Raw", "Behaviour")
    
    # Reduce the dataframe size to only the variables we are interested in
    data.df <- as.data.frame(data.df[which(is.na(data.df$Behaviour) == FALSE),c(1:3,11:14,19)])
    
    # Add individual names - for later to make cross measurement between individuals (cut in filename the individual name)
    underscoreindices <- str_locate_all(pattern ='_', daysIndFlightFiles[i])
    hyphenindices <- str_locate_all(pattern ='-', daysIndFlightFiles[i])
    data.df$Indname <- substr(daysIndFlightFiles[i], 
                              as.integer(underscoreindices[[1]][2,1])+1, 
                              as.integer(hyphenindices[[1]][2,1])-1)
    
    # Filter out all time with decimal seconds - because of tags with sampling method > 1Hz
    data.df <- data.df[which(str_sub(data.df$Time,start=-2) == "00"),]
    data.df <- data.df %>%
      filter(! duplicated(Time)) # Remove potential duplicated times
    
    
    # Concatenate into flight_data.df
    flight_data.df <- rbind(flight_data.df,data.df)
  }
  
  # Order by individuals and time - change Henri spelling name to match other data
  flight_data.df <- flight_data.df[order(flight_data.df$Indname, flight_data.df$Time),]
  flight_data.df$Indname <- str_replace(flight_data.df$Indname, "Henri", "Henry")
  
  # Save a Rdata file per day containing all individuals flight
  # save(flight_data.df,
  #      file=paste("./Output/Files/Fullday_Flights/FullFlights_",flightDates[d],".RData",sep=""))
  
  # Create a list with needed information for each days
  flightList <- list(Day_files = daysIndFlightFiles, 
                     Ind = unique(flight_data.df$Indname),
                     Day_flights = fieldData[which(fieldData$Date == unique(flight_data.df$Date)),])
  
  allFlightData[[flightDates[d]]] <- flightList
}

# Finally save the list data
# save(allFlightData, 
#      file="./Data/Rocamadour_FlightData/AllFlightData.info.RData")



#_________________________________________________________________#
###
##### Visualization of group flight, Rocher des Aigles  ####
###
#_________________________________________________________________#

# load files names
dailyFilesNames <- list.files(path = "./Output/Files/Fullday_Flights", full.names = TRUE, recursive = FALSE)

underscoreindices <- str_locate_all(pattern ='_', dailyFilesNames[1])
endindices <- str_locate_all(pattern ='.RData', dailyFilesNames[1])
dailyFilesNames <- dailyFilesNames[order(dmy(substr(dailyFilesNames, 
                                                    as.integer(underscoreindices[[1]][3,1])+1, 
                                                    as.integer(endindices[[1]][1,1])-1)))] # Ordered files names


load(dailyFilesNames[10]) # d # extract all flights -> flight_data.df

# Arrange time to match with move requirement (used in next script)
df <- flight_data.df %>% 
  mutate(Dates = as.character(as.Date(Date, "%d/%m/%Y"),"%Y-%m-%d")) %>% 
  mutate(Times = as.character(hms::as_hms(Time))) %>% 
  mutate(Timestamp = as.POSIXct(paste(Dates, Times), format="%Y-%m-%d %H:%M:%S")) %>% 
  dplyr::select(c(Timestamp,Latitude,Longitude,Height.above.msl,Ground.speed.km.h,Dates,Behaviour,Indname))


# ---- Interactive plot of full flight in 3D -----
unique(df$Behaviour)
want1 <- df[df$Behaviour == "Flight4",]

fig <- plot_ly(want1, 
               x = ~Longitude, 
               y = ~Latitude, 
               z = ~Height.above.msl, 
               type = 'scatter3d', 
               mode = 'lines', 
               color = ~Indname,
               colors = c("red", "green", "blue", "gold")) %>% 
  layout(title = paste0("All individuals - ",unique(want1$Dates), " - ",unique(want1$Behaviour)))

fig



# ---- Example of full flight in 3D - Fig1  -----
want1 <- want1 %>% 
  dplyr::group_by(Indname) %>% 
  dplyr::mutate(IndnameNb = cur_group_id())

par(mar = c(1,1,1,1))
scatter3D(x = want1$Longitude, 
          y = want1$Latitude, 
          z = want1$Height.above.msl,
          phi = 20, # Angle from which to see the graph
          theta = -60,
          colkey = FALSE,
          colvar = as.integer(want1$IndnameNb),
          col = brewer.pal(6, "Set2"),
          bty = "b",
          type = "l", # l for lines only
          #ticktype = "detailed", # Detailed axis
          pch = 20, # size of points
          #cex = c(2, 2, 2),
          lwd = 1,
          xlab = "Longitude",
          ylab = "Latitude",
          zlab = 'Altitude',
          main="Group flight")





#____________________________________________________#
###
##### Animation of group flight, Rocher des Aigles  ####
###
#____________________________________________________#


for (i in 1:length(dailyFilesNames)){ # Loop on files 
  
  # load data & add a Timestamp column
  load(dailyFilesNames[i])
  flight_data.df <- unite(flight_data.df,Date,Time,col = "Timestamp",sep=" ")
  
  
  # find number of flight that day
  nb_flight <- unique(flight_data.df$Behaviour)
  
  for (f in 1:length(nb_flight)){ # Loop on flights per day
    
    # Subsample by flight
    flight_dataSub.df <- flight_data.df[which(flight_data.df$Behaviour == nb_flight[f]),]
    
    # Convert our dataframe into a move object
    flight_data.m <- move(x=flight_dataSub.df$Longitude, y=flight_dataSub.df$Latitude, 
                          time=as.POSIXct(flight_dataSub.df$Timestamp,format="%d/%m/%Y %H:%M:%S", tz="UTC"), 
                          proj=CRS("+proj=longlat +ellps=WGS84"), 
                          animal=flight_dataSub.df$Indname)
    
    # Align tracks of each individuals to match in time
    m1 <- align_move(flight_data.m, res = 1, digit = "mean", unit = "secs") # Select on which temporal resolution you want tracks to be align (secs, mins, hours, days)
    
    # Create the corresponding spatial frames
    frames <- frames_spatial(m1, 
                             trace_show = FALSE,
                             path_mitre = 2,
                             equidistant = FALSE,
                             path_fade = TRUE,
                             path_colours = brewer.pal(6, "Set2"),
                             map_service = "mapbox",
                             map_token = "pk.eyJ1IjoieW9oYW5zYXNzaSIsImEiOiJjbGNnZWhwdzM3cmUxM3drZWZ6NnBtcmF2In0.TXBTp8u2Y618cXMOIDxnKQ",
                             map_type = "satellite") 
    
    #frames[[250]] #plots a single frame
    
    # Tune frames options
    frames <- frames %>% 
      add_labels(title = "Vulture group flight", 
                 caption = "Trajectory data collected at Rocher des Aigles, Rocamadour, France.
  Credit: Yohan Sassi 
  Projection: Geographic, WGS84", 
                 x = "Longitude", 
                 y = "Latitude") %>% 
      add_timestamps(type = "label") %>% 
      add_progress(colour = "white") %>% 
      add_northarrow(colour = "white", position = "bottomleft") %>% 
      add_scalebar(colour = "white", height = 0.022, position = "bottomright", label_margin = 1.4) #colour = "black", position = "bottomright", distance = 5, units = "km"
    
    #frames[[550]] #plots a single frame
    
    
    # Create temporary animation
    # animate_frames(frames,
    #                out_file = tempfile(fileext = ".mov"))

    # Create animation names
    underscoreindices <- str_locate_all(pattern ='_', dailyFilesNames[i])
    endindice <- str_locate_all(pattern ='.RData', dailyFilesNames[i])
    fileDate <- substr(dailyFilesNames[i], 
                      as.integer(underscoreindices[[1]][3,1])+1, 
                      as.integer(endindice[[1]][1,1])-1)
    
    # Create and save animation
    animate_frames(frames, 
                   display = FALSE,
                   out_file = paste("./Output/FlightAnimation/animGroupFlight_",fileDate,"_",nb_flight[f],"_",".mov",sep=""))
    
    # Free memory
    rm(frames,m1)

  }
}


