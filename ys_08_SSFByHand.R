#________________________________________________________________________________________________________________________________________________
###
###       
###       Project: Explore the balance between personal and social information in thermal updraft selection
###
###        Data come from captive vulture from the Rocher des Aigles, Rocamadour, France experiencing free-flights
###        in a 120m deep canyon.
### 
###       !!! Script used to recreate by hand a iSSF investigating how vultures choose their thermals !!!
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



library(plyr)
library(doParallel)
library(amt)
library(dplyr)
library(stringr)
library(raster)
library(lubridate)
library(jtools)
library(RColorBrewer)
library(ggplot2)


source("./Functions/df_to_spatial_df_utm.R") #To convert long/lat to UTM
source("./Functions/getDistance3D.R") #To estimate 3D distances between projected locations



#_________________________________________________________________#
###
##### Detect thermal selected vs available at the time  ####
###
#_________________________________________________________________#

# Load the list of thermals clustering by flights
thermalClust_fls <- list.files(path = "./Output/Files/ClusteredThermals", full.names = TRUE, recursive = FALSE)


# SBSMatrix list - both years combined
SBSMatrix.ls <- vector("list")
load("./Output/Files/SBSMatrix_2021_1m.RData") # 2021 social bond matrices & add capital letters to individual names (to match with data)
SBSMatrix.ls[["2021"]] <- SBSMatrix_estimatedreord
colnames(SBSMatrix.ls[["2021"]]) <- str_to_title(colnames(SBSMatrix.ls[["2021"]]))
rownames(SBSMatrix.ls[["2021"]]) <- str_to_title(rownames(SBSMatrix.ls[["2021"]]))

load("./Output/Files/SBSMatrix_2022_1m.RData") # 2022 social bond matrices
SBSMatrix.ls[["2022"]] <- SBSMatrix_estimatedreord
colnames(SBSMatrix.ls[["2022"]]) <- str_to_title(colnames(SBSMatrix.ls[["2022"]]))
rownames(SBSMatrix.ls[["2022"]]) <- str_to_title(rownames(SBSMatrix.ls[["2022"]]))

# Flight info
load("./Data/Rocamadour_FlightData/AllFlightData.info.RData") #(open list allFlightData)
allFlightData[["27:05:2021"]]$Day_flights[3,9:10] <- NA # Problem of GPS - removed otherwise induce problems in further loops

# Vulture info
vultInfo.df <- as.data.frame(cbind(c("Kazimir","Henry","Gregoire","Leon","Kirikou","Bulma","Mathilda","Hercule"),
                                   c(2014,1997,2011,2016,2014,2014,2018,2012))) %>% 
  dplyr::rename(., Name = V1, yearOfBirth = V2) %>% 
  mutate(Age = 2022 - as.integer(yearOfBirth))


# EloRanks - Hierarchy - both years combined
load("./Output/Files/EloRanks_ls.RData") #(Open list EloRanks.ls)


#--
# Load the files containing the thermals used by individuals across flights
load("./Output/Files/AllThermalsInfos_Filtered_1m.RData")
df <- AllThermalsInfos_Filtered


# Filter to keep all thermals that have been selected after being detected by another individual
chosenThermal_df <- df %>% 
  filter(Discovery == 0) %>% 
  na.omit() %>% 
  dplyr::mutate(
    chosen = 1, # Chosen = 1 will be the used thermals while chosen = 0 will represent the available
    IDbout = 1:length(chosen) # Association of an ID to each bout, will be used as step index in SSF
  ) %>% 
  dplyr::select(IndName, Dates, chosen, IDbout, Flight, ClusterNb, StartTime, Duration) %>% 
  filter(IDbout != 412) # Remove the event 412, for this event only the thermal use was discovered by conspecifics and they left before focal was released.
                          # Considered as socially unformed while not true, -> removed
  


# Extract the duration of life of each thermal, considering that we are sure they exist between the time the first individual enter it
# And the last individual exit
thermalLifeTime_df <- df %>%
  group_by(Dates, Flight, ClusterNb) %>% 
  dplyr::summarise(startLifeThermal = min(StartTime),
                  endLifeThermal = max(StartTime + Duration)) 

# Detect which thermal was available at the time individual took the considered one
# Parallelisation to reduce computational time
doParallel::registerDoParallel(detectCores()-2)

availableThermal_df <- llply(1:nrow(chosenThermal_df), function(n){ 
  #print(n)
  dfTransitory <- chosenThermal_df[n,]
  availableThermal <- thermalLifeTime_df %>% 
    filter(
        Dates == dfTransitory$Dates[1] &
        Flight == dfTransitory$Flight[1] &
        startLifeThermal < dfTransitory$StartTime[1] &
        endLifeThermal > (dfTransitory$StartTime[1] + dfTransitory$Duration[1]) &
        ClusterNb != dfTransitory$ClusterNb[1]
    ) %>% 
    mutate(
      IndName = dfTransitory$IndName[1],
      Flight = dfTransitory$Flight[1],
      Duration = dfTransitory$Duration[1],
      StartTime = as.character(dfTransitory$StartTime[1]),
      startLifeThermal = as.character(startLifeThermal),
      endLifeThermal = as.character(endLifeThermal),
      chosen = 0,
      IDbout = dfTransitory$IDbout[1]
    ) %>% 
    dplyr::select(
      IndName, Dates, chosen, IDbout, Flight, ClusterNb, StartTime, Duration
    )
  if(nrow(availableThermal) > 0){
    return(availableThermal)
  }else{
    returnTransitory <- c(rep(NA, times = 8)) %>% t() %>% as.data.frame()
    colnames(returnTransitory) <- colnames(chosenThermal_df)
    return(returnTransitory)
    
  }
}, .parallel=T)

availableThermal_df <- do.call("rbind", availableThermal_df)
availableThermal_df <- availableThermal_df[!is.na(availableThermal_df$Flight),]


## Reshaped data for SSF and add needed variables
dataForSSF <- rbind(chosenThermal_df, availableThermal_df) %>% 
              arrange(IDbout, chosen) %>%
              dplyr::rename(Chosen = chosen)



#_________________________________________________________________#
###
##### Update variables for available thermals ####
###
#_________________________________________________________________#

dataForSSF_fulldf <- llply(1:nrow(dataForSSF), function(i){
  #print(i)
  dfTransitory <- dataForSSF[i,]
  
  neededInfo_df <- dfTransitory %>% 
    bind_cols(AllThermalsInfos_Filtered %>% 
                filter(
                  IndName == dfTransitory$IndName[1] &
                    Dates == dfTransitory$Dates[1] & 
                    Flight == dfTransitory$Flight[1] &
                    StartTime == dfTransitory$StartTime[1]
                ) %>% 
                dplyr::select(DistThermBefore, AltExitBefore, VertSpeedMax, EloInd, IndAlreadyUsedThermCluster, SBmeanPast,
                              OlderPresent, FullAbove, NbIndUsedTherm, NbIndUsingTherm, SBmax, Dt_Elo, Dt_sincefstTakeOff)) 
  
  #If chosen == 0 some information need to be corrected
  if(dfTransitory$Chosen == 0){ 
    
    # Need to load the raw data 
    load(thermalClust_fls[grepl(paste0(dfTransitory$Dates,"_",dfTransitory$Flight), thermalClust_fls)])
    
    # Extract all location in the not chosen thermal before StartTime (focal individual can be inside!)
    notChosenThermalloc <- allInd_thermalClust.df[which(allInd_thermalClust.df$Timestamp < dfTransitory$StartTime & 
                                                          allInd_thermalClust.df$ClusterNb == dfTransitory$ClusterNb),]
    
    # Extract the last loc of the focal individual in the previous thermal
    lastLocPreviousTherm <- allInd_thermalClust.df %>% 
                            filter(Indname == dfTransitory$IndName) %>% 
                            filter(Timestamp < dfTransitory$StartTime) %>% 
                            filter(ClusterNb != allInd_thermalClust.df[which(allInd_thermalClust.df$Indname == dfTransitory$IndName & 
                                                                               allInd_thermalClust.df$Timestamp == dfTransitory$StartTime),"ClusterNb"][[1]]) %>% 
                            tail(1)
    
    if(nrow(lastLocPreviousTherm) == 0){ # If individual was not in another thermal before, this means it was at the Rocher des Aigles
      lastLocPreviousTherm <- data.frame(Longitude=1.611650, Latitude=44.801389, Height.above.msl = 264 ,
                                         Timestamp = as.POSIXct(paste(neededInfo_df$Dates , as.character(hms::as_hms(hms::as_hms(str_subset(allFlightData[[format(as_date(neededInfo_df$Dates), "%d:%m:%Y")]]$Day_flights[which(
                                           allFlightData[[format(as_date(neededInfo_df$Dates), "%d:%m:%Y")]]$Day_flights$Flight_nb == str_sub(neededInfo_df$Flight, -1) & 
                                             allFlightData[[format(as_date(neededInfo_df$Dates), "%d:%m:%Y")]]$Day_flights$Individual == neededInfo_df$IndName),"Release_UTC"], ":"))))))) # Lat/long/alt of the Rocher des Aigles
      neededInfo_df$AltExitBefore <- lastLocPreviousTherm$Height.above.msl # Height above msl of the Rocher des Aigles
    } else {
      
      # Otherwise change the AltExitBefore with the altitude at the exit time of the previous thermal
      neededInfo_df$AltExitBefore <- lastLocPreviousTherm$Height.above.msl
    }

    # Find the minimum distance between this last loc and the locs of individual in the not chosen thermal (required as we cannot use by definition the point of entrance of the focal ind)
    notChosenThermalloc_filt <- notChosenThermalloc %>% filter(Indname != dfTransitory$IndName) # Filtered in case individual was inside
    
    if(nrow(notChosenThermalloc_filt) != 0) { # In case individuals was the only one who used the available thermal (shouldn't be the case...)
      # neededInfo_df$DistThermBefore <- min(pointDistance(cbind(lastLocPreviousTherm$Longitude, lastLocPreviousTherm$Latitude),
      #                                                    cbind(notChosenThermalloc_filt$Longitude,notChosenThermalloc_filt$Latitude),
      #                                                    lonlat = TRUE)
      # -- FOR 3D DISTANCE
      # Transform points of paths into UTM
      p1 <- data.frame(x = lastLocPreviousTherm$Longitude, y = lastLocPreviousTherm$Latitude, z = lastLocPreviousTherm$Height.above.msl) # Last point in the previous thermal
      df_to_spatial_df_utm(dataframe=p1,
                           column_longitude="x",
                           column_latitude="y",
                           utm_zone=31, hemisphere="N")
      p1_utm <- as.data.frame(cbind(dataframe.sp_utm@coords, dataframe.sp_utm@data))
      colnames(p1_utm)[1:3] <- c("x", "y", "z")

      p2 <- data.frame(x = notChosenThermalloc_filt$Longitude, y = notChosenThermalloc_filt$Latitude, z = notChosenThermalloc_filt$Height.above.msl)
      df_to_spatial_df_utm(dataframe=p2,
                           column_longitude="x",
                           column_latitude="y",
                           utm_zone=31, hemisphere="N")
      p2_utm <- as.data.frame(cbind(dataframe.sp_utm@coords, dataframe.sp_utm@data))
      colnames(p2_utm)[1:3] <- c("x", "y", "z")

      # Estimate min 3D distance with points associated to the nonchosen thermal
      neededInfo_df$DistThermBefore <- getDistance(p1_utm, p2_utm, DDD = TRUE)[[1]]

      # Was the selected thermal fully above the previous one?
      neededInfo_df$FullAbove <- ifelse(lastLocPreviousTherm$Height.above.msl < min(notChosenThermalloc$Height.above.msl), 1, 0)

      # In the individual present, is there older ones?
      neededInfo_df$OlderPresent <- ifelse(!any(vultInfo.df[which(vultInfo.df$Name == dfTransitory$IndName),"Age"] <
                                    vultInfo.df[which(vultInfo.df$Name %in% unique(notChosenThermalloc$Indname)),"Age"]),0,1)

      
      #neededInfo_df$VertSpeedMax <- max(notChosenThermalloc$vertSpeed_smooth) # Maximum vertical speed reached by individuals in the not chosen thermal
      releaseInd_time <- as.POSIXct(paste(neededInfo_df$Dates , as.character(hms::as_hms(hms::as_hms(str_subset(allFlightData[[format(as_date(neededInfo_df$Dates), "%d:%m:%Y")]]$Day_flights[which(
        allFlightData[[format(as_date(neededInfo_df$Dates), "%d:%m:%Y")]]$Day_flights$Flight_nb == str_sub(neededInfo_df$Flight, -1) & 
          allFlightData[[format(as_date(neededInfo_df$Dates), "%d:%m:%Y")]]$Day_flights$Individual == neededInfo_df$IndName),"Release_UTC"], ":"))))))
      
      neededInfo_df$VertSpeedMax <- max(notChosenThermalloc[which(notChosenThermalloc$Timestamp >= releaseInd_time ),"vertSpeed_smooth"])
      
      # ifelse(nrow(notChosenThermalloc[which(notChosenThermalloc$Timestamp >= lastLocPreviousTherm$Timestamp),]) != 0,
      #                           max(notChosenThermalloc[which(notChosenThermalloc$Timestamp >= lastLocPreviousTherm$Timestamp),"vertSpeed_smooth"]),
      #                           max(notChosenThermalloc[which(notChosenThermalloc$Timestamp >= releaseInd_time ),"vertSpeed_smooth"]))
      # 
      # 
      
      neededInfo_df$IndAlreadyUsedThermCluster <- ifelse(dfTransitory$IndName %in% unique(notChosenThermalloc$Indname),1,0) # Did the individual already used this thermal previously
      neededInfo_df$NbIndUsedTherm <- length(which(unique(notChosenThermalloc$Indname) != dfTransitory$IndName)) # How many individuals used this thermal (except the considered individual)
      neededInfo_df$NbIndUsingTherm <- length(unique(allInd_thermalClust.df$Indname[which(allInd_thermalClust.df$Timestamp == dfTransitory$StartTime & 
                                                                                            allInd_thermalClust.df$ClusterNb == dfTransitory$ClusterNb)])) # Number of individual using the not chosen thermal
      
      # Weighted (by the number of visit) mean of the social bond with individual that used the thermal
      nbIndVisits <- notChosenThermalloc_filt %>% 
        group_by(Indname) %>% 
        dplyr::summarise(n = length(unique(IndTherm))) %>% 
        mutate(SBWeighted = SBSMatrix.ls[[as.character(year(dfTransitory$Dates))]][dfTransitory$IndName,.$Indname] * .$n)
      
      neededInfo_df$SBmeanPast <- sum(nbIndVisits$SBWeighted) / sum(nbIndVisits$n)
      
      # Extract identity of individual present in not chosen Thermal at the time focal individual enter the chosen thermal (StartTime)
      indInNotChosenThermal <- unique(allInd_thermalClust.df$Indname[which(allInd_thermalClust.df$Timestamp == dfTransitory$StartTime & 
                                                                             allInd_thermalClust.df$ClusterNb == dfTransitory$ClusterNb)])
      
      if(length(indInNotChosenThermal) == 0) { # In case no individuals are present in the not chosen thermal
        neededInfo_df$SBmax <- 0 # My preferred affiliates is thus absent
        neededInfo_df$Dt_Elo <- 0 # And there is no difference of hierarchy ranks
      } else {
        # If there are individuals in the not chosen thermal
        # Binary to know if my preferred affiliates is present in the not chosen thermal
        neededInfo_df$SBmax <- ifelse(max(SBSMatrix.ls[[as.character(year(dfTransitory$Dates))]][dfTransitory$IndName,indInNotChosenThermal]) ==
                                        max(SBSMatrix.ls[[as.character(year(dfTransitory$Dates))]][dfTransitory$IndName,]), 1, 0) 
        
        # Mean of the cubed difference of ranks between focal individual and those in the not chosen thermal
        EloEstimate <- as_tibble(cbind(rep(neededInfo_df$EloInd,length(indInNotChosenThermal)),EloRanks.ls[[as.character(year(dfTransitory$Dates))]][indInNotChosenThermal,"ranks"])) %>% 
          dplyr::rename(., EloInd = V1, EloIndInNotChosenThermal = V2) %>% 
          mutate(cubedDiff = (.$EloInd - .$EloIndInNotChosenThermal)^3)
        
        neededInfo_df$Dt_Elo <- sum(EloEstimate$cubedDiff) / nrow(EloEstimate) 
        
        #neededInfo_df$Dt_Elo <- neededInfo_df$EloInd - max(EloRanks.ls[[as.character(year(dfTransitory$Dates))]][indInNotChosenThermal,])
      }
    } else {
      neededInfo_df[,9:ncol(neededInfo_df)] <- NA # Add NA in every columns in this case to be able to filter out these cases after
    }
  }
    
return(neededInfo_df)
  
}, .parallel=T)

dataForSSF_fulldf <- do.call("rbind", dataForSSF_fulldf)
dataForSSF_fulldf <- dataForSSF_fulldf %>% na.omit()



#_________________________________________________________________#
###
##### Run step selection function  ####
###
#_________________________________________________________________#

# Filter to keep only event where 2 options at the minimum
NbEventToConsider <- dataForSSF_fulldf %>% 
  dplyr::group_by(IDbout) %>% 
  dplyr::summarise(., n = n())

nrow(NbEventToConsider[which(NbEventToConsider$n >= 2),2]) #Nb total of even we can consider


# List of IDbout for SSF
IDboutSSF <- NbEventToConsider %>% 
  filter(n >= 2) %>% 
  dplyr::select(IDbout) %>% 
  pull() 


# Scale variables
dataForSSF_fulldf$IndAlreadyUsedThermCluster <- as.factor(dataForSSF_fulldf$IndAlreadyUsedThermCluster)
dataForSSF_fulldf$SBmax <- as.factor(dataForSSF_fulldf$SBmax)
dataForSSF_fulldf$DistThermBefore <- scale(dataForSSF_fulldf$DistThermBefore, center = T, scale = T)
dataForSSF_fulldf$VertSpeedMax <- scale(dataForSSF_fulldf$VertSpeedMax, center = T, scale = T)
dataForSSF_fulldf$NbIndUsingTherm <- scale(dataForSSF_fulldf$NbIndUsingTherm, center = T, scale = T)
dataForSSF_fulldf$SBmeanPast <- scale(dataForSSF_fulldf$SBmeanPast, center = T, scale = T)
dataForSSF_fulldf$Dt_Elo <- scale(dataForSSF_fulldf$Dt_Elo, center = T, scale = T)




#### > SSF model ####
# Filter events and Fit the SSF model controlling step with IDbout (choice between thermals with and without individual in it possible)
modSSF <- dataForSSF_fulldf %>% 
  dplyr::filter(., IDbout %in% IDboutSSF) %>% 
  dplyr::rename(., step_id_ = IDbout) %>% 
  fit_clogit(Chosen ~ DistThermBefore + IndAlreadyUsedThermCluster + SBmax + 
            SBmeanPast + VertSpeedMax + NbIndUsingTherm + Dt_Elo +
          strata(step_id_),
          model = TRUE) # model with clogit link with IDbout used as step (strata)


summary(modSSF)




#### > Forest Plot SSF ####
plotSum_mod2 <- plot_summs(modSSF$model,
           ci_level = 0.95,
           colors = "grey29",
           point.shape = 19,
           coefs = c("Distance to previous thermal" = "DistThermBefore",
                     "Already used" = "IndAlreadyUsedThermCluster1",
                     "Preferred affiliate present" = "SBmax1",
                     "Mean social bond" = "SBmeanPast",
                     "Maximum vertical speed" = "VertSpeedMax",
                     "Number of conspecifics present" = "NbIndUsingTherm",
                     "Hierarchy difference" = "Dt_Elo"),
           omit.coefs = c("sd__(Intercept)", "(Intercept)")) +
  labs(x = "\n Estimates \n ") +
  scale_x_continuous(limits=c(-1.87, 2.4),
                     breaks = seq(-1.8,2.4,0.2)) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(hjust = 0.46),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) +
  ylab("")


plot_grid(plotSum_mod1,
          plotSum_mod2,
          labels=c("A","B"), 
          ncol = 1, nrow = 2,
          align = "v", axis = "tb") ## combine plot





#### > SSF prediction for significant variables ####
datMean_forSSF <- dataForSSF_fulldf %>% 
  dplyr::filter(., IDbout %in% IDboutSSF) %>% 
  dplyr::summarise_all(mean, na.rm = T)

x1 <- data.frame(NbIndUsingTherm = seq(0,5,0.1),
                 DistThermBefore = datMean_forSSF$DistThermBefore,
                 IndAlreadyUsedThermCluster = factor("0", levels = levels(dataForSSF_fulldf$IndAlreadyUsedThermCluster)),
                 SBmax = factor("0", levels = levels(dataForSSF_fulldf$SBmax)),
                 SBmeanPast = datMean_forSSF$SBmeanPast,
                 VertSpeedMax = datMean_forSSF$VertSpeedMax,
                 Dt_Elo = datMean_forSSF$Dt_Elo) 

x2 <- data.frame(NbIndUsingTherm = datMean_forSSF$NbIndUsingTherm,
                 DistThermBefore = datMean_forSSF$DistThermBefore,
                 IndAlreadyUsedThermCluster = factor("0", levels = levels(dataForSSF_fulldf$IndAlreadyUsedThermCluster)),
                 SBmax = factor("0", levels = levels(dataForSSF_fulldf$SBmax)),
                 SBmeanPast = datMean_forSSF$SBmeanPast,
                 VertSpeedMax = datMean_forSSF$VertSpeedMax,
                 Dt_Elo = datMean_forSSF$Dt_Elo)


logRSS <- log_rss(object = modSSF, x1 = x1, x2 = x2, ci = 'se', ci_level = 0.95)

# plot
ggplot(logRSS$df, aes(x = NbIndUsingTherm_x1, y = exp(log_rss),
                      ymin = exp(lwr), ymax = exp(upr))) +
  geom_ribbon(color = 'black', fill = 'grey92', linetype = 'dashed') +
  geom_line() +
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') +
  xlab("Number of individual using the thermal") +
  ylab("RSS vs mean number of individual using thermals") + 
  theme_bw() +
  theme(
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white")) 

# Results means e.g.: when having the choice, an individual is 25 times [range 7x, 125x] more likely to choose a 
# thermal with 5 individuals compared to a thermal with only 1. 











#### > SSF model with presence only ####
# Ssf model with choice between thermal where individuals are present only
NbEventToConsider2 <- dataForSSF_fulldf %>% 
                      dplyr::filter(NbIndUsingTherm != 0) %>% 
                      dplyr::group_by(IDbout) %>% 
                      dplyr::summarise(., n = n(), chosen = length(unique(Chosen)))

IDboutSSF2 <- NbEventToConsider2 %>% 
              filter(n >= 2 & chosen == 2) %>% 
              dplyr::select(IDbout) %>% 
              pull()
  
modSSF2 <- dataForSSF_fulldf %>% 
        dplyr::filter(., IDbout %in% IDboutSSF2) %>% # Bout where at least another option was available
        dplyr::rename(., step_id_ = IDbout) %>% 
        fit_clogit(Chosen ~ DistThermBefore + IndAlreadyUsedThermCluster + SBmax +
                          SBmeanPast + VertSpeedMax + NbIndUsingTherm + Dt_Elo +
                          strata(step_id_),
                          model = TRUE) # model with clogit link with IDbout used as step (strata)

  
summary(modSSF2)



### Fit checks ---
library(survival)
library(survminer)
library(performance)

check_collinearity(modSSF$model)
check_model(modSSF$model, fittedWith = "glmmTMB")

# Same model but fitted with coxph function to be able to use the output for ssf verifications
res.cox <- dataForSSF_fulldf %>% 
  dplyr::filter(., IDbout %in% IDboutSSF) %>% 
  survival::coxph(formula = Surv(rep(1, 409), Chosen) ~ DistThermBefore + 
                    IndAlreadyUsedThermCluster + SBmax + SBmeanPast + VertSpeedMax + 
                    NbIndUsingTherm + Dt_Elo + strata(IDbout), data = ., model = TRUE)

summary(res.cox)
confint(res.cox)

    -0.19# 
test.ph <- survival::cox.zph(res.cox) # Not working yet... 
ggcoxzph(test.ph)


# Diagnose influential observations
# Specifying the argument type = “dfbeta”, plots the estimated changes in the regression coefficients upon deleting each observation in turn
ggcoxdiagnostics(res.cox, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
# Plots show that comparing the magnitudes of the largest dfbeta values to the regression coefficients suggests that none of the observations is terribly influential individually

# Check outliers by visualizing the deviance residuals. The deviance residual is a normalized transform of the martingale residual. 
# These residuals should be roughtly symmetrically distributed about zero with a standard deviation of 1.
ggcoxdiagnostics(res.cox, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())
# The pattern looks fairly symmetric around 0.


# Test for non-linearity
# Nonlinearity is not an issue for categorical variables, so we only examine plots of martingale residuals and partial residuals against a continuous variable.
# res.cox2 <- dataForSSF_fulldf %>% 
#   dplyr::filter(., IDbout %in% IDboutSSF) %>% 
#   survival::coxph(formula = Surv(rep(1, 412L), Chosen) ~ DistThermBefore + SBmeanPast + VertSpeedMax + Dt_Elo + NbIndUsingTherm, data = .)
# 
# dataForSSF_fulldf %>% 
#   dplyr::filter(., IDbout %in% IDboutSSF) %>% 
#   survminer::ggcoxfunctional(res.cox2, data = .) #SBmeanPast + VertSpeedMax + NbIndUsingTherm + Dt_Elo + strata(IDbout)




#### > Forest plot SSF sensitivity analysis on events ####
plot_summs(modSSF$model,
           modSSF2$model,
           ci_level = 0.95,
           point.shape = 19,
           colors = brewer.pal(3, "Set2"),
           model.names = c("All events (n = 178)","With presence only (n = 61)"),
           coefs = c("Distance to previous thermal" = "DistThermBefore",
                     "Already used" = "IndAlreadyUsedThermCluster1",
                     "Preferred affiliate present" = "SBmax1",
                     "Mean social bond" = "SBmeanPast",
                     "Maximum vertical speed" = "VertSpeedMax",
                     "Number of conspecifics present" = "NbIndUsingTherm",
                     "Hierarchy difference" = "Dt_Elo")) + 
  labs(x = "\n Estimates \n ") +
  scale_x_continuous(limits=c(-1.7,+1.7),
                     breaks = seq(-1.5,1.5,0.5)) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) +
  ylab("")








#### > Forest plot SSF sensitivity analysis on social bond strengths ####
plot_summs(modSSF$model,
           modSSF_1.3m$model,
           modSSF_1m$model,
           ci_level = 0.95,
           point.shape = 19,
           colors = brewer.pal(3, "Set2"),
           model.names = c("Threshold of 1.55m","Threshold of 1.3m","Threshold of 1m"),
           coefs = c("Distance to previous thermal" = "DistThermBefore",
                     "Already used" = "IndAlreadyUsedThermCluster1",
                     "Preferred affiliate present" = "SBmax1",
                     "Mean social bond" = "SBmeanPast",
                     "Maximum vertical speed" = "VertSpeedMax",
                     "Number of conspecifics present" = "NbIndUsingTherm",
                     "Hierarchy difference" = "Dt_Elo")) + 
  labs(x = "\n Estimates \n ") +
  scale_x_continuous(limits=c(-1.7,+1.7),
                     breaks = seq(-1.5,1.5,0.5)) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) +
  ylab("")


