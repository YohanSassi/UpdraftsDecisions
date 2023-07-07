#________________________________________________________________________________________________________________________________________________
###
###       
###       Project: Explore the balance between personal and social information in thermal updraft selection
###
###        Data come from captive vulture from the Rocher des Aigles, Rocamadour, France experiencing free-flights
###        in a 120m deep canyon.
### 
###       !!! Script used to fit models investigating how vultures choose their thermals !!!
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

library(glmmTMB)
library(DHARMa)
library(car)
library(MuMIn)
library(dplyr)
library(sjPlot)
library(ggplot2)
library(plotly)
library(lubridate)
library(cowplot)
library(ggeffects)
library(jtools)
library(emmeans)


source("./Functions/Check_Model.R")




#____________________________________________________________#
###
##### Exploration of the data ####
###
#____________________________________________________________#


load("./Output/Files/AllThermalsInfos_Filtered.RData")

str(AllThermalsInfos_Filtered) # Structure of the variables


# Change variable types
AllThermalsInfos_Filtered$Duration <- as.integer(AllThermalsInfos_Filtered$Duration)
AllThermalsInfos_Filtered$Dt_Discovery <- abs(AllThermalsInfos_Filtered$Dt_Discovery)
AllThermalsInfos_Filtered$Age <- as.numeric(AllThermalsInfos_Filtered$Age)
AllThermalsInfos_Filtered$Temp <- as.numeric(AllThermalsInfos_Filtered$Temp)
AllThermalsInfos_Filtered$Nebulosity <- as.numeric(AllThermalsInfos_Filtered$Nebulosity)
AllThermalsInfos_Filtered$ReleaseOrder <- as.factor(AllThermalsInfos_Filtered$ReleaseOrder)
AllThermalsInfos_Filtered$AltEntrance <- as.numeric(AllThermalsInfos_Filtered$AltEntrance)
AllThermalsInfos_Filtered$Dt_sincefstTakeOff <- as.numeric(AllThermalsInfos_Filtered$Dt_sincefstTakeOff)
AllThermalsInfos_Filtered$IndAlreadyUsedThermCluster <- as.factor(AllThermalsInfos_Filtered$IndAlreadyUsedThermCluster)



# for(i in 1:nrow(AllThermalsInfos_Filtered)){
#   AllThermalsInfos_Filtered$EloInd[i] <- EloRanks.ls[[as.character(year(AllThermalsInfos_Filtered[i,"Dates"][[1]]))]][AllThermalsInfos_Filtered[i,"IndName"][[1]],"ranks"]
# }


# Number of thermals discovered with personal information VS used thanks to social information (TOTAL)
AllThermalsInfos_Filtered %>% 
  group_by(year(.$Dates)) %>% 
  summarise(nbThermDiscovered = sum(Discovery == 1),
            nbThermSocial = sum(Discovery == 0),
            nTot = n())


# Frequency at which individuals detect thermals by their own
AllThermalsInfos_Filtered  %>% 
  group_by(IndName) %>% 
  summarise(nDiscovery = sum(Discovery == 0),
                   ntot = n(),
                   freqDiscovery = nDiscovery / n())


# Hierarchy by years
load("./Output/Files/EloRanks_ls.RData")
EloRanks.ls


# Exploration of Finesse
load("./Output/Files/GlideToKeep_df.RData")
head(glideToKeep_Allvariables.df)

glideFiltered <- glideToKeep_Allvariables.df %>% 
  mutate(Date_Flight = paste0(Date,"_",Flight)) %>% 
  filter(Straightness > 0.95) %>% # Consider only straight lines
  group_by(Date_Flight, IndName) %>%  #Date, Flight, IndName
  summarise(n = n(),
            Max_Finesse = max(Finesse),
            Min_Finesse = min(Finesse),
            Mean_Finesse = mean(Finesse))

plot1 <- ggplot(glideFiltered, aes(glideFiltered$Date_Flight, glideFiltered$Max_Finesse, group = glideFiltered$IndName)) + 
  geom_point(aes(col = glideFiltered$IndName)) +
  geom_line(aes(col = glideFiltered$IndName)) +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(y = "Max Finesse", x = "Date_flight", color = "Individuals")

ggplotly(plot1)



# Thermal duration distribution by individuals
plot2 <- ggplot(AllThermalsInfos_Filtered, aes(x = Duration, fill = IndName)) + 
  geom_density(alpha = 0.4) +
  labs(x = "Time spent in thermal [sec]", fill = "Individuals") +
  geom_vline(xintercept = 42, linetype = "dashed", alpha = 0.6)

# test <- AllThermalsInfos_Filtered %>% 
#   group_by(Duration) %>% 
#   summarise(n = n()) %>% 
#   filter(n == max(n))

ggplotly(plot2)

plot3 <- ggplot(AllThermalsInfos_Filtered, aes(x=IndName, y=Duration, fill=IndName)) +
  geom_violin(width=1.2, alpha = 0.7) +
  geom_boxplot(width=0.1, color="black", alpha=0.8) +
  theme_bw() +
  theme(
    legend.position="none",
    axis.text.x = element_text(size=12),
  ) + 
  xlab("class") +
  labs(y = "Time spent in thermal [sec]", x = "")
ggplotly(plot3)



# Thermal duration distribution in function of Age
plot4 <- ggplot(AllThermalsInfos_Filtered, aes(x=as.factor(Age), y=Duration, fill = as.factor(Age))) +
  geom_violin(width=1.2, alpha = 0.7) +
  geom_boxplot(width=0.1, color="black", alpha=0.8) +
  theme_bw() +
  theme(
    legend.position="none",
    axis.text.x = element_text(size=12),
  ) + 
  xlab("class") +
  labs(y = "Time spent in thermal [sec]", x = "Age")


# Time spent in a thermal in function of the maximal vertical speed that have been reached in it 
# (by all individual that used it before the focal enter)
AllThermalsInfos_Filtered %>% 
  select(Duration,VertSpeedMax) %>% 
  na.omit() %>% 
  ggplot(.,aes(x = VertSpeedMax, y = Duration)) + 
  geom_point() + 
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  theme_bw() +
  labs(y = "Time spent in thermal [sec]", x = "Max vertical speed reached in thermal")
  



# Time spent in a thermal in function of the number of individual using it at the time the individual joined it
AllThermalsInfos_Filtered %>% 
  select(Duration,NbIndUsingTherm) %>% 
  na.omit() %>% 
  ggplot(.,aes(x = NbIndUsingTherm, y = Duration)) + 
  geom_point() + 
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  theme_bw() +
  labs(y = "Time spent in thermal [sec]", x = "NbIndUsingTherm")


AllThermalsInfos_Filtered %>% 
  select(Duration,NbIndUsingTherm) %>% 
  na.omit() %>%
  ggplot(.,aes(x = NbIndUsingTherm, y = Duration, fill=as.factor(NbIndUsingTherm))) +
  geom_violin(width=1.2, alpha = 0.7) +
  geom_boxplot(width=0.1, color="black", alpha=0.8) +
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  theme_bw() +
  theme(
    legend.position="none",
    axis.text.x = element_text(size=12),
  ) + 
  labs(y = "Time spent in thermal [sec]", x = "Number of individual present in thermal when joining")


# test <- AllThermalsInfos_Filtered %>% 
#   select(Duration,NbIndUsingTherm) %>% 
#   na.omit()
# 
# summary(lm(test$Duration ~ test$NbIndUsingTherm))


# Time spent in the thermal depending on the altitude at which the individual entered
AllThermalsInfos_Filtered %>% 
  select(Duration,AltEntrance, IndName) %>% 
  na.omit() %>% 
  #filter(AltEntrance < 500) %>% 
  ggplot(.,aes(x = AltEntrance, y = Duration, color = IndName)) + 
  geom_point() + 
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  theme_bw() +
  labs(y = "Time spent in thermal [sec]", x = "Altitude entrance in thermal")


# Check correlation between altEntrance & AltExitBefore
AllThermalsInfos_Filtered %>% 
  select(AltExitBefore,AltEntrance, IndName) %>% 
  na.omit() %>% 
  ggplot(.,aes(x = AltEntrance, y = AltExitBefore)) + 
  geom_point() + 
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_bw() +
  labs(y = "Altitude at the exit of the previous thermal", x = "Altitude entrance in thermal")


# Check correlation between Nebulosity and temperature
AllThermalsInfos_Filtered %>% 
  select(Nebulosity,Temp) %>% 
  na.omit() %>% 
  ggplot(.,aes(x = Temp, y = Nebulosity)) + 
  geom_point() + 
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_bw() +
  labs(y = "Nebulosity", x = "Temperature")

cor(AllThermalsInfos_Filtered$Nebulosity, AllThermalsInfos_Filtered$Temp)

##--




# histogram of variables of interests
library(Hmisc)
hist.data.frame(AllThermalsInfos_FilteredForSocial[,c("DistThermBefore","AltExitBefore","NbIndUsingTherm","EloInd","Finesse","Age","",)])









#____________________________________________________________________________________________________#
###
#### Model to predict the probability of discovering thermal - both year combined ####
###
#____________________________________________________________________________________________________#


# Scale all numeric variable
AllThermalsInfos_Filtered$DistThermBefore <- scale(AllThermalsInfos_Filtered$DistThermBefore, center = TRUE, scale= TRUE)
AllThermalsInfos_Filtered$AltExitBefore <- scale(log(AllThermalsInfos_Filtered$AltExitBefore), center = TRUE, scale= TRUE)
AllThermalsInfos_Filtered$NbIndUsingTherm <- scale(AllThermalsInfos_Filtered$NbIndUsingTherm, center = TRUE, scale= TRUE)
AllThermalsInfos_Filtered$EloInd <- scale(AllThermalsInfos_Filtered$EloInd, center = TRUE, scale= TRUE)
AllThermalsInfos_Filtered$Finesse <- scale(AllThermalsInfos_Filtered$Finesse, center = TRUE, scale= TRUE)
AllThermalsInfos_Filtered$Age <- scale(AllThermalsInfos_Filtered$Age, center = TRUE, scale= TRUE)
AllThermalsInfos_Filtered$Nebulosity <- scale(AllThermalsInfos_Filtered$Nebulosity, center = TRUE, scale= TRUE)
AllThermalsInfos_Filtered$Temp <- scale(AllThermalsInfos_Filtered$Temp, center = TRUE, scale= TRUE)
AllThermalsInfos_Filtered$Dt_sincefstTakeOff <- scale(AllThermalsInfos_Filtered$Dt_sincefstTakeOff, center = TRUE, scale= TRUE)



## MODEL
AllThermalsInfos_Filtered <- AllThermalsInfos_Filtered %>% 
                              mutate(Year = year(Dates)) %>% # Add year for random effect
                              mutate(Discovery = ifelse(Discovery == 1, 0, 1)) %>% # Change for discovery == 1 when individual join other -> simplification for the paper
                              mutate(WindSpeed = relevel(as.factor(WindSpeed), "Nul"))

mod_DiscoverTherm <- glmmTMB(Discovery ~ Nebulosity + WindSpeed + Temp + EloInd + Age + Finesse + 
                      AltExitBefore + DistThermBefore + ReleaseOrder + Dt_sincefstTakeOff + (1|IndName),
                    family = "binomial",
                    data = AllThermalsInfos_Filtered)

summary(mod_DiscoverTherm)


### Fit check ---
modDrop <- drop1(mod_DiscoverTherm, test = "Chisq")

modelDiag <- diagnostics.plot.dharma(mod_DiscoverTherm)
modelCheck <- check_model(mod_DiscoverTherm,
                           fittedWith = "glmmTMB")

### Summary tables ---
sjPlot::tab_model(mod_DiscoverTherm,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")


### Estimation of effect size [ + ESTIMATE ODD RATIO]
#-- 
ef_temp <- emmeans(mod_DiscoverTherm, "Temp", type = "response", at = list(Temp = c(min(AllThermalsInfos_Filtered$Temp),
                                                                 max(AllThermalsInfos_Filtered$Temp))))
pairs(ef_temp, type="response")

# --
ef_dist <- emmeans(mod_DiscoverTherm, "DistThermBefore", type = "response", at = list(DistThermBefore = c(min(AllThermalsInfos_Filtered$DistThermBefore),
                                                                                                          max(AllThermalsInfos_Filtered$DistThermBefore))))
pairs(ef_dist, type="response")

# --
ef_alt <- emmeans(mod_DiscoverTherm, "AltExitBefore", type = "response", at = list(AltExitBefore = c(min(AllThermalsInfos_Filtered$AltExitBefore),
                                                                                                     max(AllThermalsInfos_Filtered$AltExitBefore))))
pairs(ef_alt, type="response")


# --
ef_elo <- emmeans(mod_DiscoverTherm, "EloInd", type = "response", at = list(EloInd = c(min(AllThermalsInfos_Filtered$EloInd),
                                                                                       max(AllThermalsInfos_Filtered$EloInd))))
pairs(ef_elo, type="response")



#### > Model forest plot ####
plotSum_mod1 <- plot_summs(mod_DiscoverTherm,
           ci_level = 0.95,
           colors = "grey29",
           point.shape = 19,
           coefs = c("Cloudiness" = "Nebulosity",
                     "Wind speed high" = "WindSpeedFort",
                     "Wind speed medium" = "WindSpeedMoyen",
                     "Wind speed low" = "WindSpeedFaible",
                     "Temperature" = "Temp",
                     "Rank" = "EloInd",
                     "Age" = "Age",
                     "Glide-ratio" = "Finesse",
                     "Exit altitude from previous thermal" = "AltExitBefore",
                     "Distance to previous thermal" = "DistThermBefore",
                     "Release order" = "ReleaseOrder2",
                     "Time since 1st take-off" = "Dt_sincefstTakeOff")) +
  scale_x_continuous(limits=c(-1.87, 2.4),
                     breaks = seq(-1.8,2.4,0.2)) +
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
  ylab("") + xlab("")






#### > Significant variables prediction plot ####
forProba <- data.frame(bind = c(1000,2000, 3000, 4000, 5000))
forProba$proba1 <- rep(NA,nrow(forProba))
forProba$ntot1 <-rep(NA,nrow(forProba))
forProba$Alt <- c(500, 750, 1000, 1250, 1500)
forProba$proba2 <- rep(NA,nrow(forProba))
forProba$ntot2 <-rep(NA,nrow(forProba))
forProba$Temp <- c(20,23,26,29,31)
forProba$proba3 <- rep(NA,nrow(forProba))
forProba$ntot3 <-rep(NA,nrow(forProba))
  
for(i in 1:nrow(forProba)){
  if (i == 1){
    dat <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$DistThermBefore <= 1000),"Discovery"]
    forProba$proba1[i] <- nrow(dat[which(dat$Discovery == 1),]) / nrow(dat)
    forProba$ntot1[i] <- nrow(dat)
    
    dat2 <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$AltExitBefore <= 500),"Discovery"]
    forProba$proba2[i] <- nrow(dat2[which(dat2$Discovery == 1),]) / nrow(dat2)
    forProba$ntot2[i] <-nrow(dat2)
    
    dat3 <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$Temp <= 20),"Discovery"]
    forProba$proba3[i] <- nrow(dat3[which(dat3$Discovery == 1),]) / nrow(dat3)
    forProba$ntot3[i] <-nrow(dat3)
    
  } else if (i == nrow(forProba)){
    dat <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$DistThermBefore >= 5000),"Discovery"]
    forProba$proba1[i] <- nrow(dat[which(dat$Discovery == 1),]) / nrow(dat)
    forProba$ntot1[i] <- nrow(dat)
    
    dat2 <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$AltExitBefore >= 1500),"Discovery"]
    forProba$proba2[i] <- nrow(dat2[which(dat2$Discovery == 1),]) / nrow(dat2)
    forProba$ntot2[i] <- nrow(dat2)
    
    dat3 <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$Temp > 29),"Discovery"]
    forProba$proba3[i] <- nrow(dat3[which(dat3$Discovery == 1),]) / nrow(dat3)
    forProba$ntot3[i] <- nrow(dat3)
    
  } else {
    dat <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$DistThermBefore > forProba$bind[i-1] & 
                                             AllThermalsInfos_Filtered$DistThermBefore < forProba$bind[i]),"Discovery"]
    forProba$proba1[i] <- nrow(dat[which(dat$Discovery == 1),]) / nrow(dat)
    forProba$ntot1[i] <- nrow(dat)
    
    dat2 <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$AltExitBefore > forProba$Alt[i-1]& 
                                              AllThermalsInfos_Filtered$AltExitBefore < forProba$Alt[i]),"Discovery"]
    forProba$proba2[i] <- nrow(dat2[which(dat2$Discovery == 1),]) / nrow(dat2)
    forProba$ntot2[i] <- nrow(dat2)
    
    dat3 <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$Temp > forProba$Temp[i-1]& 
                                              AllThermalsInfos_Filtered$Temp <= forProba$Temp[i]),"Discovery"]
    forProba$proba3[i] <- nrow(dat3[which(dat3$Discovery == 1),]) / nrow(dat3)
    forProba$ntot3[i] <- nrow(dat3)
  }
}
    


# Predict Temperature --- 
mydf <- ggemmeans(mod_DiscoverTherm, ci.lvl = 0.95, terms = "Temp")

plotTemp <- ggplot() +
  geom_line(data = mydf, aes(x, predicted)) + 
  geom_ribbon(data = mydf, aes(x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(aes(x = Temp, y = proba3, size = ntot3), data = forProba, color = "grey49") +
  xlab("Temperature [°C]") + 
  ylab("Probability to use a thermal 
previously discovered") +  
  expand_limits(x = 15, y = 0) + 
  ylim(0,1) + 
  theme_light() +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) 





# Predict DistThermBefore ---
mydf <- ggemmeans(mod_DiscoverTherm, ci.lvl = 0.95, terms = "DistThermBefore[all]")

plotDistToThermal <- ggplot() +
  geom_line(data = mydf, aes(x, predicted)) + 
  geom_ribbon(data = mydf, aes(x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(aes(x = bind, y = proba1, size = ntot1), data = forProba, color = "grey49") +
  xlab("Distance to the previous thermal [m]") + 
  ylab("") +  
  expand_limits(x = 0, y = 0) + 
  ylim(0,1) + 
  theme_light() +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) 



# Predict AltExitBefore ---
mydf <- ggemmeans(mod_DiscoverTherm, ci.lvl = 0.95, terms = "AltExitBefore[all]")

plotAltPrevThermal <- ggplot() +
  geom_line(data = mydf, aes(x, predicted)) +
  geom_ribbon(data = mydf, aes(x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(aes(x = Alt, y = proba2, size = ntot2), data = forProba, color = "grey49") + 
  xlab("Previous thermal exit altitude [m]") + 
  ylab("Probability to use a thermal 
previously discovered") +  
  expand_limits(x = 0, y = 0) + 
  ylim(0,1) + 
  theme_light() +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) 






forRanks <- data.frame(rank = 1:6)
forRanks$proba1 <- rep(NA,nrow(forRanks))
forRanks$ntot1 <- rep(NA,nrow(forRanks))
for(i in 1:nrow(forRanks)){
    dat <- AllThermalsInfos_Filtered[which(AllThermalsInfos_Filtered$EloInd == forRanks$rank[i]),"Discovery"]
    forRanks$proba1[i] <- nrow(dat[which(dat$Discovery == 1),]) / nrow(dat)
    forRanks$ntot1[i] <- nrow(dat)
}

# Predict EloInd ---
mydf <- ggemmeans(mod_DiscoverTherm, ci.lvl = 0.95, terms = "EloInd")

plotHierarchy <- ggplot() +
  geom_line(data = mydf, aes(x, predicted)) + 
  geom_ribbon(data = mydf, aes(x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_point(aes(x = rank, y = proba1, size = ntot1), data = forRanks, color = "grey49") + 
  xlab("Dominance hierarchy rank") + 
  ylab("") +  
  expand_limits(x = 0, y = 0) + 
  ylim(0,1) + 
  theme_light() +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) 




### Combined plot of prediction ---
plot_grid(plotTemp, plotDistToThermal, 
          plotAltPrevThermal, plotHierarchy,
          labels=c("A","B","C","D"), 
          ncol = 2, nrow = 2) ## combine plot







#____________________________________________________________________________________________________#
###
#### Model to predict the probability of discovering thermal - developed on 2021 and fitted on 2022  ####
###
#____________________________________________________________________________________________________#

# Split data by years
df_2021 <- AllThermalsInfos_Filtered %>%
          filter(year(Dates) == 2021)

df_2022 <- AllThermalsInfos_Filtered %>%
  filter(year(Dates) == 2022)


# Models
mod_DiscoverTherm2021 <- glmmTMB(Discovery ~ Nebulosity + WindSpeed + Temp + EloInd + Age + Finesse + 
                                   AltExitBefore + DistThermBefore + ReleaseOrder + Dt_sincefstTakeOff + (1|IndName),
                                 family = "binomial",
                                 data = df_2021)

summary(mod_DiscoverTherm2021)



mod_DiscoverTherm2022 <- glmmTMB(Discovery ~ Nebulosity + WindSpeed + Temp + EloInd + Age + Finesse + 
                                   AltExitBefore + DistThermBefore + ReleaseOrder + Dt_sincefstTakeOff + (1|IndName),
                                 family = "binomial",
                                 data = df_2022)

summary(mod_DiscoverTherm2022)


### Fit check ---
modDrop <- drop1(mod_DiscoverTherm2021, test = "Chisq")
modDrop <- drop1(mod_DiscoverTherm2022, test = "Chisq")



modelDiag <- diagnostics.plot.dharma(mod_DiscoverTherm2021)
modelCheck <- check_model(mod_DiscoverTherm2021,
                          fittedWith = "glmmTMB")

modelDiag <- diagnostics.plot.dharma(mod_DiscoverTherm2022)
modelCheck <- check_model(mod_DiscoverTherm2022,
                          fittedWith = "glmmTMB")

### Summary tables ---
sjPlot::tab_model(mod_DiscoverTherm2021,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")

sjPlot::tab_model(mod_DiscoverTherm2022,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")


#### > Forest plot model 1 fitted on 2021 / 2022 data ####

plot_summs(mod_DiscoverTherm,
           mod_DiscoverTherm2021,
           mod_DiscoverTherm2022,
           ci_level = 0.95,
           point.shape = 19,
           colors = brewer.pal(3, "Set2"),
           omit.coefs = c("sd__(Intercept)", "(Intercept)"),
           model.names = c("Both years combined","2021 data only","2022 data only"),
           coefs = c("Cloudiness" = "Nebulosity",
                     "Wind speed high" = "WindSpeedFort",
                     "Wind speed medium" = "WindSpeedMoyen",
                     "Wind speed low" = "WindSpeedFaible",
                     "Temperature" = "Temp",
                     "Rank" = "EloInd",
                     "Age" = "Age",
                     "Glide-ratio" = "Finesse",
                     "Exit altitude from previous thermal" = "AltExitBefore",
                     "Distance to previous thermal" = "DistThermBefore",
                     "Release Order" = "ReleaseOrder2",
                     "Time since 1st take-off" = "Dt_sincefstTakeOff")) + 
  labs(x = "\n Estimates \n ") +
  scale_x_continuous(limits=c(-1.8,+2.4),
                     breaks = seq(-2.6,2.6,0.2)) +
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





# plot_model(mod_SocialTherm,
#            colors = "bw",
#            vline.color = "black",
#            sort.est = TRUE,
#            show.intercept = FALSE,
#            show.p = TRUE,
#            show.values = TRUE, 
#            value.offset = .15,
#            p.style = "asterisk",
#            p.threshold = c(0.05, 0.01, 0.001),
#            line.size = 0.5) +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", size = 16),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 12, face = "bold"),
#         panel.grid.minor = element_line(colour = "grey93"),
#         panel.grid.major = element_line(colour = "grey93"),
#         strip.background = element_rect(colour = "white",
#                                         fill = "white"),
#         strip.text = element_text(face = "bold", size = 14))


