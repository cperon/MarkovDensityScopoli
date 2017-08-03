
rm(list=ls())

setwd('C:/Users/Clara PERON/Documents/PELAGIC2/1-Methods/Shearwaters/gps_integration-master/markov_density/')

require(Matrix) # to create sparse matrices
library(matlab)
library(raster)
library(RColorBrewer)
library(CLmapping)
library(fields)
library(gdata)
source("markovdens.R")
source("secondtest.R")

# Open datasets
  setwd('C:/Users/Clara PERON/Documents/PELAGIC2/1-Methods/Shearwaters/Tracking/')

  load('results/GAM_models/2012/presence_absence/GPS_SIM_12min.Rdata')

# Select colonies  
  GPS <- GPS_SIM_12min[GPS_SIM_12min$Site %in% c('Riou'),]

# Keep only observed tracks
  GPS <- GPS[is.na(GPS$ID_Sim)==T,]
  
# Simplify dataset  
  dat <- GPS[, c('lon', 'lat', 'trip_ID')]
  dat$trip_ID <- as.factor(as.character(dat$trip_ID))

  data <- NULL  
  for(i in levels(dat$trip_ID)){
  sub <- dat[dat$trip_ID==i,]
  sub$ind <- 1:nrow(sub)
  data <- rbind(data, sub)
  }
  
  ID <- unique(data$trip_ID)
  
  ID <- as.character(ID[1:10])
  
  data <- data[data$trip_ID %in% ID,]
  
  
  names(data) <- c('x', 'y', 'id', 'ind')
  head(data)
  save(data, file='C:/Users/Clara PERON/Documents/PELAGIC2/1-Methods/Shearwaters/gps_integration-master/markov_density/Shearw_data.Rdata')
  
 