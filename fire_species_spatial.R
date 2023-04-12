# Created: April 11, 2023
# Modified: 

#### This script plots fire-experiencing species on a map to visualize their distributions in the study area. Are they spatially correlated? How do distributions differ on western vs eastern slope? 


### Load libraries
library(tidyverse)
library(raster)
library(maptools)
library(rgdal)
library(rgeos)
library(sf)
library(maps)
library(mapdata)
library(mapproj)
library(ggstar)
library(geodata)

### Read in data
# Plot survey data
dat <- read_csv("data/1_presence_fires_unrarefied.csv") %>% #unfortunately the plot names in this file in column Plot differ between legacy and resurvey time points
  mutate(Plot.Name.2015 = ifelse(Data.Type=="Resurvey", Plot, NA),
         Plot.Name.1980 = ifelse(Data.Type=="Legacy", Plot, NA))
dat.leg <- dat %>% filter(Data.Type=="Legacy")
dat.res <- dat %>% filter(Data.Type=="Resurvey")

# Plot locations and fire status
plots <- read_csv("data/0_Lat.Long.csv") %>% 
  drop_na() %>% 
  mutate(Longitude=ifelse(Longitude>0, -Longitude, Longitude)) #input file accidentally has some longitudes in E instead of W
colnames(plots)[2] = "Plot.Name.2015"
colnames(plots)[3] = "Plot.Name.1980"
plots$Plot.Name.1980 = as.character(plots$Plot.Name.1980)

# Join plot location and fire info to resurvey data
dat.leg <- left_join(dat.leg, plots, by=c("Plot"="Plot.Name.1980")) %>% #NOTE: this join is adding 390 rows -- TROUBLESHOOT! 
  dplyr::select(-Plot, -Plot.Name.2015.x)
colnames(dat.leg)[8] = "Plot.Name.2015"
dat.res <- left_join(dat.res, plots, by=c("Plot"="Plot.Name.2015")) %>% #this join is good
  dplyr::select(-Plot, -Plot.Name.1980.x) 
colnames(dat.res)[8] = "Plot.Name.1980"

# Re-merge time points into one frame
dat.all <- bind_rows(dat.leg, dat.res) %>% # this now has plot locations and plot names translated from legacy <--> resurvey
  drop_na() # Latitude=NA for legacy 4042, 4063, resurvey Dia234, Thor221 because these are missing from plots file

dat.pres <- filter(dat.all, Pres.Abs==1)

ACMI <- filter(dat.pres, Species.Code=="ACMI")

### Set up maps
prj.wgs <- "+proj=longlat + type=crs"
prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs"

coordinates(dat.pres) <- ~Longitude+Latitude #convert to spatial data
projection(dat.pres) <- CRS('+proj=longlat') #define projection
dat.pres <- spTransform(dat.pres, CRS=CRS(prj.wgs)) #transform projection 
dat.pres.lcc <- spTransform(dat.pres, CRS=CRS(prj.lcc)) #transform projection 

plot(subset(dat.pres, Species.Code=="ACMI"))
plot(subset(dat.pres, Species.Code=="ARUV"))
plot(subset(dat.pres, Species.Code=="VAME"))
