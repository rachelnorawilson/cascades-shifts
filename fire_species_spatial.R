# Created: April 11, 2023
# Modified: 

#### This script plots fire-experiencing species on a map to visualize their distributions in the study area. Are they spatially correlated? Limited to the eastern slope? 


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
dat <- read_csv("data/1_presence_fires_unrarefied.csv")

# Plot locations
plots <- read_csv("data/0_Lat.Long.csv") %>% 
  drop_na() %>% 
  mutate(Longitude=ifelse(Longitude>0, -Longitude, Longitude)) #input file accidentally has some longitudes in E instead of W
colnames(plots)[2] = "Plot.Name.2015"
colnames(plots)[3] = "Plot.Name.1980"
plots$Plot.Name.1980 = as.character(plots$Plot.Name.1980)

# Plot fire history
fire.points <- read_csv("data/0_All_Plots_Wildfire_Join.csv") %>% 
  mutate(FireHistory = ifelse(CAL_YEAR>=1983, "Burned", "Unburned")) %>% 
  mutate(FireHistory = replace_na(FireHistory, "Unburned"))

# Merge locations and fire category with survey data
dat2 <- left_join(dat, plots, by=c("Plot"="Plot.Name.1980")) 
# NOTE: this join is adding several hundred rows - troubleshoot later if it affects relevant species

dat.all <- left_join(dat2, fire.points, by=c("Plot.Name.2015"="Name"))


### Set up maps
prj.wgs <- "+proj=longlat + type=crs"
prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs"

coordinates(dat.all) <- ~Longitude+Latitude #convert to spatial data
projection(dat.all) <- CRS('+proj=longlat') #define projection
dat.all <- spTransform(dat.all, CRS=CRS(prj.wgs)) #transform projection 
dat.all.lcc <- spTransform(dat.all, CRS=CRS(prj.lcc)) #transform projection 

