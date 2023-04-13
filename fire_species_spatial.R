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


### Set up maps

# Define projections
prj.wgs <- "+proj=longlat + type=crs"
prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs"

# Transform to spatial data
coordinates(dat.pres) <- ~Longitude+Latitude #convert to spatial data
projection(dat.pres) <- CRS('+proj=longlat') #define projection
dat.pres <- spTransform(dat.pres, CRS=CRS(prj.wgs)) #transform projection 
dat.pres.lcc <- spTransform(dat.pres, CRS=CRS(prj.lcc)) #transform projection 

# Define extent of study area
ext <- extent(min(dat.pres$Longitude)-0.5, max(dat.pres$Longitude)+0.5, min(dat.pres$Latitude)-0.5, max(dat.pres$Latitude)+0.5)
bbox <- as(ext, "SpatialPolygons") #convert coordinates to a bounding box

# USA polygons (used for lat/lon gridlines) 
sta = readOGR("data/shapefiles/states/gz_2010_us_040_00_500k.shp")
projection(sta) = CRS(prj.wgs)
# Crop state lines to study area
sta.crop <- crop(sta, bbox)

# Park boundary
park <- readOGR("data/shapefiles/park/NOCA_Park_boundary.shp")
park <- spTransform(park, CRS=CRS(prj.wgs))
park.lcc <- spTransform(park, CRS=CRS(prj.lcc))

# Wildfire polygons
fires <- readOGR("data/shapefiles/fires/NOCA_Wildfire_History.shp")
fires <- spTransform(fires, CRS=CRS(prj.wgs))
fires <- st_as_sf(fires) %>% 
  filter(CAL_YEAR>=1983)
fires.sp <- as(fires, "Spatial")
fires.lcc <- spTransform(fires.sp, CRS=CRS(prj.lcc))

# Prescribed burn polygons
burns <- readOGR("data/shapefiles/fires/Prescribed_burn_history.shp")
burns <- spTransform(burns, CRS=CRS(prj.wgs))
burns <- st_as_sf(burns) %>% 
  filter(CAL_YEAR>=1983)
burns.sp <- as(burns, "Spatial")
burns.lcc <- spTransform(burns.sp, CRS=CRS(prj.lcc))

# Fire treatment polygons (how are these different from prescribed burns?)
trtmts <- readOGR("data/shapefiles/fires/Fire_treatment_history.shp")
trtmts <- spTransform(trtmts, CRS=CRS(prj.wgs))
trtmts <- st_as_sf(trtmts) %>% 
  filter(TreatYear>=1983)
trtmts.sp <- as(trtmts, "Spatial")
trtmts.lcc <- spTransform(trtmts.sp, CRS=CRS(prj.lcc))

# Set up gridlines & lat/lon labels	
frame.grd <- gridlines(sta.crop)
frame.grd.lcc <- spTransform(frame.grd, CRS=CRS(prj.lcc))
gridatt <- gridat(frame.grd, side="EN")
gridat.lcc <- spTransform(gridatt, CRS=CRS(prj.lcc))


### Maps of fire-experiencing species

# ACMI
#pdf(file="figures/map_ACMI_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="ACMI" & Data.Type=="Legacy"), pch=1, cex=1.5, col="black", add=T) # historical presences of ACMI
plot(subset(dat.pres.lcc, Species.Code=="ACMI"& Data.Type=="Resurvey"), pch=4, cex=1.5, col="black", add=T) # contemporary presences of ACMI
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
#dev.off()

# ARUV
#pdf(file="figures/map_VAME_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="ARUV" & Data.Type=="Legacy"), pch=1, cex=1.5, col="black", add=T) # historical presences of ARUV
plot(subset(dat.pres.lcc, Species.Code=="ARUV"& Data.Type=="Resurvey"), pch=4, cex=1.5, col="black", add=T) # contemporary presences of ARUV
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
#dev.off()

# VAME
#pdf(file="figures/map_VAME_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="VAME" & Data.Type=="Legacy"), pch=1, cex=1.5, col="black", add=T) # historical presences of VAME
plot(subset(dat.pres.lcc, Species.Code=="VAME"& Data.Type=="Resurvey"), pch=4, cex=1.5, col="black", add=T) # contemporary presences of VAME
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
#dev.off()

# CARU
#pdf(file="figures/map_CARU_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="CARU" & Data.Type=="Legacy"), pch=1, cex=1.5, col="black", add=T) # historical presences of CARU
plot(subset(dat.pres.lcc, Species.Code=="CARU"& Data.Type=="Resurvey"), pch=4, cex=1.5, col="black", add=T) # contemporary presences of CARU
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
#dev.off()

# CEVE
#pdf(file="figures/map_CEVE_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="CEVE" & Data.Type=="Legacy"), pch=1, cex=1.5, col="black", add=T) # historical presences of CEVE
plot(subset(dat.pres.lcc, Species.Code=="CEVE"& Data.Type=="Resurvey"), pch=4, cex=1.5, col="black", add=T) # contemporary presences of CEVE
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
#dev.off()

# EPAN (CHAN)
#pdf(file="figures/map_EPAN_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="EPAN" & Data.Type=="Legacy"), pch=1, cex=1.5, col="black", add=T) # historical presences of EPAN
plot(subset(dat.pres.lcc, Species.Code=="EPAN"& Data.Type=="Resurvey"), pch=4, cex=1.5, col="black", add=T) # contemporary presences of EPAN
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
#dev.off()

# PAMY
#pdf(file="figures/map_PAMY_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="PAMY" & Data.Type=="Legacy"), pch=1, cex=1.5, col="black", add=T) # historical presences of PAMY
plot(subset(dat.pres.lcc, Species.Code=="PAMY"& Data.Type=="Resurvey"), pch=4, cex=1.5, col="black", add=T) # contemporary presences of PAMY
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
#dev.off()

