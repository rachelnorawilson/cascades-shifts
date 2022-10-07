# Created: Jan. 11, 2021
# Updated: Jan. 25, 2021

# IMPORTANT: Do NOT overwrite original rarefied dataset! This script is for posterity only.

# This script will be used to created a RAREFIED set of datasets, whittled to shared species (those appearing in both surveys). Rarefaction is used to correct for higher plot-level species richness at the resurvey time point, likely driven by  greater search and identification effort.

# Packages needed:

library(plyr)



#### STEP 1: Import data ####

und.cover <- read.csv("data/1_cover_with_fires.csv", header = TRUE, na.strings = "")
und.cover$Fires <- as.factor(und.cover$Fires)
und.cover$Data.Type <- as.factor(und.cover$Data.Type)
load("data/Species.List.Rda") #Note that this file was made in an undocumented step
species.list <- shifts$Species.Code[!shifts$Species.Code=="MOSS"] #removing "MOSS"
species.list <- factor(species.list)
load("data/plot.names.Rda")
load("data/list.fires.Rda")
plot.names$Plot.2015 <- factor(plot.names$Plot.2015)


# How many species were found per plot, on average, after reducing to common species?
mean(table(und.cover$Plot[und.cover$Data.Type == "Legacy"])) # 5.688347 species/plot
mean(table(und.cover$Plot[und.cover$Data.Type == "Resurvey"])) # 7.612466 species/plot

# In other words, the resurvey found ~2 more species per plot.




#### STEP 2: Create 100 rarefied presence datasets; store in list ####

rare.ALL <- list()

for(j in 1:100) {
  
  rare.resur.PLOTS.list <- list()
  
  # For resurvey data: on plot-by-plot basis remove 2 species, weighted by % cover
  for(i in 1:nrow(plot.names)) {
    und.cover.PLOT <- und.cover[und.cover$Plot == levels(plot.names$Plot.2015)[i], ]
    und.cover.PLOT$Species.Code <- factor(und.cover.PLOT$Species.Code)
    if(length(und.cover.PLOT$Species.Code) > 2) {
      Species.Code.Rare <- sample(x = und.cover.PLOT$Species.Code, 
                                  size = length(und.cover.PLOT$Species.Code)-2, 
                                  prob = und.cover.PLOT$Percent.Cover, 
                                  replace=FALSE)} else {
                                    Species.Code.Rare <- und.cover.PLOT$Species.Code}
    Species.Code.Rare <- factor(Species.Code.Rare)
    rare.resur.PLOTS.list[[i]] <- und.cover.PLOT[und.cover.PLOT$Species.Code 
                                                 %in% Species.Code.Rare, ]
  }
  
  rare.resur.PLOTS <- ldply(rare.resur.PLOTS.list, data.frame)
  leg.PLOTS <- und.cover[und.cover$Year=="1980", ]
  
  # Join rarefied resurvey dataset to legacy dataset
  und.rare <- rbind(rare.resur.PLOTS, leg.PLOTS)
  
  # Creating a new data frame from und.rare of 1s/0s presence/absence in each plot
  # (Same steps as in 1_tidydata.R)
  pres.abs.rare <- table(und.rare$Plot, und.rare$Species.Code)
  rare.presence.small <- melt(pres.abs.rare, id.vars=c("Plot", "Species.Code"))
  names(rare.presence.small) <- c("Plot", "Species.Code", "Pres.Abs")
  rare.pres <- merge(rare.presence.small, list.fires, by = "Plot")
  rare.pres$Data.Type <- as.factor(rare.pres$Data.Type)
  rare.pres$Pres.Abs <- ifelse(rare.pres$Pres.Abs >= 1, 1, 0)
  
  # Write presence/absence dataset to list
  
  rare.ALL[[j]] <- rare.pres
  if(j == 100) system("say Your loop is done")
} 

# Saving list of rarefied datasets as R object, to be imported during analysis step

# save(rare.ALL, file = "data/rare.ALL.Rda")








