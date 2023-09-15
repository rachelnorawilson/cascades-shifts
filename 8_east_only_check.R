# Created: June 3, 2023

# This script will NOT produce the final version of the analysis.
# It is solely for determining whether subsetting to east-side plots affects
# fire-experiencing species patterns. It includes versions of:
# - 3b* (to generate error logs),
# - 3c (to run models and create coefficients CSV used in vis),
# - 6 (to visualize model predictions)
# - 5 - (to visualize arithmetic changes)
## *this script differs from 3b in that it forces all species to conform to fire models.

# Part of this script (3b) creates a large file and places it in a folder
# "warning_logs"; these are not stored on GitHub but can be shared upon request and/or
# reproduced locally. User will need to adjust the destination directory (lines 102 and 104).

# No CIs are generated in this script.

# IMPORTANT NOTE: unless otherwise indicated, always use Understory_All.csv for these analyses as it is the ONLY file with up-to-date corrections.

################# SCRIPT 3B ANALOGUE ###############

###### NOTE: GitHub will die if you try to push 700 error logs. Just delete them locally - they aren't needed after this section is completed (i.e., before running script 3C)

# Packages needed:

library(MuMIn)
library(dplyr)
source("3a_dredge_log_to_df.R")

# Function to create data frame of coefficients, AIC and model ID:
df.fun <- function(ModID) {
  df <- as.data.frame(t(coef(dredge.list[[ModID]]))) # Transposed DF of coefs
  df$new_Model_id <- paste(ModID) # Store unique model ID
  df$AIC <- dredge.list[[ModID]]$aic # Store model AIC
  df$logLik <- as.numeric(logLik(dredge.list[[ModID]])) # Store logLik
  return(df)
}

# List of coefficients between fire and non-fire model framework
coeff.all <-c( "Elevation.m", 
               "Elevation.m2", 
               "New.Data.TypeResurvey.Burned", 
               "New.Data.TypeResurvey.Unburned", 
               "Elevation.m:New.Data.TypeResurvey.Burned", 
               "Elevation.m:New.Data.TypeResurvey.Unburned", 
               "Elevation.m2:New.Data.TypeResurvey.Burned", 
               "Elevation.m2:New.Data.TypeResurvey.Unburned",
               "Data.TypeResurvey",
               "Data.TypeResurvey:Elevation.m",
               "Data.TypeResurvey:Elevation.m2")

#### STEP 1: Import data #### Takes ~ 5 sec to run

# To analyze un-rarefied data (exclude 100x loop):
# und.presence <- read.csv("data/1_presence_fires_unrarefied", header = TRUE, na.strings = "")

# To analyze rarefied data (as list of dataframes named rare.ALL):
load("data/rare.ALL.Rda") # Run time 5 sec

# Import for both data types:
load("data/0.Species.List.Rda")

# Load list of plot regions
regions <- read.csv("data/plot_regions.csv")
east.regions <- regions[regions$NAME == "East Cascades", 3:4]
east.list <- c(east.regions$Plot.Name.1980, east.regions$Plot.Name.2015)

# Modify species list to include ONLY fire-experiencing species (overwrite file name)
species.list <- c("ACMI", "ARUV", "CARU", "CEVE", "EPAN", "PAMY", "VAME")
species.list <- factor(species.list)

#TODO Formula to correct erroneous burn coding. Run time 3 sec
for(D in 1:100) {
  rare.ALL[[D]]$Fires[rare.ALL[[D]]$Data.Type == "Legacy"] <- rep("Unburned")
  rare.ALL[[D]]$New.Data.Type <- paste(rare.ALL[[D]]$Data.Type, rare.ALL[[D]]$Fires, sep = ".")
}


#### STEP 2: Loop to analyze presence data - FIRE SPECIES ONLY 

# Exclude D loop to run for only one dataset (e.g. unrarefied data)

coeff.ALLDAT <- list()
warn.ALLDAT <- list()

for(D in 1:100) { #RUN TIME: ~5 min
  
  coeff.ALLSPEC <- list()
  warn.ALLSPEC <- list()
  und.presence.init <- rare.ALL[[D]]
  und.presence <- und.presence.init[und.presence.init$Plot %in% east.list, ]
  und.presence$Plot <- factor(und.presence$Plot)
  und.presence$Elevation.m <- as.numeric(und.presence$Elevation.m)
  und.presence$New.Data.Type <- factor(und.presence$New)
  und.presence$Fires <- as.factor(und.presence$Fires)
  und.presence$Data.Type <- as.factor(und.presence$Data.Type)
  
  for(S in 1:length(species.list)) {
    
    # Create subset for species of interest S
    und.presence.SPEC <- subset(und.presence, Species.Code == levels(species.list)[S])
    und.presence.SPEC$New.Data.Type <- relevel(und.presence.SPEC$New.Data.Type, 
                                               ref="Legacy.Unburned")
    und.presence.SPEC$Fires <- relevel(und.presence.SPEC$Fires, ref="Unburned")
    und.presence.SPEC$Data.Type <- relevel(und.presence.SPEC$Data.Type, ref="Legacy")
    und.presence.SPEC$Elevation.m <- poly(und.presence.SPEC$Elevation.m, 2)[ , 1] #different
    und.presence.SPEC$Elevation.m2 <- poly(und.presence.SPEC$Elevation.m, 2)[ , 2] #different
    und.presence.SPEC <- und.presence.SPEC[complete.cases(und.presence.SPEC), ] #Just in case
    
    # Emptying out previous model objects
    mod.globfi <- NULL
    mod.globnofi <- NULL
    
    # Did the species occur in burned plots 5+ times?
    num.burns <- 
      table(und.presence.SPEC$Pres.Abs, und.presence.SPEC$Fires, und.presence.SPEC$Data.Type)
    
    
    # Create new error-logging file specific to permutation & species
    log.file <- file(paste("data/temp_warnings/warnings", 
                           D, levels(species.list)[S], ".txt", sep = "_"), open = "wt")
    log.file.path <- paste("data/temp_warnings/warnings", 
                           D, levels(species.list)[S], ".txt", sep = "_")
    sink(log.file, append = TRUE, type = "output") # Sink to log file
    sink(log.file, append = TRUE, type = "message")
    
    
    # If burn > 5, include fire as a predictor in global model (New.Data.Type):

      mod.globfi <- glm(Pres.Abs ~ (Elevation.m + Elevation.m2) * New.Data.Type, 
                        data = und.presence.SPEC, family = "binomial", na.action = na.fail)
      options(warn = -1) # Ignore warnings - not for logging.
      dredge.list <- lapply(dredge(mod.globfi, rank = AIC, subset = 
                                     dc(Elevation.m, Elevation.m2) &&
                                     dc(New.Data.Type:Elevation.m, New.Data.Type:Elevation.m2) &&
                                     dc(Elevation.m:New.Data.Type, Elevation.m2:New.Data.Type), 
                                   trace = FALSE, evaluate = FALSE), eval)
      names(dredge.list) <- paste("Mod", 
                                  as.numeric(names(dredge.list)) - 1, 
                                  sep = ".") # Converting to model ID
      options(warn = 1) # Tell me if a model throws an error - for logging.
      dredge.globfi <- dredge(mod.globfi, rank = AIC, subset = 
                                dc(Elevation.m, Elevation.m2) &&
                                dc(New.Data.Type:Elevation.m, New.Data.Type:Elevation.m2) &&
                                dc(Elevation.m:New.Data.Type, Elevation.m2:New.Data.Type), 
                              trace = 1)
    
    
    # Stop writing to .txt, import dataframe of warning logs
    closeAllConnections()
    log.warn <- dredge_log_to_df(log.file.path) # See 3a_dredge_log_to_df.R
    log.warn$new_Model_id <- paste("Mod", log.warn$Model_id, sep = ".")
    
    # Run df.fun to pull out coefficients, AIC, model ID from dredge list
    coeff.df <- ldply(lapply(log.warn$new_Model_id, df.fun))
    
    # Join warning log to dataframe of coefficients and AIC
    coeff.warn <- merge(coeff.df, log.warn, by = "new_Model_id", all.x = TRUE)
    
    # Exclude models for which there was a warning
    coeff.nowarn <- coeff.warn[coeff.warn$Has_warning == FALSE, ]
    
    # Calculate delta AIC based on warning-less models, reduce to delta <=2, add weights
    coeff.nowarn$delta <- coeff.nowarn$AIC - min(coeff.nowarn$AIC)
    top.mods.coeff <- coeff.nowarn[coeff.nowarn$delta <= 2, ]
    top.mods.coeff$weight <- Weights(top.mods.coeff$AIC)
    
    # Null model (for Psuedo-R-squared calculation later)
    mod.NULL <- glm(Pres.Abs ~ 1, 
                    data = und.presence.SPEC, family = "binomial", na.action = na.fail)
    
    # Adding in missing coefficient headers
    for(C in 1:length(coeff.all)) {
      if(!coeff.all[C] %in% colnames(top.mods.coeff)) {
        top.mods.coeff[, coeff.all[C]] <- rep(NA, times = nrow(top.mods.coeff))
      }
    }
    
    ## Storing top model coefficients
    
    Mods.list <- list()
    
    for(i in 1:nrow(top.mods.coeff)) {
      Mods.list[[i]] <- data.frame(
        Species = levels(species.list)[S], 
        Dataset = D,
        L.Occ = sum(num.burns["1", , "Legacy"]), 
        R.Occ = sum(num.burns["1", , "Resurvey"]), 
        Fire.Included = ifelse(is.null(mod.globnofi) == TRUE, "Yes", "No"),
        Type = "Unavg", 
        Mod.ID = top.mods.coeff$new_Model_id[i],
        deltaAIC = top.mods.coeff$delta[i], 
        Weight = top.mods.coeff$weight[i], 
        Rsquared = 1 - top.mods.coeff$logLik[i] / as.numeric(logLik(mod.NULL)),
        Intercept = top.mods.coeff$`(Intercept)`[i],
        Elevation.m = top.mods.coeff$Elevation.m[i],
        Elevation.m2 = top.mods.coeff$Elevation.m2[i],
        Data.Type.nofi = top.mods.coeff$Data.TypeResurvey[i],
        Data.Type.Elevation.m.nofi = top.mods.coeff$`Data.TypeResurvey:Elevation.m`[i],
        Data.Type.Elevation.m2.nofi = top.mods.coeff$`Data.TypeResurvey:Elevation.m2`[i],
        Resurvey.Burned.fi = 
          top.mods.coeff$New.Data.TypeResurvey.Burned[i],
        Resurvey.Unburned.fi = 
          top.mods.coeff$New.Data.TypeResurvey.Unburned[i],
        Elevation.m.Res.Burn.fi = 
          top.mods.coeff$`Elevation.m:New.Data.TypeResurvey.Burned`[i],
        Elevation.m.Res.Unburn.fi = 
          top.mods.coeff$`Elevation.m:New.Data.TypeResurvey.Unburned`[i],
        Elevation.m2.Res.Burn.fi = 
          top.mods.coeff$`Elevation.m2:New.Data.TypeResurvey.Burned`[i],
        Elevation.m2.Res.Unburn.fi = 
          top.mods.coeff$`Elevation.m2:New.Data.TypeResurvey.Unburned`[i],
        row.names = NULL)
    }
    
    Mods <- ldply(Mods.list, data.frame)
    
    ## Storing model-averaged parameters
    
    # Re-code NA coefficients as 0
    top.mods.zeroes <- top.mods.coeff
    top.mods.zeroes[is.na(top.mods.zeroes)] <- 0
    
    Avg <- data.frame(
      Species = levels(species.list)[S],
      Dataset = D,
      L.Occ = sum(num.burns["1", , "Legacy"]), 
      R.Occ = sum(num.burns["1", , "Resurvey"]), 
      Fire.Included = ifelse(is.null(mod.globnofi) == TRUE, "Yes", "No"),
      Type = "Avg", 
      Mod.ID = NA,
      deltaAIC = NA, 
      Weight = NA,
      Rsquared = NA,
      Intercept = sum(top.mods.zeroes$`(Intercept)` * top.mods.zeroes$weight),
      Elevation.m = sum(top.mods.zeroes$`Elevation.m` * top.mods.zeroes$weight),
      Elevation.m2 = sum(top.mods.zeroes$`Elevation.m2` * top.mods.zeroes$weight),
      Data.Type.nofi = sum(top.mods.zeroes$`Data.TypeResurvey` * top.mods.zeroes$weight),
      Data.Type.Elevation.m.nofi = sum(top.mods.zeroes$`Data.TypeResurvey:Elevation.m` * 
                                         top.mods.zeroes$weight),
      Data.Type.Elevation.m2.nofi = sum(top.mods.zeroes$`Data.TypeResurvey:Elevation.m2` * 
                                          top.mods.zeroes$weight),
      Resurvey.Burned.fi = 
        sum(top.mods.zeroes$`New.Data.TypeResurvey.Burned` * top.mods.zeroes$weight),
      Resurvey.Unburned.fi = 
        sum(top.mods.zeroes$`New.Data.TypeResurvey.Unburned` * top.mods.zeroes$weight),
      Elevation.m.Res.Burn.fi = 
        sum(top.mods.zeroes$`Elevation.m:New.Data.TypeResurvey.Burned` * top.mods.zeroes$weight),
      Elevation.m.Res.Unburn.fi = 
        sum(top.mods.zeroes$`Elevation.m:New.Data.TypeResurvey.Unburned` * top.mods.zeroes$weight),
      Elevation.m2.Res.Burn.fi = 
        sum(top.mods.zeroes$`Elevation.m2:New.Data.TypeResurvey.Burned` * top.mods.zeroes$weight),
      Elevation.m2.Res.Unburn.fi = 
        sum(top.mods.zeroes$`Elevation.m2:New.Data.TypeResurvey.Unburned` * 
              top.mods.zeroes$weight),
      row.names = NULL)
    
    
    coeff.ALLSPEC[[S]] <- rbind(Mods, Avg)
    
    ## Storing errors
    
    coeff.warn$Species <- levels(species.list)[S]
    coeff.warn$Dataset <- D
    coeff.warn$L.Occ <- sum(num.burns["1", , "Legacy"])
    coeff.warn$R.Occ <- sum(num.burns["1", , "Resurvey"])
    coeff.warn$Fire.Included <- ifelse(is.null(mod.globnofi) == TRUE, "Yes", "No")
    coeff.warn$Type <- "Unavg.w.warn"
    
    warn.ALLSPEC[[S]] <- coeff.warn
    
  }  
  
  #Collapse ALLSPEC list into single dataframe, store that df as part of ALLDAT list
  
  coeff.ALLDAT[[D]] <- ldply(coeff.ALLSPEC, data.frame)
  warn.ALLDAT[[D]] <- ldply(warn.ALLSPEC, data.frame)
  
  #### END OF SPECIES LOOP
  
  if(D == 100) system("say Your loop is done")
  
}

#### END OF DATASET LOOP

# Collate ALLDAT lists into one big DF

coeff.ALLDAT.finaldf <- ldply(coeff.ALLDAT, data.frame)
warn.ALLDAT.finaldf <- ldply(warn.ALLDAT, data.frame)


# Store output as CSV.

write.csv(coeff.ALLDAT.finaldf, 
          file = "data/east_only_check/eastonly_3_presence_ALLDAT_ALLSPEC_coefficients.csv", 
          row.names = FALSE) 
write.csv(warn.ALLDAT.finaldf, 
          file = "data/east_only_check/eastonly_3_presence_ALLDAT_ALLSPEC_warnings.csv", 
          row.names = FALSE)


###################### SCRIPT 3C ANALOGUE #########################

# This script will be used to undertake part c of the PRESENCE analyses (modeling)
# Use this script to produce the FINAL version of the analysis, plus CIs. No error logging.

# This script accommodates the following species-specific issues:
# --> Exclude the most-complex model (X * elev^2) for species in which this model threw warnings
# # ---> VAME
# --> For species that flip-flopped between yes/no fire, use "majority rules" to decide
#     which framework to use
# # ---> Fire: ARUV, CARU, VAME*
# As well as the following general changes:
# --> Switch elevation^2 term to poly(elevation)

# Packages needed:

library(MuMIn)
library(plyr)


# List of coefficients between fire and non-fire model framework
coeff.all <-c("(Intercept)",
              "Elevation.m.poly", 
              "Elevation.m2.poly", 
              "New.Data.TypeResurvey.Burned", 
              "New.Data.TypeResurvey.Unburned", 
              "Elevation.m.poly:New.Data.TypeResurvey.Burned", 
              "Elevation.m.poly:New.Data.TypeResurvey.Unburned", 
              "Elevation.m2.poly:New.Data.TypeResurvey.Burned", 
              "Elevation.m2.poly:New.Data.TypeResurvey.Unburned",
              "Data.TypeResurvey",
              "Data.TypeResurvey:Elevation.m.poly",
              "Data.TypeResurvey:Elevation.m2.poly")

#### STEP 1: Import data #### Takes ~ 5 sec to run

# To analyze un-rarefied data (exclude 100x loop):
# und.presence <- read.csv("data/1_presence_fires_unrarefied", header = TRUE, na.strings = "")

# To analyze rarefied data (as list of dataframes named rare.ALL):
load("data/rare.ALL.Rda") # Run time 5 sec

# Import for both data types:
load("data/Species.List.Rda") #TODO this file was made in an undocumented step
# Modify species list to include ONLY fire-experiencing species (overwrite file name)
species.list <- c("ACMI", "ARUV", "CARU", "CEVE", "EPAN", "PAMY", "VAME")
species.list <- factor(species.list)

# Load list of plot regions
regions <- read.csv("data/plot_regions.csv")
east.regions <- regions[regions$NAME == "East Cascades", 3:4]
east.list <- c(east.regions$Plot.Name.1980, east.regions$Plot.Name.2015)

# Porting in CSV of warnings
warn.ALLDAT <- read.csv("data/east_only_check/eastonly_3_presence_ALLDAT_ALLSPEC_warnings.csv", header = TRUE)
# Bug correction - force VAME and HODI to run under correct framework
warn.ALLDAT$Has_warning[warn.ALLDAT$Species == "VAME" 
                        & warn.ALLDAT$Dataset == 36] <- paste("TRUE")

# Formula to correct erroneous burn coding. Run time 3 sec
for(D in 1:100) {
  rare.ALL[[D]]$Fires[rare.ALL[[D]]$Data.Type == "Legacy"] <- rep("Unburned")
  rare.ALL[[D]]$New.Data.Type <- paste(rare.ALL[[D]]$Data.Type, rare.ALL[[D]]$Fires, sep = ".")
}


#### STEP 2: Loop to analyze presence data

# Can be run as a loop outputting all species, or S can be modified to isolated specific species. Check number here:
(numbered.species <- data.frame(Species=species.list, No.=rep(1:length(species.list))))
species.with.fire <- numbered.species[numbered.species$Species == "ACMI" |
                                        numbered.species$Species == "ARUV" |
                                        numbered.species$Species == "CARU" |
                                        numbered.species$Species == "CEVE" |
                                        numbered.species$Species == "EPAN" |
                                        numbered.species$Species == "PAMY" |
                                        numbered.species$Species == "VAME", ]
species.without.fire <- numbered.species[!numbered.species$Species 
                                         %in% species.with.fire$Species, ]
# Exclude D loop to run for only one dataset (e.g. unrarefied data)

coeff.ALLDAT <- list() # Store coefficient outputs
avg.confint.ALLDAT <- list() # Store model-averaged coefficients
framework.ALLDAT <- list () # Store record of which decision framework was used

for(D in 1:100) { #RUN TIME: 5 min
  
  coeff.ALLSPEC <- list()
  avg.confint.ALLSPEC <- list() 
  framework.ALLSPEC <- list()
  und.presence.init <- rare.ALL[[D]]
  und.presence <- und.presence.init[und.presence.init$Plot %in% east.list, ]
  und.presence$Plot <- factor(und.presence$Plot)
  und.presence$Elevation.m <- as.numeric(und.presence$Elevation.m)
  und.presence$New.Data.Type <- factor(und.presence$New.Data.Type)
  und.presence$Fires <- as.factor(und.presence$Fires)
  und.presence$Data.Type <- as.factor(und.presence$Data.Type)
  
  # Create dataframe of warnings for this dataset
  warn.dataset <- subset(warn.ALLDAT, Dataset == D,
                         select = c(Species, Dataset, Fire.Included, Has_warning))
  
  for(S in 1:length(species.list)) {
    
    # Create subset for species of interest S
    und.presence.SPEC <- subset(und.presence, Species.Code == levels(species.list)[S])
    und.presence.SPEC$New.Data.Type <- relevel(und.presence.SPEC$New.Data.Type, 
                                               ref="Legacy.Unburned")
    und.presence.SPEC$Fires <- relevel(und.presence.SPEC$Fires, ref="Unburned")
    und.presence.SPEC$Data.Type <- relevel(und.presence.SPEC$Data.Type, ref="Legacy")
    #und.presence.SPEC$Elevation.m2 <- und.presence.SPEC$Elevation.m^2
    # Added Mar. 1: Replace Elevation.m with poly()
    und.presence.SPEC$Elevation.m.poly <- poly(und.presence.SPEC$Elevation.m, 2)[ , 1]
    und.presence.SPEC$Elevation.m2.poly <- poly(und.presence.SPEC$Elevation.m, 2)[ , 2]
    und.presence.SPEC <- und.presence.SPEC[complete.cases(und.presence.SPEC), ] #Just in case
    
    # Create subset of warnings for species of interest S
    warn.SPEC <- subset(warn.dataset, Species == levels(species.list)[S])
    
    # Create empty data frame to record decision framework
    framework.SPEC <- data.frame(Dataset = paste(D), 
                                 Species = paste(levels(species.list)[S]),
                                 Forced.Fire = paste(NA),
                                 Forced.No.Fire = paste(NA),
                                 Forced.Simpler.Mod = paste(NA),
                                 Discard.Later = paste(NA),
                                 One.Top.Mod = paste(NA))
    
    # Emptying out previous model objects
    mod.globfi <- NULL
    mod.globnofi <- NULL
    mod.globfi.reduced <- NULL
    mod.globnofi.reduced <- NULL
    dredge.globfi <- NULL
    dredge.globnofi <- NULL
    avg.mods.coeff <- NULL
    avg.mods.confint <- NULL
    top.mods.coeff <- NULL
    
    # Did the species occur in burned plots 5+ times?
    num.burns <- 
      table(und.presence.SPEC$Pres.Abs, und.presence.SPEC$Fires, und.presence.SPEC$Data.Type)
    
    #### START OF MODEL FRAMEWORK #### 
    
    # If burn > 5 in >50% of datasets, include fire as a predictor in global model (New.Data.Type):
    if(levels(species.list)[S] %in% species.with.fire$Species == TRUE) {
      
      # Did we force this D & S to run in the fire framework?
      if(levels(factor(warn.SPEC$Fire.Included)) == "No") {
        framework.SPEC$Forced.Fire <- paste("Yes") # Record
      } else {
        framework.SPEC$Forced.Fire <- paste("No") # Record
      }
      
      # Was there a warning associated with this D & S?
      if("TRUE" %in% levels(factor(warn.dataset$Has_warning[warn.dataset$Species 
                                                            == levels(species.list)[S]]))) {
        
        if(levels(species.list)[S] == "VAME") { # No elev^2 * year-burn
          mod.globfi.reduced <- glm(Pres.Abs ~ (Elevation.m.poly + Elevation.m2.poly) + 
                                      New.Data.Type + Elevation.m.poly:New.Data.Type, 
                                    data = und.presence.SPEC, family = "binomial", na.action = na.fail)
          dredge.globfi.reduced <- dredge(mod.globfi.reduced, rank = AIC, subset = 
                                            dc(Elevation.m.poly, Elevation.m2.poly))
          
          # Model averaging:
          if(nrow(subset(dredge.globfi.reduced, delta <= 2)) == 1) { # Error workaround
            avg.mods.coeff <- as.data.frame(coef(subset(dredge.globfi.reduced, delta <= 2)))
            avg.mods.confint <- #TODO Hideously clunky
              as.data.frame(t(as.data.frame(lapply(get.models(dredge.globfi.reduced, 
                                                              subset = delta <= 2), confint))))
            row.names(avg.mods.confint) <- c("2.5 %", "97.5 %")
            framework.SPEC$One.Top.Mod <- paste("Yes")
          } else { # Normal
            avg.mods <- model.avg(dredge.globfi.reduced, subset = delta <= 2)
            avg.mods.coeff <- as.data.frame(t(avg.mods$coefficients["full",]))
            avg.mods.confint <- as.data.frame(t(confint(avg.mods, full = TRUE)))
            framework.SPEC$One.Top.Mod <- paste("No")
          }
          
          # Coefficients:
          top.mods.coeff <- as.data.frame(coef(subset(dredge.globfi.reduced, delta <= 2)))
          top.mods.coeff$delta <- dredge.globfi.reduced$delta[dredge.globfi.reduced$delta <= 2]
          top.mods.coeff$weight <- Weights(top.mods.coeff$delta)
          top.mods.coeff$logLik <- dredge.globfi.reduced$logLik[dredge.globfi.reduced$delta <= 2]
          
          framework.SPEC$Forced.Simpler.Mod <- paste("Yes") # Record
          
        } else {
          framework.SPEC$Discard.Later <- paste("Yes") # Record
          # Create empty DFs
          top.mods.coeff <- data.frame(delta = rep(NA, 1), 
                                       weight = rep(NA, 1),
                                       logLik = rep(NA, 1))
          avg.mods.coeff <- data.frame(Elevation.m.poly = rep(NA, 1))
          avg.mods.confint <- data.frame(Elevation.m.poly = rep(NA, 2), 
                                         row.names = c("2.5 %", "97.5 %"))
        }
        
        
      } else { # Run as normal
        mod.globfi <- glm(Pres.Abs ~ (Elevation.m.poly + Elevation.m2.poly) * New.Data.Type, 
                          data = und.presence.SPEC, family = "binomial", na.action = na.fail)
        dredge.globfi <- dredge(mod.globfi, rank = AIC, subset = 
                                  dc(Elevation.m.poly, Elevation.m2.poly) &&
                                  dc(New.Data.Type:Elevation.m.poly, 
                                     New.Data.Type:Elevation.m2.poly) &&
                                  dc(Elevation.m.poly:New.Data.Type, 
                                     Elevation.m2.poly:New.Data.Type))
        
        # Model averaging:
        if(nrow(subset(dredge.globfi, delta <= 2)) == 1) { # Error workaround
          avg.mods.coeff <- as.data.frame(coef(subset(dredge.globfi, delta <= 2)))
          avg.mods.confint <- #TODO Hideously clunky
            as.data.frame(t(as.data.frame(lapply(get.models(dredge.globfi, 
                                                            subset = delta <= 2), confint))))
          row.names(avg.mods.confint) <- c("2.5 %", "97.5 %")
          framework.SPEC$One.Top.Mod <- paste("Yes")
        } else { # Normal
          avg.mods <- model.avg(dredge.globfi, subset = delta <= 2)
          avg.mods.coeff <- as.data.frame(t(avg.mods$coefficients["full",]))
          avg.mods.confint <- as.data.frame(t(confint(avg.mods, full = TRUE)))
          framework.SPEC$One.Top.Mod <- paste("No")
        }
        
        # Storing coefficients:
        top.mods.coeff <- as.data.frame(coef(subset(dredge.globfi, delta <= 2)))
        top.mods.coeff$delta <- dredge.globfi$delta[dredge.globfi$delta <= 2]
        top.mods.coeff$weight <- Weights(top.mods.coeff$delta)
        top.mods.coeff$logLik <- dredge.globfi$logLik[dredge.globfi$delta <= 2]
        
      }
    }
    
    #If burn < 5 in >50% of datasets, exclude fire from global model:
    if(levels(species.list)[S] %in% species.with.fire$Species == FALSE) {
      
      # Did we force this D & S to run in the no-fire framework?
      if(levels(factor(warn.SPEC$Fire.Included)) == "Yes") {
        framework.SPEC$Forced.No.Fire <- paste("Yes") # Record
      } else {
        framework.SPEC$Forced.No.Fire <- paste("No") # Record
      }
      
      # Was there a warning associated with this D & S?
      if("TRUE" %in% levels(factor(warn.dataset$Has_warning[warn.dataset$Species 
                                                            == levels(species.list)[S]]))) {
        if(levels(species.list)[S] == "HODI") { # No elev^2 * year
          mod.globnofi.reduced <- glm(Pres.Abs ~ Elevation.m.poly + Elevation.m2.poly + 
                                        Data.Type + Elevation.m.poly:Data.Type, 
                                      data = und.presence.SPEC, family = "binomial", 
                                      na.action = na.fail)
          dredge.globnofi.reduced <- dredge(mod.globnofi.reduced, rank = AIC, subset = 
                                              dc(Elevation.m.poly, Elevation.m2.poly))
          
          # Model averaging:
          if(nrow(subset(dredge.globnofi.reduced, delta <= 2)) == 1) { # Error workaround
            avg.mods.coeff <- as.data.frame(coef(subset(dredge.globnofi.reduced, delta <= 2)))
            avg.mods.confint <- #TODO Hideously clunky
              as.data.frame(t(as.data.frame(lapply(get.models(dredge.globnofi.reduced, 
                                                              subset = delta <= 2), confint))))
            row.names(avg.mods.confint) <- c("2.5 %", "97.5 %")
            framework.SPEC$One.Top.Mod <- paste("Yes")
          } else { # Normal
            avg.mods <- model.avg(dredge.globnofi.reduced, subset = delta <= 2)
            avg.mods.coeff <- as.data.frame(t(avg.mods$coefficients["full",]))
            avg.mods.confint <- as.data.frame(t(confint(avg.mods, full = TRUE)))
            framework.SPEC$One.Top.Mod <- paste("No")
          }
          
          # Coefficients:
          top.mods.coeff <- as.data.frame(coef(subset(dredge.globnofi.reduced, delta <= 2)))
          top.mods.coeff$delta <- dredge.globnofi.reduced$delta[dredge.globnofi.reduced$delta <= 2]
          top.mods.coeff$weight <- Weights(top.mods.coeff$delta)
          top.mods.coeff$logLik <- 
            dredge.globnofi.reduced$logLik[dredge.globnofi.reduced$delta <= 2]
          
          framework.SPEC$Forced.Simpler.Mod <- paste("Yes") # Record
          
        } else {
          framework.SPEC$Discard.Later <- paste("Yes") # Record
          # Create empty DFs
          top.mods.coeff <- data.frame(delta = rep(NA, 1), 
                                       weight = rep(NA, 1), 
                                       logLik = rep(NA, 1))
          avg.mods.coeff <- data.frame(Elevation.m.poly = rep(NA, 1))
          avg.mods.confint <- data.frame(Elevation.m.poly = rep(NA, 2), 
                                         row.names = c("2.5 %", "97.5 %"))
        }
        
      } else { # Run as normal
        mod.globnofi <- glm(Pres.Abs ~ Data.Type * (Elevation.m.poly + Elevation.m2.poly), 
                            data = und.presence.SPEC, family = "binomial", na.action = na.fail) 
        dredge.globnofi <- dredge(mod.globnofi, rank = AIC, subset = 
                                    dc(Elevation.m.poly, Elevation.m2.poly) &&
                                    dc(Data.Type:Elevation.m.poly, Data.Type:Elevation.m2.poly))
        
        # Model averaging
        if(nrow(subset(dredge.globnofi, delta <= 2)) == 1) { # Error workaround
          avg.mods.coeff <- as.data.frame(coef(subset(dredge.globnofi, delta <= 2)))
          avg.mods.confint <- #TODO Hideously clunky
            as.data.frame(t(as.data.frame(lapply(get.models(dredge.globnofi, 
                                                            subset = delta <= 2), confint))))
          row.names(avg.mods.confint) <- c("2.5 %", "97.5 %")
          framework.SPEC$One.Top.Mod <- paste("Yes")
        } else { # Normal
          avg.mods <- model.avg(dredge.globnofi, subset = delta <= 2)
          avg.mods.coeff <- as.data.frame(t(avg.mods$coefficients["full",]))
          avg.mods.confint <- as.data.frame(t(confint(avg.mods, full = TRUE)))
          framework.SPEC$One.Top.Mod <- paste("No")
        }
        
        # Coefficients:
        top.mods.coeff <- as.data.frame(coef(subset(dredge.globnofi, delta <= 2)))
        top.mods.coeff$delta <- dredge.globnofi$delta[dredge.globnofi$delta <= 2]
        top.mods.coeff$weight <- Weights(top.mods.coeff$delta)
        top.mods.coeff$logLik <- dredge.globnofi$logLik[dredge.globnofi$delta <= 2]
      }
    }
    
    
    
    #### END OF MODEL FRAMEWORK ####
    # Output: 3 DFs named top.mods.coeff, avg.mods.coeff, avg.mods.confint
    
    
    # Null model (for Psuedo-R-squared calculation later)
    mod.NULL <- glm(Pres.Abs ~ 1, 
                    data = und.presence.SPEC, family = "binomial", na.action = na.fail)
    
    # Adding in missing coefficient headers
    for(C in 1:length(coeff.all)) {
      if(!coeff.all[C] %in% colnames(top.mods.coeff)) {
        top.mods.coeff[, coeff.all[C]] <- rep(NA, times = nrow(top.mods.coeff))
      }
      if(!coeff.all[C] %in% colnames(avg.mods.coeff)) {
        avg.mods.coeff[, coeff.all[C]] <- rep(NA, times = nrow(avg.mods.coeff))
      }
      if(!coeff.all[C] %in% colnames(avg.mods.confint)) {
        avg.mods.confint[, coeff.all[C]] <- rep(NA, times = nrow(avg.mods.confint))
      }
    }
    
    ## Storing top model coefficients
    
    Mods.list <- list()
    
    for(i in 1:nrow(top.mods.coeff)) {
      Mods.list[[i]] <- data.frame(
        Species = levels(species.list)[S], 
        Dataset = D,
        L.Occ = sum(num.burns["1", , "Legacy"]), 
        R.Occ = sum(num.burns["1", , "Resurvey"]), 
        Fire.Included = ifelse(levels(species.list)[S] 
                               %in% species.with.fire$Species == TRUE, "Yes", "No"),
        Type = "Unavg", 
        deltaAIC = top.mods.coeff$delta[i], 
        Weight = top.mods.coeff$weight[i], 
        Rsquared = 1 - top.mods.coeff$logLik[i] / as.numeric(logLik(mod.NULL)),
        Intercept = top.mods.coeff$`(Intercept)`[i],
        Elevation.m = top.mods.coeff$Elevation.m.poly[i],
        Elevation.m2 = top.mods.coeff$Elevation.m2.poly[i],
        Data.Type.nofi = top.mods.coeff$Data.TypeResurvey[i],
        Data.Type.Elevation.m.nofi = top.mods.coeff$`Data.TypeResurvey:Elevation.m.poly`[i],
        Data.Type.Elevation.m2.nofi = top.mods.coeff$`Data.TypeResurvey:Elevation.m2.poly`[i],
        Resurvey.Burned.fi = 
          top.mods.coeff$New.Data.TypeResurvey.Burned[i],
        Resurvey.Unburned.fi = 
          top.mods.coeff$New.Data.TypeResurvey.Unburned[i],
        Elevation.m.Res.Burn.fi = 
          top.mods.coeff$`Elevation.m.poly:New.Data.TypeResurvey.Burned`[i],
        Elevation.m.Res.Unburn.fi = 
          top.mods.coeff$`Elevation.m.poly:New.Data.TypeResurvey.Unburned`[i],
        Elevation.m2.Res.Burn.fi = 
          top.mods.coeff$`Elevation.m2.poly:New.Data.TypeResurvey.Burned`[i],
        Elevation.m2.Res.Unburn.fi = 
          top.mods.coeff$`Elevation.m2.poly:New.Data.TypeResurvey.Unburned`[i],
        row.names = NULL)
    }
    
    Mods <- ldply(Mods.list, data.frame)
    
    ## Storing model-averaged parameters
    
    Avg <- data.frame(
      Species = levels(species.list)[S],
      Dataset = D,
      L.Occ = sum(num.burns["1", , "Legacy"]), 
      R.Occ = sum(num.burns["1", , "Resurvey"]), 
      Fire.Included = ifelse(levels(species.list)[S] 
                             %in% species.with.fire$Species == TRUE, "Yes", "No"),
      Type = "Avg", 
      deltaAIC = NA, 
      Weight = NA,
      Rsquared = NA,
      Intercept = avg.mods.coeff$`(Intercept)`,
      Elevation.m = avg.mods.coeff$Elevation.m.poly,
      Elevation.m2 = avg.mods.coeff$Elevation.m2.poly,
      Data.Type.nofi = avg.mods.coeff$Data.TypeResurvey,
      Data.Type.Elevation.m.nofi = avg.mods.coeff$`Data.TypeResurvey:Elevation.m.poly`,
      Data.Type.Elevation.m2.nofi = avg.mods.coeff$`Data.TypeResurvey:Elevation.m2.poly`,
      Resurvey.Burned.fi = 
        avg.mods.coeff$New.Data.TypeResurvey.Burned,
      Resurvey.Unburned.fi = 
        avg.mods.coeff$New.Data.TypeResurvey.Unburned,
      Elevation.m.Res.Burn.fi = 
        avg.mods.coeff$`Elevation.m.poly:New.Data.TypeResurvey.Burned`,
      Elevation.m.Res.Unburn.fi = 
        avg.mods.coeff$`Elevation.m.poly:New.Data.TypeResurvey.Unburned`,
      Elevation.m2.Res.Burn.fi = 
        avg.mods.coeff$`Elevation.m2.poly:New.Data.TypeResurvey.Burned`,
      Elevation.m2.Res.Unburn.fi = 
        avg.mods.coeff$`Elevation.m2.poly:New.Data.TypeResurvey.Unburned`,
      row.names = NULL)
    
    coeff.ALLSPEC[[S]] <- rbind(Mods, Avg)
    
    
    ## Storing model-averaged CIs
    
    Avg.Confint <- data.frame(
      Species = levels(species.list)[S],
      Dataset = D,
      L.Occ = sum(num.burns["1", , "Legacy"]), 
      R.Occ = sum(num.burns["1", , "Resurvey"]), 
      Fire.Included = ifelse(levels(species.list)[S] 
                             %in% species.with.fire$Species == TRUE, "Yes", "No"),
      Type = "Avg.Confidence",
      # Lower CIs:
      Intercept.CI.Lower = avg.mods.confint$`(Intercept)`[1],
      Elevation.m.CI.Lower = avg.mods.confint$Elevation.m.poly[1],
      Elevation.m2.CI.Lower = avg.mods.confint$Elevation.m2.poly[1],
      Data.Type.nofi.CI.Lower = avg.mods.confint$Data.TypeResurvey[1],
      Data.Type.Elevation.m.nofi.CI.Lower = 
        avg.mods.confint$`Data.TypeResurvey:Elevation.m.poly`[1],
      Data.Type.Elevation.m2.nofi.CI.Lower = 
        avg.mods.confint$`Data.TypeResurvey:Elevation.m2.poly`[1],
      Resurvey.Burned.fi.CI.Lower = 
        avg.mods.confint$New.Data.TypeResurvey.Burned[1],
      Resurvey.Unburned.fi.CI.Lower = 
        avg.mods.confint$New.Data.TypeResurvey.Unburned[1],
      Elevation.m.Res.Burn.fi.CI.Lower = 
        avg.mods.confint$`Elevation.m.poly:New.Data.TypeResurvey.Burned`[1],
      Elevation.m.Res.Unburn.fi.CI.Lower = 
        avg.mods.confint$`Elevation.m.poly:New.Data.TypeResurvey.Unburned`[1],
      Elevation.m2.Res.Burn.fi.CI.Lower = 
        avg.mods.confint$`Elevation.m2.poly:New.Data.TypeResurvey.Burned`[1],
      Elevation.m2.Res.Unburn.fi.CI.Lower = 
        avg.mods.confint$`Elevation.m2.poly:New.Data.TypeResurvey.Unburned`[1],
      # Upper CIs: 
      Intercept.CI.Upper = avg.mods.confint$`(Intercept)`[2],
      Elevation.m.CI.Upper = avg.mods.confint$Elevation.m.poly[2],
      Elevation.m2.CI.Upper = avg.mods.confint$Elevation.m2.poly[2],
      Data.Type.nofi.CI.Upper = avg.mods.confint$Data.TypeResurvey[2],
      Data.Type.Elevation.m.nofi.CI.Upper = 
        avg.mods.confint$`Data.TypeResurvey:Elevation.m.poly`[2],
      Data.Type.Elevation.m2.nofi.CI.Upper = 
        avg.mods.confint$`Data.TypeResurvey:Elevation.m2.poly`[2],
      Resurvey.Burned.fi.CI.Upper = 
        avg.mods.confint$New.Data.TypeResurvey.Burned[2],
      Resurvey.Unburned.fi.CI.Upper = 
        avg.mods.confint$New.Data.TypeResurvey.Unburned[2],
      Elevation.m.Res.Burn.fi.CI.Upper = 
        avg.mods.confint$`Elevation.m.poly:New.Data.TypeResurvey.Burned`[2],
      Elevation.m.Res.Unburn.fi.CI.Upper = 
        avg.mods.confint$`Elevation.m.poly:New.Data.TypeResurvey.Unburned`[2],
      Elevation.m2.Res.Burn.fi.CI.Upper = 
        avg.mods.confint$`Elevation.m2.poly:New.Data.TypeResurvey.Burned`[2],
      Elevation.m2.Res.Unburn.fi.CI.Upper = 
        avg.mods.confint$`Elevation.m2.poly:New.Data.TypeResurvey.Unburned`[2],
      row.names = NULL)
    
    avg.confint.ALLSPEC[[S]] <- Avg.Confint
    
    ## Storing record of framework
    
    framework.ALLSPEC[[S]] <- framework.SPEC
    
  }  
  
  #Collapse ALLSPEC lists into single dataframe, store that df as part of ALLDAT list
  
  coeff.ALLDAT[[D]] <- ldply(coeff.ALLSPEC, data.frame)
  avg.confint.ALLDAT[[D]] <- ldply(avg.confint.ALLSPEC, data.frame)
  framework.ALLDAT[[D]] <- ldply(framework.ALLSPEC, data.frame)
  
  #### END OF SPECIES LOOP
  
  if(D == 100) system("say Your loop is done")
  
}

#### END OF DATASET LOOP

# Collate ALLDAT lists into one big DF

coeff.ALLDAT.allsets <- ldply(coeff.ALLDAT, data.frame)
avg.confint.ALLDAT.allsets <- ldply(avg.confint.ALLDAT, data.frame)
framework.ALLDAT.allsets <- ldply(framework.ALLDAT, data.frame)

## Remove error-laden datasets

# Create new column for joining
framework.ALLDAT.allsets$Join <- 
  paste(framework.ALLDAT.allsets$Dataset, framework.ALLDAT.allsets$Species, sep = ".")
discard.me <- framework.ALLDAT.allsets[ , c(6, 8)]
coeff.ALLDAT.allsets$Join <- 
  paste(coeff.ALLDAT.allsets$Dataset, coeff.ALLDAT.allsets$Species, sep = ".")
avg.confint.ALLDAT.allsets$Join <- 
  paste(avg.confint.ALLDAT.allsets$Dataset, avg.confint.ALLDAT.allsets$Species, sep = ".")

# Merge framework DF with coeff and confint DF
coeff.ALLDAT.allsets.join <- merge(coeff.ALLDAT.allsets, discard.me, 
                                   by = "Join")
avg.confint.ALLDAT.allsets.join <- merge(avg.confint.ALLDAT.allsets, discard.me, 
                                         by = "Join")

# Remove datasets associated with errors (Discard.Later == Yes)
coeff.ALLDAT.finaldf.big <- 
  coeff.ALLDAT.allsets.join[!coeff.ALLDAT.allsets.join$Discard.Later == "Yes", ]
avg.confint.ALLDAT.finaldf.big <- 
  avg.confint.ALLDAT.allsets.join[!avg.confint.ALLDAT.allsets.join$Discard.Later == "Yes", ]
coeff.ALLDAT.finaldf <- coeff.ALLDAT.finaldf.big[, c(2:22)]
avg.confint.ALLDAT.finaldf <- avg.confint.ALLDAT.finaldf.big[, c(2:31)]

# Store output as CSV

write.csv(coeff.ALLDAT.finaldf, 
          file = "data/east_only_check/eastonly_3c_new_coefficients.csv", 
          row.names = FALSE)
write.csv(avg.confint.ALLDAT.finaldf, 
          file = "data/east_only_check/eastonly_3c_new_confint.csv", 
          row.names = FALSE)
write.csv(framework.ALLDAT.allsets, 
          file = "data/east_only_check/eastonly_3c_new_framework_logs.csv", 
          row.names = FALSE)



################### SCRIPT 6 ANALOGUE ####################

#### This script visualizes model-averaged predictions for each rarefied dataset

library(tidyverse)
library(cowplot)

#### Read in and prepare tables of coefficients

# coefficients from all top models for each rarefied dataset that ran without warnings
# these models use the poly() formulaton for orthogonal elevation^2 terms
coeff.ALLDAT <- read.csv("data/east_only_check/eastonly_3c_new_coefficients.csv", header = TRUE)

# filter to just model average for each rarefied dataset
coeff.avgs <- coeff.ALLDAT %>% filter(Type=="Avg") # now we have up to 100 model averages per species (<100 for the species for which some rarefied datasets threw warnings)

# split into species modeled with fire vs without
coeffs.fire <- coeff.avgs %>% filter(Fire.Included=="Yes")
coeffs.nofire <- coeff.avgs %>% filter(Fire.Included=="No")

# species lists for ugly loops below
species.list.fire <- coeffs.fire %>% 
  group_by(Species) %>% 
  summarise(Species=first(Species))
species.list.nofire <- coeffs.nofire %>% 
  group_by(Species) %>% 
  summarise(Species=first(Species))


#### Elevation vectors for multiplying by coefficients

# values in poly-transformed units
dat <- read_csv("data/3d_transformed_polynomials.csv")

# linear term vector
# range based on min/max values 
elev.vec.lin = as.numeric(seq(min(dat$Elevation.m.poly), max(dat$Elevation.m.poly), by=0.0001)) 

# quadratic term vector
# quadratic function is given by this model
poly.mod <- lm(Elevation.m2.poly ~ Elevation.m.poly + I(Elevation.m.poly^2), data=dat)

elev.vec.quad = poly.mod$coefficients[1] + 
  elev.vec.lin*poly.mod$coefficients[2] + 
  elev.vec.lin*elev.vec.lin*poly.mod$coefficients[3]

# elevation list for back-transformed axis labels
# needs to be in raw units (m)
back.mod <- lm(Elevation.m.poly ~ Elevation.m, data=dat)
raw.ticks = c(100, 600, 1100, 1600, 2100)
poly.ticks = back.mod$coefficients[1] + raw.ticks*back.mod$coefficients[2]


#### other prep work

# function for converting predictions to 0-1 response scale
response = function(y) {
  exp(as.numeric(y))/(1+exp(as.numeric(y)))
}

# color palettes 
# for fire species
col.pal.fire <- c("turquoise4", "red3", "goldenrod1")
# for no-fire species
col.pal.nofire <- c("turquoise4", "goldenrod1")


#### FIRE SPECIES: big loop to calculate predicted values across each rarefaction for each species

# empty matrices for writing best-fit lines into
pred.leg.reps = matrix(nrow=length(elev.vec.lin),ncol=100)
pred.res.unburn.reps = matrix(nrow=length(elev.vec.lin),ncol=100)
pred.res.burn.reps = matrix(nrow=length(elev.vec.lin),ncol=100)

for (i in 1:dim(species.list.fire)[1]) {
  sp = as.list(species.list.fire[i,1]) 
  mods <- coeffs.fire %>% 
    filter(Species==sp) %>% 
    select(Int=Intercept, 
           Elev=Elevation.m, 
           Elev2=Elevation.m2, 
           ResurvBurn = Resurvey.Burned.fi,
           ResurvUnburn = Resurvey.Unburned.fi,
           ResurvBurnxElev = Elevation.m.Res.Burn.fi,
           ResurvBurnxElev2 = Elevation.m2.Res.Burn.fi,
           ResurvUnburnxElev = Elevation.m2.Res.Unburn.fi,
           ResurvUnburnXElev2 = Elevation.m2.Res.Unburn.fi) %>% 
    mutate_if(is.numeric, ~replace(., is.na(.), 0))
  
  for (j in 1:dim(mods)[1]) {
    pred.leg.reps[,j] = mods$Int[j] + mods$Elev[j]*elev.vec.lin + mods$Elev2[j]*elev.vec.quad
    
    pred.res.unburn.reps[,j] = mods$Int[j] + mods$Elev[j]*elev.vec.lin + mods$Elev2[j]*elev.vec.quad + mods$ResurvUnburn[j] + mods$ResurvUnburnxElev[j]*elev.vec.lin + mods$ResurvUnburnXElev2[j]*elev.vec.quad
    
    pred.res.burn.reps[,j] = mods$Int[j] + mods$Elev[j]*elev.vec.lin + mods$Elev2[j]*elev.vec.quad + mods$ResurvBurn[j] + mods$ResurvBurnxElev[j]*elev.vec.lin 
    + mods$ResurvBurnxElev2[j]*elev.vec.quad
  }
  
  t1.unburn.reps <- as.data.frame(cbind(elev.vec.lin, 'legacy', pred.leg.reps)) %>% 
    mutate(across(c(3:102), response))
  t1.unburn.reps.tall <- gather(t1.unburn.reps, "rep", "preds", 3:102) %>% 
    mutate(elev.vec.lin = as.numeric(elev.vec.lin))
  t1.unburn.summary <- t1.unburn.reps.tall %>% 
    group_by(V2, elev.vec.lin) %>% 
    summarise(mean.resp = mean(preds),
              lower.resp = unname(quantile(preds, c(0.05))),
              upper.resp = unname(quantile(preds, c(0.95))))
  t2.unburn.reps <- as.data.frame(cbind(elev.vec.lin, 'res.unburn', pred.res.unburn.reps)) %>% 
    mutate(across(c(3:102), response))
  t2.unburn.reps.tall <- gather(t2.unburn.reps, "rep", "preds", 3:102) %>% 
    mutate(elev.vec.lin = as.numeric(elev.vec.lin))
  t2.unburn.summary <- t2.unburn.reps.tall %>% 
    group_by(V2, elev.vec.lin) %>% 
    summarise(mean.resp = mean(preds),
              lower.resp = unname(quantile(preds, c(0.05))),
              upper.resp = unname(quantile(preds, c(0.95))))
  t2.burn.reps <- as.data.frame(cbind(elev.vec.lin, 'res.burn', pred.res.burn.reps)) %>% 
    mutate(across(c(3:102), response))
  t2.burn.reps.tall <- gather(t2.burn.reps, "rep", "preds", 3:102) %>% 
    mutate(elev.vec.lin = as.numeric(elev.vec.lin))
  t2.burn.summary <- t2.burn.reps.tall %>% 
    group_by(V2, elev.vec.lin) %>% 
    summarise(mean.resp = mean(preds),
              lower.resp = unname(quantile(preds, c(0.05))),
              upper.resp = unname(quantile(preds, c(0.95))))
  graph.dat.tall <- bind_rows(t1.unburn.reps.tall, t2.unburn.reps.tall, t2.burn.reps.tall)  
  graph.dat.means <- bind_rows(t1.unburn.summary, t2.unburn.summary, t2.burn.summary) %>% 
    mutate(preds = mean.resp)
  
  gg <- ggplot(graph.dat.means, aes(x = elev.vec.lin, y = preds, color = V2)) + 
    geom_line(data=graph.dat.tall, aes(group=interaction(V2, rep), color=V2), alpha=0.08, show.legend = FALSE) +
    geom_line(size=4, linetype="dotted", show.legend=FALSE) +
    theme_classic() +
    scale_color_manual("Time x fire", values=col.pal.fire, labels=c("legacy", "resurvey, burned", "resurvey, unburned")) +
    scale_x_continuous(breaks=poly.ticks, labels=raw.ticks) +
    xlab("") + #Elevation (m)
    ylab("") #Probability of presence
  
  ggsave(paste("figures/FigureS6_",sp,".pdf",sep=""), gg, width=5, height=5)
  
  assign(paste0("preds_graph_",sp), gg)
} 

### Loop writes graphs to PDF. Multi-panel figure compiled in Illustrator. View graphs below:

preds_graph_ACMI
preds_graph_ARUV
preds_graph_CARU
preds_graph_CEVE
preds_graph_EPAN
preds_graph_PAMY
preds_graph_VAME





















################# SCRIPT 5 ANALOGUE ###############


#### This script summarizes, analyzes, and visualizes arithmetic patterns in rarefied data; see script 6 for visualizing model outputs

### Load libraries
library(tidyverse)
library(Hmisc) #for mean_sdl function to put mean+/-SD on violin plots
library(cowplot) #for multipanel figures
library(viridis) #for color-vision-friendly palettes

### Load data
# rarefied data (as list of dataframes named rare.ALL):
load("data/rare.ALL.Rda")

### Species list

# start with all species
species.list <- c("ACMI", "ARUV", "CARU", "CEVE", "EPAN", "PAMY", "VAME")

# Load list of plot regions
regions <- read.csv("data/plot_regions.csv")
east.regions <- regions[regions$NAME == "East Cascades", 3:4]
east.list <- c(east.regions$Plot.Name.1980, east.regions$Plot.Name.2015)

# separate out fire species
species.list.fire <- read.csv("data/east_only_check/eastonly_3c_new_coefficients.csv", header = TRUE) %>% 
  filter(Type=="Avg")  %>% 
  filter(Fire.Included=="Yes") %>% 
  group_by(Species) %>% 
  summarise(Species=first(Species))


### Set up empty matrices to store values
# As above, for fire species
el.mins.025.leg.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.maxs.975.leg.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.meds.leg.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.mins.025.res.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.maxs.975.res.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.meds.res.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.highshift.975.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.lowshift.025.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])

### Clunky for-loops across rarefied datasets by species

## fire species
for(D in 1:100) { #Run time: ~10 seconds
  und <- rare.ALL[[D]]
  und <- und[und$Plot %in% east.list, ]
  und$Data.Type <- as.factor(und$Data.Type)
  und$Elevation.m <- as.numeric(und$Elevation.m)
  for(S in 1:dim(species.list.fire)[1]) {
    species <- as.list(species.list.fire[S,1])
    und.SPEC <- und %>% 
      filter(Species.Code == species) %>% 
      droplevels()
    und.presence.SPEC.leg <- und.SPEC %>% 
      filter(Data.Type=="Legacy" & Pres.Abs==1) %>% 
      mutate(el.med.leg.fire = median(Elevation.m),
             el.min.025.leg.fire = quantile(Elevation.m, probs=0.025),
             el.max.975.leg.fire = quantile(Elevation.m, probs=0.975))
    und.presence.SPEC.res <- und.SPEC %>% 
      filter(Data.Type=="Resurvey" & Pres.Abs==1) %>% 
      mutate(el.med.res.fire = median(Elevation.m),
             el.min.025.res.fire = quantile(Elevation.m, probs=0.025),
             el.max.975.res.fire = quantile(Elevation.m, probs=0.975))
    el.mins.025.leg.fire[D,S] <- und.presence.SPEC.leg$el.min.025.leg.fire[1]
    el.maxs.975.leg.fire[D,S] <- und.presence.SPEC.leg$el.max.975.leg.fire[1]
    el.meds.leg.fire[D,S] <- und.presence.SPEC.leg$el.med.leg.fire[1]
    el.mins.025.res.fire[D,S] <- und.presence.SPEC.res$el.min.025.res.fire[1]
    el.maxs.975.res.fire[D,S] <- und.presence.SPEC.res$el.max.975.res.fire[1]
    el.meds.res.fire[D,S] <- und.presence.SPEC.res$el.med.res.fire[1]
    el.highshift.975.fire[D,S] <- und.presence.SPEC.res$el.max.975.res.fire[1] - und.presence.SPEC.leg$el.max.975.leg.fire[1]
    el.lowshift.025.fire[D,S] <- und.presence.SPEC.res$el.min.025.res.fire[1] - und.presence.SPEC.leg$el.min.025.leg.fire[1]
  }
}

## collapse to means across rarefied replicates

# fire species
el.mins.025.leg.means.fire <- as.data.frame(el.mins.025.leg.fire) %>% summarise(across(starts_with("V"),mean))
el.mins.025.res.means.fire <- as.data.frame(el.mins.025.res.fire) %>% summarise(across(starts_with("V"),mean))
el.maxs.975.leg.means.fire <- as.data.frame(el.maxs.975.leg.fire) %>% summarise(across(starts_with("V"),mean))
el.maxs.975.res.means.fire <- as.data.frame(el.maxs.975.res.fire) %>% summarise(across(starts_with("V"),mean))
el.meds.leg.means.fire <- as.data.frame(el.meds.leg.fire) %>% summarise(across(starts_with("V"),mean))
el.meds.res.means.fire <- as.data.frame(el.meds.res.fire) %>% summarise(across(starts_with("V"),mean))
el.highshift.means.fire <- as.data.frame(el.highshift.975.fire) %>% summarise(across(starts_with("V"),mean))
el.highshift.sd.fire <- as.data.frame(el.highshift.975.fire) %>% summarise(across(starts_with("V"),sd))
el.lowshift.means.fire <- as.data.frame(el.lowshift.025.fire) %>% summarise(across(starts_with("V"),mean))
el.lowshift.sd.fire <- as.data.frame(el.lowshift.025.fire) %>% summarise(across(starts_with("V"),sd))

## reshape and join into one frame

# fire species
el.mins.025.leg.tall.fire <- gather(el.mins.025.leg.means.fire, "species", "min.025.leg", 1:dim(species.list.fire)[1])
el.mins.025.res.tall.fire <- gather(el.mins.025.res.means.fire, "species", "min.025.res", 1:dim(species.list.fire)[1])
el.maxs.975.leg.tall.fire <- gather(el.maxs.975.leg.means.fire, "species", "max.975.leg", 1:dim(species.list.fire)[1])
el.maxs.975.res.tall.fire <- gather(el.maxs.975.res.means.fire, "species", "max.975.res", 1:dim(species.list.fire)[1])
el.meds.leg.tall.fire <- gather(el.meds.leg.means.fire, "species", "med.leg", 1:dim(species.list.fire)[1])
el.meds.res.tall.fire <- gather(el.meds.res.means.fire, "species", "med.res", 1:dim(species.list.fire)[1])
el.mean.highshift.95.tall.fire <- gather(el.highshift.means.fire, "species", "mean.high", 1:dim(species.list.fire)[1])
el.sd.highshift.95.tall.fire <- gather(el.highshift.sd.fire, "species", "sd.high", 1:dim(species.list.fire)[1])
el.mean.lowshift.95.tall.fire <- gather(el.lowshift.means.fire, "species", "mean.low", 1:dim(species.list.fire)[1])
el.sd.lowshift.95.tall.fire <- gather(el.lowshift.sd.fire, "species", "sd.low", 1:dim(species.list.fire)[1])

# for plotting and t-tests
rarefied.change.fire <- left_join(left_join(left_join(left_join(left_join(el.meds.leg.tall.fire, el.meds.res.tall.fire), el.mins.025.leg.tall.fire), el.mins.025.res.tall.fire), el.maxs.975.leg.tall.fire), el.maxs.975.res.tall.fire)
rarefied.change.fire <- cbind(rarefied.change.fire, species.list.fire)
rarefied.change.fire$fire <- "yes"

# for supp table
rarefied.change.fire.supp <- left_join(left_join(left_join(el.mean.highshift.95.tall.fire, el.sd.highshift.95.tall.fire),el.mean.lowshift.95.tall.fire), el.sd.lowshift.95.tall.fire)
rarefied.change.fire.supp <- cbind(rarefied.change.fire.supp, species.list.fire)
rarefied.change.fire.supp$fire <- "yes"

## calculate range changes

# fire species
rarefied.change.fire <- rarefied.change.fire %>% 
  mutate(rear.change.perc = min.025.res - min.025.leg,
         med.change = med.res - med.leg,
         lead.change.perc = max.975.res - max.975.leg)

rarefied.change.calcs <- rarefied.change.fire

rarefied.change.calcs.supp <- rarefied.change.fire.supp


### FIGURE 3: VIOLIN PLOTS

# fire species
rarefied.change.tall.fire.perc <- rarefied.change.calcs %>% 
  filter(fire=="yes") %>% 
  select(rear.change.perc, med.change, lead.change.perc) %>% 
  gather("edge", "change", 1:3) %>% 
  add_column(fire="yes")

rarefied.change.tall.perc <- rarefied.change.tall.fire.perc

level_order.perc = c("rear.change.perc", "med.change", "lead.change.perc")

# Violins by range position - spits errors but display is ok
violin.plot.perc <- ggplot(rarefied.change.tall.perc, aes(x=factor(edge, level=level_order.perc), y=change, fill=fire)) + 
  geom_violin(aes(fill=fire), trim=FALSE, position="dodge", fill = "lightyellow", color = "black") +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange", position=position_dodge(0.9), cex=0.7)  +
  theme_classic() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, position=position_jitterdodge(jitter.width=0.1, jitter.height=10, dodge.width=0.9)) +
  geom_hline(yintercept=0, lty="dashed") +
  xlab("") +
  scale_x_discrete(labels=c("Lower\nedge", "Range\ncenter", "Upper\nedge")) +
  theme(axis.text.x=element_text(size=12)) +
  ylim(-600,700) +
  ylab(c("Elevational change (m)\n1983-2015")) + 
  guides(fill = FALSE)
violin.plot.perc

ggsave("figures/FigureS5_eastonly_violin_1panel_perc.pdf", violin.plot.perc, device="pdf", width=8, height=5) 

### FIGURE 2: Freeman-style elevation ranges

# fire and no-fire species separated
rarefied.change.calcs.fire <- rarefied.change.calcs %>% 
  filter(fire=="yes") %>% 
  mutate(species.rank.med = dense_rank(med.leg)+35, #silly workaround to get unique 4-letter species codes on fire species
         both.min.perc = pmax(min.025.leg, min.025.res),
         both.max.perc = pmin(max.975.leg, max.975.res))

rarefied.change.calcs.facet <- rarefied.change.calcs.fire

species.labels <- rarefied.change.calcs.facet$Species[order(rarefied.change.calcs.facet$fire, rarefied.change.calcs.facet$species.rank.med)]

fire.labs <- c("non-fire-experiencing species", "fire-experiencing (east only)")
names(fire.labs) <- c("no", "yes")

p.facet <- ggplot(rarefied.change.calcs.facet) + 
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=min.025.leg, ymax=max.975.leg), fill = "#F8766D") + # historic range in red; will show areas of range contractions
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=min.025.res, ymax=max.975.res), fill = "#00BFC4") + # modern range in blue; will show areas of range expansions
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=both.min.perc, ymax=both.max.perc), fill = "#bdbdbd") + # areas common to both in grey
  geom_text(aes(x = species.rank.med, y = 0, label = Species), vjust = -0.8, size = 4, angle = 45) +  # Add this line to add labels to the bars
  facet_grid(. ~ fire, scale="free", space="free", labeller=labeller(fire=fire.labs)) +
  ylim(0,2200) +
  xlab("Species") +
  ylab("Elevation (m)") +
  theme_bw() +
  theme(text=element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=90, hjust=0, vjust=0.5)) +
  scale_x_discrete()
p.facet

ggsave("figures/FigureS4_elevation_ranges.pdf", p.facet, device="pdf", width=3.5, height=6) 






