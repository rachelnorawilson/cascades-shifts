# Created: June 3, 2023

# This script will NOT produce the final version of the analysis.
# It is solely for determining whether subsetting to east-side plots affects
# fire-experiencing species patterns. It includes versions of 3b* (to generate error logs)
# and script 3c (to run models and create coefficients CSV used in vis)
## *this script differs from 3b in that it forces all species to conform to fire models.

# Part of this script (3b) creates a large file and places it in a folder
# "warning_logs"; these are not stored on GitHub but can be shared upon request and/or
# reproduced locally. User will need to adjust the destination directory (lines 102 and 104).

# No CIs are generated in this script.

# IMPORTANT NOTE: unless otherwise indicated, always use Understory_All.csv for these analyses as it is the ONLY file with up-to-date corrections.

################# SCRIPT 3B ANALOGUE ###############

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


# Store output as CSV

write.csv(coeff.ALLDAT.finaldf, 
          file = "data/TEMP_eastonly_3_presence_ALLDAT_ALLSPEC_coefficients.csv", 
          row.names = FALSE)
write.csv(warn.ALLDAT.finaldf, 
          file = "data/TEMP_eastonly_3_presence_ALLDAT_ALLSPEC_warnings.csv", 
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

# Porting in old CSV of warnings
warn.ALLDAT <- read.csv("data/TEMP_eastonly_3_presence_ALLDAT_ALLSPEC_warnings.csv", header = TRUE)
# Bug correction - force VAME and HODI to run under correct framework
warn.ALLDAT$Has_warning[warn.ALLDAT$Species == "VAME" 
                        & warn.ALLDAT$Dataset == 36] <- paste("TRUE")

# Formula to correct erroneous burn coding. Run time 3 sec
for(D in 1:100) {
  rare.ALL[[D]]$Fires[rare.ALL[[D]]$Data.Type == "Legacy"] <- rep("Unburned")
  rare.ALL[[D]]$New.Data.Type <- paste(rare.ALL[[D]]$Data.Type, rare.ALL[[D]]$Fires, sep = ".")
}


#### STEP 2: Loop to analyze presence data ####

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
  und.presence <- rare.ALL[[D]]
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
          file = "data/3c_new_coefficients.csv", 
          row.names = FALSE)
write.csv(avg.confint.ALLDAT.finaldf, 
          file = "data/3c_new_confint.csv", 
          row.names = FALSE)
write.csv(framework.ALLDAT.allsets, 
          file = "data/3c_new_framework_logs.csv", 
          row.names = FALSE)



###################### SCRIPT 4A ANALOGUE ########################################



# Everything below is essentially a carbon copy of script 4a, used for reviewing model fitting errors

# Functions needed:

# Turning averaged coefficients into +/-/0
simplify.fun <- function(varib) {
  coeff.summary.SPEC <- coeff.summary.empty[varib]
  dataset.vect <- as.vector(coeff.SPEC.avg$Dataset)
  for(d in dataset.vect) {
    if(coeff.SPEC.avg[coeff.SPEC.avg$Dataset == d, varib] > 0) {
      coeff.summary.SPEC[d, varib] <- paste("+")
    }
    if(coeff.SPEC.avg[coeff.SPEC.avg$Dataset == d, varib] == 0) {
      coeff.summary.SPEC[d, varib] <- paste("0")
    }
    if(coeff.SPEC.avg[coeff.SPEC.avg$Dataset == d, varib] < 0) {
      coeff.summary.SPEC[d, varib] <- paste("-")
    }
  }
  return(coeff.summary.SPEC)
}

# Summing numbers of +/-/0 for each coeff for one species
# Empty data frame
coeff.count.empty <- data.frame(Elevation.m = rep(NA, times = 4), 
                                Elevation.m2 = rep(NA, times = 4), 
                                Resurvey.Burned.fi = rep(NA, times = 4), 
                                Resurvey.Unburned.fi = rep(NA, times = 4), 
                                Elevation.m.Res.Burn.fi = rep(NA, times = 4), 
                                Elevation.m.Res.Unburn.fi = rep(NA, times = 4), 
                                Elevation.m2.Res.Burn.fi = rep(NA, times = 4), 
                                Elevation.m2.Res.Unburn.fi = rep(NA, times = 4), 
                                Data.Type.nofi = rep(NA, times = 4), 
                                Data.Type.Elevation.m.nofi = rep(NA, times = 4), 
                                Data.Type.Elevation.m2.nofi = rep(NA, times = 4),
                                row.names = c("+", "-", 0, "Ignore"))

# Collapsing 100 rows of +/-/0 into summary
count.fun <- function(varib) {
  df.count <- coeff.count.empty[varib]
  if(is.na(table(simple.coeffs.SPEC[varib])["+"]) == FALSE) {
    df.count["+", varib] <- table(simple.coeffs.SPEC[varib])["+"]
  } 
  if(is.na(table(simple.coeffs.SPEC[varib])["-"]) == FALSE) {
    df.count["-", varib] <- table(simple.coeffs.SPEC[varib])["-"]
  } 
  if(is.na(table(simple.coeffs.SPEC[varib])["0"]) == FALSE) {
    df.count["0", varib] <- table(simple.coeffs.SPEC[varib])["0"]
  } 
  if(is.na(table(simple.coeffs.SPEC[varib])["Ignore"]) == FALSE) {
    df.count["Ignore", varib] <- table(simple.coeffs.SPEC[varib])["Ignore"]
  } 
  return(df.count)
}


#### STEP 1: Import data

warn.ALLDAT <- warn.ALLDAT.finaldf

# Can use either top_mod input (script stored outside of GitHub) or new_coefficients input
coeff.ALLDAT <- coeff.ALLDAT.finaldf
coeff.ALLDAT[is.na(coeff.ALLDAT)] <- paste(0)


#### STEP 2 (OPTIONAL): Exploratory visualizations of warnings

(numbered.species <- data.frame(Species=species.list, No.=rep(1:length(species.list))))

# Choose species of interest
warn.SPEC <- warn.ALLDAT[warn.ALLDAT$Species == "ACMI", ]

# Assess  warnings
table(warn.SPEC$Has_warning, warn.SPEC$Dataset)
table(table(warn.SPEC$Has_warning, warn.SPEC$Dataset))
head(warn.SPEC[warn.SPEC$Has_warning == TRUE, ])

# (Optional): How did rarefaction vary between species? #
#par(mar=c(5, 1.5, 1, 0), oma=c(0,4,0,0)) # 3-paneled graph
#hist(coeff.SPEC$R.Occ, main = "RHAL", xlab = "No. resurvey occ", ylab = "Frequency")


#### STEP 3: Exploring coefficients

coeff.count.LIST <- list()
coeff.perc.LIST <- list()

for(S in 1:nrow(numbered.species)) { # Run time: 5 sec
  #coeff.SPEC.avg <- coeff.ALLDAT[coeff.ALLDAT$Species == levels(numbered.species$Species)[S], ]
  coeff.SPEC.avg <- coeff.ALLDAT[coeff.ALLDAT$Species == levels(numbered.species$Species)[S] & 
                                   coeff.ALLDAT$Type == "Avg", ]
  
  coeff.summary.empty <- data.frame(Elevation.m = rep(NA, times = 100), 
                                    Elevation.m2 = rep(NA, times = 100), 
                                    Resurvey.Burned.fi = rep(NA, times = 100), 
                                    Resurvey.Unburned.fi = rep(NA, times = 100), 
                                    Elevation.m.Res.Burn.fi = rep(NA, times = 100), 
                                    Elevation.m.Res.Unburn.fi = rep(NA, times = 100), 
                                    Elevation.m2.Res.Burn.fi = rep(NA, times = 100), 
                                    Elevation.m2.Res.Unburn.fi = rep(NA, times = 100), 
                                    Data.Type.nofi = rep(NA, times = 100), 
                                    Data.Type.Elevation.m.nofi = rep(NA, times = 100), 
                                    Data.Type.Elevation.m2.nofi = rep(NA, times = 100))
  
  simple.coeffs.SPEC <- bind_cols(lapply(names(coeff.summary.empty), simplify.fun))
  simple.coeffs.SPEC[is.na(simple.coeffs.SPEC)] <- paste("Ignore")
  
  # Summarizing the summary
  
  coeff.count.SPEC <- bind_cols((lapply(names(coeff.count.empty), count.fun)))
  coeff.count.SPEC$Species <- levels(factor(coeff.SPEC.avg$Species))
  if(length(levels(factor(coeff.SPEC.avg$Fire.Included))) == 1) {
    coeff.count.SPEC$Fire.Included <- levels(factor(coeff.SPEC.avg$Fire.Included))
  } else(coeff.count.SPEC$Fire.Included <- "Sometimes")
  coeff.count.SPEC$Effect <- row.names(coeff.count.SPEC)
  coeff.count.LIST[[S]] <- coeff.count.SPEC
  
  # Converting to % of good datasets
  good.sets.SPEC <- 
    sum(coeff.count.SPEC$Elevation.m[!coeff.count.SPEC$Effect == "Ignore"], na.rm = TRUE)
  coeff.perc.SPEC.A <- data.frame(Species = coeff.count.SPEC$Species,
                                  Fire.Included = coeff.count.SPEC$Fire.Included,
                                  Total.Sets = rep(good.sets.SPEC, times = nrow(coeff.count.SPEC)),
                                  Effect = coeff.count.SPEC$Effect, 
                                  Elevation.m = coeff.count.SPEC$Elevation.m / good.sets.SPEC,
                                  Elevation.m2 = coeff.count.SPEC$Elevation.m2 / good.sets.SPEC, 
                                  Resurvey.Burned.fi = 
                                    coeff.count.SPEC$Resurvey.Burned.fi / good.sets.SPEC, 
                                  Resurvey.Unburned.fi = 
                                    coeff.count.SPEC$Resurvey.Unburned.fi / good.sets.SPEC, 
                                  Elevation.m.Res.Burn.fi = 
                                    coeff.count.SPEC$Elevation.m.Res.Burn.fi / good.sets.SPEC, 
                                  Elevation.m.Res.Unburn.fi = 
                                    coeff.count.SPEC$Elevation.m.Res.Unburn.fi / good.sets.SPEC, 
                                  Elevation.m2.Res.Burn.fi = 
                                    coeff.count.SPEC$Elevation.m2.Res.Burn.fi / good.sets.SPEC, 
                                  Elevation.m2.Res.Unburn.fi = 
                                    coeff.count.SPEC$Elevation.m2.Res.Unburn.fi / good.sets.SPEC, 
                                  Data.Type.nofi =
                                    coeff.count.SPEC$Data.Type.nofi / good.sets.SPEC, 
                                  Data.Type.Elevation.m.nofi =
                                    coeff.count.SPEC$Data.Type.Elevation.m.nofi / good.sets.SPEC,  
                                  Data.Type.Elevation.m2.nofi = 
                                    coeff.count.SPEC$Data.Type.Elevation.m2.nofi / good.sets.SPEC)
  coeff.perc.SPEC.B <- coeff.perc.SPEC.A[!coeff.perc.SPEC.A$Effect == "Ignore", ]
  coeff.perc.LIST[[S]] <- coeff.perc.SPEC.B
  
} # End of loop

coeff.count <- ldply(coeff.count.LIST, data.frame)
coeff.perc <- ldply(coeff.perc.LIST, data.frame)

#write.csv(coeff.count, 
 #         file = "data/4a_presence_coefficients_count.csv", 
  #        row.names = FALSE)


#write.csv(coeff.perc.fire, 
 #         file = "data/4a_pres_coefficients_percent_fire.csv", row.names = FALSE)
#write.csv(coeff.perc.nofire, 
      #    file = "data/4a_pres_coefficients_percent_NOfire.csv", row.names = FALSE)





##################### STOP #######################

# TBD: visualizations

