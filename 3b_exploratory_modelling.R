# Created: Dec. 11, 2020
# Updated: Feb. 22, 2021

# This script will be used to undertake part b of the PRESENCE analyses (modeling)
# It will NOT produce the final version of the analysis, but is instead used for logging 
# model warnings and framework-switching (fire vs. no fire). 

# Part of this script creates a ~16 MB worth of 4,200+ text logs and places them in a folder
# "warning_logs"; these are not stored on GitHub but can be shared upon request and/or
# reproduced locally. User will need to adjust the destination directory (lines 102 and 104).

# No CIs are generated in this script.

# IMPORTANT NOTE: unless otherwise indicated, always use Understory_All.csv for these analyses as it is the ONLY file with up-to-date corrections.

# Packages needed:

library(MuMIn)
library(plyr)
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
species.list <- shifts$Species.Code[!shifts$Species.Code=="MOSS"] #removing "MOSS"
species.list <- factor(species.list)

#TODO Formula to correct erroneous burn coding. Run time 3 sec
for(D in 1:100) {
  rare.ALL[[D]]$Fires[rare.ALL[[D]]$Data.Type == "Legacy"] <- rep("Unburned")
  rare.ALL[[D]]$New.Data.Type <- paste(rare.ALL[[D]]$Data.Type, rare.ALL[[D]]$Fires, sep = ".")
}


#### STEP 2: Loop to analyze presence data ####

# Can be run as a loop outputting all species, or S can be modified to isolate specific species. Check number here:
(numbered.species <- data.frame(Species=species.list, No.=rep(1:42)))
# Exclude D loop to run for only one dataset (e.g. unrarefied data)

coeff.ALLDAT <- list()
warn.ALLDAT <- list()

for(D in 1:100) { #RUN TIME: 32 minutes 51 sec
  
  coeff.ALLSPEC <- list()
  warn.ALLSPEC <- list()
  und.presence <- rare.ALL[[D]]
  und.presence$Elevation.m <- as.numeric(und.presence$Elevation.m)
  und.presence$New.Data.Type <- factor(und.presence$New)
  und.presence$Fires <- as.factor(und.presence$Fires)
  und.presence$Data.Type <- as.factor(und.presence$Data.Type)
  
  for(S in 1:length(species.list)) { # RUN TIME: ~15-20 sec for 1 dataset
    
    # Create subset for species of interest S
    und.presence.SPEC <- subset(und.presence, Species.Code == levels(species.list)[S])
    und.presence.SPEC$New.Data.Type <- relevel(und.presence.SPEC$New.Data.Type, 
                                               ref="Legacy.Unburned")
    und.presence.SPEC$Fires <- relevel(und.presence.SPEC$Fires, ref="Unburned")
    und.presence.SPEC$Data.Type <- relevel(und.presence.SPEC$Data.Type, ref="Legacy")
    und.presence.SPEC$Elevation.m2 <- und.presence.SPEC$Elevation.m^2
    und.presence.SPEC <- und.presence.SPEC[complete.cases(und.presence.SPEC), ] #Just in case
    
    # Emptying out previous model objects
    mod.globfi <- NULL
    mod.globnofi <- NULL
    
    # Did the species occur in burned plots 5+ times?
    num.burns <- 
      table(und.presence.SPEC$Pres.Abs, und.presence.SPEC$Fires, und.presence.SPEC$Data.Type)
    
    
    # Create new error-logging file specific to permutation & species
    log.file <- file(paste("data/warning_logs/warnings", 
                           D, levels(species.list)[S], ".txt", sep = "_"), open = "wt")
    log.file.path <- paste("data/warning_logs/warnings", 
                           D, levels(species.list)[S], ".txt", sep = "_")
    sink(log.file, append = TRUE, type = "output") # Sink to log file
    sink(log.file, append = TRUE, type = "message")
    
    
    # If burn > 5, include fire as a predictor in global model (New.Data.Type):
    if(num.burns["1", "Burned", "Resurvey"] >= 5) {
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
    }
    
    
    #If burn < 5 in both, exclude fire from global model:
    if(num.burns["1", "Burned", "Resurvey"] < 5) {
      mod.globnofi <- glm(Pres.Abs ~ Data.Type * (Elevation.m + Elevation.m2), 
                          data = und.presence.SPEC, family = "binomial", na.action = na.fail) 
      options(warn = -1) # Ignore warnings - not for logging.
      dredge.list <- lapply(dredge(mod.globnofi, rank = AIC, subset = 
                                     dc(Elevation.m, Elevation.m2) &&
                                     dc(Data.Type:Elevation.m, Data.Type:Elevation.m2), 
                                   trace = FALSE, evaluate = FALSE), eval)
      names(dredge.list) <- paste("Mod", 
                                  as.numeric(names(dredge.list)) - 1, 
                                  sep = ".") # Converting to model ID
      options(warn = 1) # Tell me if a model throws an error - for logging.
      dredge.globnofi <- dredge(mod.globnofi, rank = AIC, subset = 
                                  dc(Elevation.m, Elevation.m2) &&
                                  dc(Data.Type:Elevation.m, Data.Type:Elevation.m2), 
                                trace = 1)
    }
    
    
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

# write.csv(coeff.ALLDAT.finaldf, 
#          file = "data/3_presence_ALLDAT_ALLSPEC_coefficients.csv", 
#          row.names = FALSE)
# write.csv(warn.ALLDAT.finaldf, 
#          file = "data/3_presence_ALLDAT_ALLSPEC_warnings.csv", 
#          row.names = FALSE)





