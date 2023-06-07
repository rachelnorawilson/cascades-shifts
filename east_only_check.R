# Created: June 3, 2023

# This script will NOT produce the final version of the analysis.
# It is solely for determining whether subsetting to east-side plots affects
# fire-experiencing species patterns. It is close to a copy of script 3b.
# However, this script differs from 3b in that it forces all species to conform to fire models.

# Part of this script creates a large file and places it in a folder
# "warning_logs"; these are not stored on GitHub but can be shared upon request and/or
# reproduced locally. User will need to adjust the destination directory (lines 102 and 104).

# No CIs are generated in this script.

# IMPORTANT NOTE: unless otherwise indicated, always use Understory_All.csv for these analyses as it is the ONLY file with up-to-date corrections.

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


#### STEP 2: Loop to analyze presence data - FIRE SPECIES ONLY ####

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

# write.csv(coeff.ALLDAT.finaldf, 
#          file = "data/3_presence_ALLDAT_ALLSPEC_coefficients.csv", 
#          row.names = FALSE)
# write.csv(warn.ALLDAT.finaldf, 
#          file = "data/3_presence_ALLDAT_ALLSPEC_warnings.csv", 
#          row.names = FALSE)



###################### STOP ########################################3



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


#### STEP 1: Import data ####

warn.ALLDAT <- warn.ALLDAT.finaldf

# Can use either top_mod input (script stored outside of GitHub) or new_coefficients input
coeff.ALLDAT <- coeff.ALLDAT.finaldf
coeff.ALLDAT[is.na(coeff.ALLDAT)] <- paste(0)


#### STEP 2 (OPTIONAL): Exploratory visualizations of warnings ####

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


#### STEP 3: Exploring coefficients ####

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

