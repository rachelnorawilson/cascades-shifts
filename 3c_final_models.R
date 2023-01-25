# Created: Feb. 22, 2021
# Updated: Mar. 15, 2021

# This script will be used to undertake part c of the PRESENCE analyses (modeling)
# Use this script to produce the FINAL version of the analysis, plus CIs. No error logging.

# This script accommodates the following species-specific issues:
# --> Discard species for whom 2 or more models threw warnings
# # ---> COST, LUPE, PHEM, RHAL, VAAL, VADE
# --> Exclude the most-complex model (X * elev^2) for species in which this model threw warnings
# # ---> HODI, VAME*
# --> For species that flip-flopped between yes/no fire, use "majority rules" to decide
#     which framework to use
# # ---> No fire: AMAL, SPBE
# # ---> Fire: ARUV, CARU, VAME*
# --> Exclude problematic sets when remainder threw no warnings 
# # ---> CAME, GAOV, HIAL, RULA, TRBO
# As well as the following general changes:
# --> Switch elevation^2 term to poly(elevation)

# IMPORTANT NOTE: unless otherwise indicated, always use Understory_All.csv for these analyses as it is the ONLY file with up-to-date corrections.

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
species.list <- shifts$Species.Code[!shifts$Species.Code=="MOSS" &
                                      !shifts$Species.Code=="COST" &
                                      !shifts$Species.Code=="LUPE" &
                                      !shifts$Species.Code=="PHEM" &
                                      !shifts$Species.Code=="RHAL" &
                                      !shifts$Species.Code=="VAAL" &
                                      !shifts$Species.Code=="VADE"] #removing problematic species
species.list <- factor(species.list)
warn.ALLDAT <- read.csv("data/3_presence_ALLDAT_ALLSPEC_warnings.csv", header = TRUE)
# Bug correction - force VAME and HODI to run under correct framework
warn.ALLDAT$Has_warning[warn.ALLDAT$Species == "VAME" 
                        & warn.ALLDAT$Dataset == 36] <- paste("TRUE")
warn.ALLDAT$Has_warning[warn.ALLDAT$Species == "HODI" 
                        & warn.ALLDAT$Dataset == c(31, 54, 61, 87)] <- paste("TRUE")

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







