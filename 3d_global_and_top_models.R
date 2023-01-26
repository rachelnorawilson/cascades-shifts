# Created: Mar. 8, 2021
# Updated: Apr. 12, 2021

# This script will be used to undertake part d of the PRESENCE analyses (modeling)
# Use this script to produce the top model and global model outputs
# This script also stores all model outputs as Rda files, amounting to 300+ MB.
# These outputs will not be stored on GitHub but can be provided upon request.
# Users can also generate the outputs and store them locally (change directory in line 480)

# NOTE: This script has been modified to use ORTHO polynomials

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
coeff.all.P <- paste( "P", coeff.all, sep = ".")

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

coeff.ALLDAT <- list() # Store coefficient outputs from top model
framework.ALLDAT <- list() # Store record of which decision framework was used
global.wP.ALLDAT <- list() # Store coefficient outputs from global model

for(D in 1:100) { #RUN TIME: 4 min
  
  coeff.ALLSPEC <- list()
  framework.ALLSPEC <- list()
  global.wP.ALLSPEC <- list()
  mods.ALLSPEC <- list()
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
    # und.presence.SPEC$Elevation.m2 <- und.presence.SPEC$Elevation.m^2
    # Updated Mar. 27: Replace Elevation.m w/ poly()
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
    top.mods.coeff <- NULL
    global.coeff.wP <- NULL
    global.confint <- NULL
    
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
          mods.ALLSPEC[[levels(species.list)[S]]] <- mod.globfi.reduced
          dredge.globfi.reduced <- dredge(mod.globfi.reduced, rank = AIC, subset = 
                                            dc(Elevation.m.poly, Elevation.m2.poly))
          
          # Coefficients:
          top.mods.coeff <- as.data.frame(coef(subset(dredge.globfi.reduced, delta == 0)))
          top.mods.coeff$logLik <- dredge.globfi.reduced$logLik[dredge.globfi.reduced$delta == 0]
          # Global coefficients:
          global.coeff <- as.data.frame(t(coef(mod.globfi.reduced)))
          global.coeff$logLik <- as.numeric(logLik(mod.globfi.reduced))
          global.coeff.P <- as.data.frame(t(summary(mod.globfi.reduced)$coefficients[,4]))
          #TODO - fix confint issue!
          global.confint <- as.data.frame(t(confint(mod.globfi.reduced)))
          # Attaching global.coeff.P to global.coeff
          names(global.coeff.P) <- paste( "P", names(global.coeff.P), sep = ".")
          global.coeff.wP <- cbind(global.coeff, global.coeff.P)
          
          framework.SPEC$Forced.Simpler.Mod <- paste("Yes") # Record
          
        } else {
          framework.SPEC$Discard.Later <- paste("Yes") # Record
          # Create empty DFs
          top.mods.coeff <- data.frame(logLik = rep(NA, 1))
          global.coeff.wP <- data.frame(logLik = rep(NA, 1))
          global.confint <- data.frame(Elevation.m.poly = rep(NA, 2), 
                                       row.names = c("2.5 %", "97.5 %"))
        }
        
        
      } else { # Run as normal
        mod.globfi <- glm(Pres.Abs ~ (Elevation.m.poly + Elevation.m2.poly) * New.Data.Type, 
                          data = und.presence.SPEC, family = "binomial", na.action = na.fail)
        mods.ALLSPEC[[levels(species.list)[S]]] <- mod.globfi
        dredge.globfi <- dredge(mod.globfi, rank = AIC, subset = 
                                  dc(Elevation.m.poly, Elevation.m2.poly) &&
                                  dc(New.Data.Type:Elevation.m.poly, 
                                     New.Data.Type:Elevation.m2.poly) &&
                                  dc(Elevation.m.poly:New.Data.Type, 
                                     Elevation.m2.poly:New.Data.Type))
        
        # Storing coefficients:
        top.mods.coeff <- as.data.frame(coef(subset(dredge.globfi, delta == 0)))
        top.mods.coeff$logLik <- dredge.globfi$logLik[dredge.globfi$delta == 0]
        # Global coefficients:
        global.coeff <- as.data.frame(t(coef(mod.globfi)))
        global.coeff$logLik <- as.numeric(logLik(mod.globfi))
        global.coeff.P <- as.data.frame(t(summary(mod.globfi)$coefficients[,4]))
        #TODO - fix confint issue!
        global.confint <- as.data.frame(t(confint(mod.globfi)))
        # Attaching global.coeff.P to global.coeff
        names(global.coeff.P) <- paste( "P", names(global.coeff.P), sep = ".")
        global.coeff.wP <- cbind(global.coeff, global.coeff.P)
        
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
                                        Data.Type + Data.Type:Elevation.m.poly, 
                                      data = und.presence.SPEC, family = "binomial", 
                                      na.action = na.fail)
          mods.ALLSPEC[[levels(species.list)[S]]] <- mod.globnofi.reduced
          dredge.globnofi.reduced <- dredge(mod.globnofi.reduced, rank = AIC, subset = 
                                              dc(Elevation.m.poly, Elevation.m2.poly))
          
          # Coefficients:
          top.mods.coeff <- as.data.frame(coef(subset(dredge.globnofi.reduced, delta == 0)))
          top.mods.coeff$logLik <- 
            dredge.globnofi.reduced$logLik[dredge.globnofi.reduced$delta == 0]
          # Global coefficients:
          global.coeff <- as.data.frame(t(coef(mod.globnofi.reduced)))
          global.coeff$logLik <- as.numeric(logLik(mod.globnofi.reduced))
          global.coeff.P <- as.data.frame(t(summary(mod.globnofi.reduced)$coefficients[,4]))
          #TODO - fix confint issue!
          global.confint <- as.data.frame(t(confint(mod.globnofi.reduced)))
          # Attaching global.coeff.P to global.coeff
          names(global.coeff.P) <- paste( "P", names(global.coeff.P), sep = ".")
          global.coeff.wP <- cbind(global.coeff, global.coeff.P)
          # Correcting a naming bug - appears for poly() only
          global.coeff.wP$`Data.TypeResurvey:Elevation.m.poly` <- 
            global.coeff.wP$`Elevation.m.poly:Data.TypeResurvey`
          global.coeff.wP$`P.Data.TypeResurvey:Elevation.m.poly` <- 
            global.coeff.wP$`P.Elevation.m.poly:Data.TypeResurvey`
          
          framework.SPEC$Forced.Simpler.Mod <- paste("Yes") # Record
          
        } else {
          framework.SPEC$Discard.Later <- paste("Yes") # Record
          # Create empty DFs
          top.mods.coeff <- data.frame(logLik = rep(NA, 1))
          global.coeff.wP <- data.frame(logLik = rep(NA, 1))
          global.confint <- data.frame(Elevation.m.poly = rep(NA, 2), 
                                       row.names = c("2.5 %", "97.5 %"))
        }
        
      } else { # Run as normal
        mod.globnofi <- glm(Pres.Abs ~ Data.Type * (Elevation.m.poly + Elevation.m2.poly), 
                            data = und.presence.SPEC, family = "binomial", na.action = na.fail) 
        mods.ALLSPEC[[levels(species.list)[S]]] <- mod.globnofi
        dredge.globnofi <- dredge(mod.globnofi, rank = AIC, subset = 
                                    dc(Elevation.m.poly, Elevation.m2.poly) &&
                                    dc(Data.Type:Elevation.m.poly, Data.Type:Elevation.m2.poly))
        
        # Coefficients:
        top.mods.coeff <- as.data.frame(coef(subset(dredge.globnofi, delta == 0)))
        top.mods.coeff$logLik <- dredge.globnofi$logLik[dredge.globnofi$delta == 0]
        # Global coefficients:
        global.coeff <- as.data.frame(t(coef(mod.globnofi)))
        global.coeff$logLik <- as.numeric(logLik(mod.globnofi))
        global.coeff.P <- as.data.frame(t(summary(mod.globnofi)$coefficients[,4]))
        #TODO - fix confint issue!
        global.confint <- as.data.frame(t(confint(mod.globnofi)))
        # Attaching global.coeff.P to global.coeff
        names(global.coeff.P) <- paste( "P", names(global.coeff.P), sep = ".")
        global.coeff.wP <- cbind(global.coeff, global.coeff.P)
        
      }
      
    }
    
    #### END OF MODEL FRAMEWORK ####
    # Output: 3 DFs named top.mods.coeff, global.coeffwP
    
    
    # Null model (for Psuedo-R-squared calculation later)
    mod.NULL <- glm(Pres.Abs ~ 1, 
                    data = und.presence.SPEC, family = "binomial", na.action = na.fail)
    
    # Adding in missing coefficient headers
    for(C in 1:length(coeff.all)) {
      if(!coeff.all[C] %in% colnames(top.mods.coeff)) {
        top.mods.coeff[, coeff.all[C]] <- rep(NA, times = nrow(top.mods.coeff))
      }
    }
    # Adding missing coeff headers for global model
    for(C in 1:length(coeff.all)) {
      if(!coeff.all[C] %in% colnames(global.coeff.wP)) {
        global.coeff.wP[, coeff.all[C]] <- rep(NA, times = nrow(global.coeff.wP))
      }
      if(!coeff.all.P[C] %in% colnames(global.coeff.wP)) {
        global.coeff.wP[, coeff.all.P[C]] <- rep(NA, times = nrow(global.coeff.wP))
      }
      if(!coeff.all[C] %in% colnames(global.confint)) {
        global.confint[, coeff.all[C]] <- rep(NA, times = nrow(global.confint))
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
    
    coeff.ALLSPEC[[S]] <- ldply(Mods.list, data.frame)
    
    # Storing global model output - GIANT DATA FRAME!!
    
    global.wP <- data.frame(Species = levels(species.list)[S], 
                            # The basics:
                            Dataset = D,
                            L.Occ = sum(num.burns["1", , "Legacy"]), 
                            R.Occ = sum(num.burns["1", , "Resurvey"]), 
                            Fire.Included = ifelse(levels(species.list)[S] 
                                                   %in% species.with.fire$Species == TRUE, "Yes", "No"),
                            Type = "Global_unavg", 
                            Rsquared = 1 - global.coeff.wP$logLik / as.numeric(logLik(mod.NULL)),
                            # Coefficients:
                            Intercept = global.coeff.wP$`(Intercept)`,
                            Elevation.m = global.coeff.wP$Elevation.m.poly,
                            Elevation.m2 = global.coeff.wP$Elevation.m2.poly,
                            Data.Type.nofi = global.coeff.wP$Data.TypeResurvey,
                            Data.Type.Elevation.m.nofi = global.coeff.wP$`Data.TypeResurvey:Elevation.m.poly`,
                            Data.Type.Elevation.m2.nofi = global.coeff.wP$`Data.TypeResurvey:Elevation.m2.poly`,
                            Resurvey.Burned.fi = 
                              global.coeff.wP$New.Data.TypeResurvey.Burned,
                            Resurvey.Unburned.fi = 
                              global.coeff.wP$New.Data.TypeResurvey.Unburned,
                            Elevation.m.Res.Burn.fi = 
                              global.coeff.wP$`Elevation.m.poly:New.Data.TypeResurvey.Burned`,
                            Elevation.m.Res.Unburn.fi = 
                              global.coeff.wP$`Elevation.m.poly:New.Data.TypeResurvey.Unburned`,
                            Elevation.m2.Res.Burn.fi = 
                              global.coeff.wP$`Elevation.m2.poly:New.Data.TypeResurvey.Burned`,
                            Elevation.m2.Res.Unburn.fi = 
                              global.coeff.wP$`Elevation.m2.poly:New.Data.TypeResurvey.Unburned`,
                            # P-values:
                            P.Intercept = global.coeff.wP$`P.(Intercept)`,
                            P.Elevation.m = global.coeff.wP$P.Elevation.m.poly,
                            P.Elevation.m2 = global.coeff.wP$P.Elevation.m2.poly,
                            P.Data.Type.nofi = global.coeff.wP$P.Data.TypeResurvey,
                            P.Data.Type.Elevation.m.nofi = global.coeff.wP$`P.Data.TypeResurvey:Elevation.m.poly`,
                            P.Data.Type.Elevation.m2.nofi = global.coeff.wP$`P.Data.TypeResurvey:Elevation.m2.poly`,
                            P.Resurvey.Burned.fi = 
                              global.coeff.wP$P.New.Data.TypeResurvey.Burned,
                            P.Resurvey.Unburned.fi = 
                              global.coeff.wP$P.New.Data.TypeResurvey.Unburned,
                            P.Elevation.m.Res.Burn.fi = 
                              global.coeff.wP$`P.Elevation.m.poly:New.Data.TypeResurvey.Burned`,
                            P.Elevation.m.Res.Unburn.fi = 
                              global.coeff.wP$`P.Elevation.m.poly:New.Data.TypeResurvey.Unburned`,
                            P.Elevation.m2.Res.Burn.fi = 
                              global.coeff.wP$`P.Elevation.m2.poly:New.Data.TypeResurvey.Burned`,
                            P.Elevation.m2.Res.Unburn.fi = 
                              global.coeff.wP$`P.Elevation.m2.poly:New.Data.TypeResurvey.Unburned`,
                            # Lower CIs:
                            Intercept.CI.Lower = global.confint$`(Intercept)`[1],
                            Elevation.m.CI.Lower = global.confint$Elevation.m.poly[1],
                            Elevation.m2.CI.Lower = global.confint$Elevation.m2.poly[1],
                            Data.Type.nofi.CI.Lower = global.confint$Data.TypeResurvey[1],
                            Data.Type.Elevation.m.nofi.CI.Lower = 
                              global.confint$`Data.TypeResurvey:Elevation.m.poly`[1],
                            Data.Type.Elevation.m2.nofi.CI.Lower = 
                              global.confint$`Data.TypeResurvey:Elevation.m2.poly`[1],
                            Resurvey.Burned.fi.CI.Lower = 
                              global.confint$New.Data.TypeResurvey.Burned[1],
                            Resurvey.Unburned.fi.CI.Lower = 
                              global.confint$New.Data.TypeResurvey.Unburned[1],
                            Elevation.m.Res.Burn.fi.CI.Lower = 
                              global.confint$`Elevation.m.poly:New.Data.TypeResurvey.Burned`[1],
                            Elevation.m.Res.Unburn.fi.CI.Lower = 
                              global.confint$`Elevation.m.poly:New.Data.TypeResurvey.Unburned`[1],
                            Elevation.m2.Res.Burn.fi.CI.Lower = 
                              global.confint$`Elevation.m2.poly:New.Data.TypeResurvey.Burned`[1],
                            Elevation.m2.Res.Unburn.fi.CI.Lower = 
                              global.confint$`Elevation.m2.poly:New.Data.TypeResurvey.Unburned`[1],
                            # Upper CIs:
                            Intercept.CI.Upper = global.confint$`(Intercept)`[2],
                            Elevation.m.CI.Upper = global.confint$Elevation.m.poly[2],
                            Elevation.m2.CI.Upper = global.confint$Elevation.m2.poly[2],
                            Data.Type.nofi.CI.Upper = global.confint$Data.TypeResurvey[2],
                            Data.Type.Elevation.m.nofi.CI.Upper = 
                              global.confint$`Data.TypeResurvey:Elevation.m.poly`[2],
                            Data.Type.Elevation.m2.nofi.CI.Upper = 
                              global.confint$`Data.TypeResurvey:Elevation.m2.poly`[2],
                            Resurvey.Burned.fi.CI.Upper = 
                              global.confint$New.Data.TypeResurvey.Burned[2],
                            Resurvey.Unburned.fi.CI.Upper = 
                              global.confint$New.Data.TypeResurvey.Unburned[2],
                            Elevation.m.Res.Burn.fi.CI.Upper = 
                              global.confint$`Elevation.m.poly:New.Data.TypeResurvey.Burned`[2],
                            Elevation.m.Res.Unburn.fi.CI.Upper = 
                              global.confint$`Elevation.m.poly:New.Data.TypeResurvey.Unburned`[2],
                            Elevation.m2.Res.Burn.fi.CI.Upper = 
                              global.confint$`Elevation.m2.poly:New.Data.TypeResurvey.Burned`[2],
                            Elevation.m2.Res.Unburn.fi.CI.Upper = 
                              global.confint$`Elevation.m2.poly:New.Data.TypeResurvey.Unburned`[2],
                            row.names = NULL)
    
    global.wP.ALLSPEC[[S]] <- global.wP
    
    ## Storing record of framework
    
    framework.ALLSPEC[[S]] <- framework.SPEC
    
  }  
  
  #Collapse ALLSPEC lists into single dataframe, store that df as part of ALLDAT list
  
  coeff.ALLDAT[[D]] <- ldply(coeff.ALLSPEC, data.frame)
  global.wP.ALLDAT[[D]] <- ldply(global.wP.ALLSPEC, data.frame)
  framework.ALLDAT[[D]] <- ldply(framework.ALLSPEC, data.frame)
  # save(mods.ALLSPEC, file = paste("Model_data/Model.Lists", D, "Rda", sep = "."))
  
  #### END OF SPECIES LOOP
  
  if(D == 100) system("say Your loop is done")
  
}

#### END OF DATASET LOOP

# Collate ALLDAT lists into one big DF

coeff.ALLDAT.allsets <- ldply(coeff.ALLDAT, data.frame)
framework.ALLDAT.allsets <- ldply(framework.ALLDAT, data.frame)
global.wP.ALLDAT.allsets <- ldply(global.wP.ALLDAT, data.frame)


#### Removing error-laden datasets ####

# Create new column for joining
framework.ALLDAT.allsets$Join <- 
  paste(framework.ALLDAT.allsets$Dataset, framework.ALLDAT.allsets$Species, sep = ".")
discard.me <- framework.ALLDAT.allsets[ , c(6, 8)]
coeff.ALLDAT.allsets$Join <- 
  paste(coeff.ALLDAT.allsets$Dataset, coeff.ALLDAT.allsets$Species, sep = ".")
global.wP.ALLDAT.allsets$Join <- 
  paste(global.wP.ALLDAT.allsets$Dataset, global.wP.ALLDAT.allsets$Species, sep = ".")

# Merge framework DF with coeff and confint DF
coeff.ALLDAT.allsets.join <- merge(coeff.ALLDAT.allsets, discard.me, 
                                   by = "Join")
global.wP.ALLDAT.allsets.join <- merge(global.wP.ALLDAT.allsets, discard.me, 
                                       by = "Join")

# Remove datasets associated with errors (Discard.Later == Yes)
coeff.ALLDAT.finaldf.big <- 
  coeff.ALLDAT.allsets.join[!coeff.ALLDAT.allsets.join$Discard.Later == "Yes", ]
global.wP.ALLDAT.finaldf.big <- 
  global.wP.ALLDAT.allsets.join[!global.wP.ALLDAT.allsets.join$Discard.Later == "Yes", ]
coeff.ALLDAT.finaldf <- coeff.ALLDAT.finaldf.big[, c(2:20)]
global.wP.ALLDAT.finaldf <- global.wP.ALLDAT.finaldf.big[, c(2:56)]

# Separate into yes / no fire (currently just for global model output)
global.wP.ALLDAT.finaldf.fire <- 
  global.wP.ALLDAT.finaldf[global.wP.ALLDAT.finaldf$Fire.Included == "Yes", 
                           c(1:10, 14:19, 20:22, 26:31, 32:34, 38:43, 44:46, 50:55)]
global.wP.ALLDAT.finaldf.nofire <- 
  global.wP.ALLDAT.finaldf[global.wP.ALLDAT.finaldf$Fire.Included == "No", 
                           c(1:10, 11:13, 20:22, 23:25, 32:34, 35:37, 44:46, 47:49)]


# Store output as CSV

write.csv(coeff.ALLDAT.finaldf, 
          file = "data/3c_top_mod_coefficients.csv", 
          row.names = FALSE)

write.csv(global.wP.ALLDAT.finaldf.fire, 
          file = "data/3c_ORTHO_global_mod_coefficients_with_P_FIRE.csv", 
          row.names = FALSE)

write.csv(global.wP.ALLDAT.finaldf.nofire, 
          file = "data/3c_ORTHO_global_mod_coefficients_with_P_NOFIRE.csv", 
          row.names = FALSE)

# Adding CSV of transformed vs untransformed variables - should be the same between species
transformed.df <- data.frame(Elevation.m = und.presence.SPEC$Elevation.m, 
                             Elevation.m.poly = und.presence.SPEC$Elevation.m.poly,
                             Elevation.m2 = und.presence.SPEC$Elevation.m2,
                             Elevation.m2.poly = und.presence.SPEC$Elevation.m2.poly,
                             row.names = NULL)

write.csv(transformed.df, 
         file = "data/3c_transformed_polynomials.csv",
        row.names = FALSE)







