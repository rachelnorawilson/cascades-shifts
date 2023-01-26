
# Created: Dec. 15, 2021 from script 8
# Updated: May 11, 2022

# IMPORTANT NOTE: unless otherwise indicated, always use Understory_All.csv for these analyses as it is the ONLY file with up-to-date corrections.

# Packages needed:
library(tidyverse)

## Step 1: Load coefficient data

coeff.ALLDAT <- read.csv("data/3c_new_coefficients.csv", header = TRUE) %>% 
  select(-Dataset, -L.Occ, -R.Occ, -deltaAIC, -Weight, -Rsquared)
coeff.ALLDAT[coeff.ALLDAT$Species == "EPAN", 1] <- paste("CHAN") # Correcting taxonomy issue
coeff.ALLDAT[is.na(coeff.ALLDAT)] <- 0 # Hard code absent coeffs as 0 before averaging
coeff.ALLDAT$Elevation.m2.Res.Burn.fi[coeff.ALLDAT$Species == "VAME"] <- NA
coeff.ALLDAT$Elevation.m2.Res.Unburn.fi[coeff.ALLDAT$Species == "VAME"] <- NA
coeff.ALLDAT$Data.Type.Elevation.m2.nofi[coeff.ALLDAT$Species == "HODI"] <- NA #remove terms that weren't fit

coeff.fire <- #split by fire vs no-fire
  coeff.ALLDAT[coeff.ALLDAT$Fire.Included == "Yes" & coeff.ALLDAT$Type == "Avg", c(1, 5:6, 10:15)] 

coeff.nofire <- coeff.ALLDAT[coeff.ALLDAT$Fire.Included == "No" & coeff.ALLDAT$Type == "Avg", c(1, 5:9)] 

## Step 2: Summarize mean, lower and upper CI of each coefficient across the rarefactions (FIRE SPECIES)
means.fire <- coeff.fire %>% 
  group_by(Species) %>% 
  summarise(across(.cols = Elevation.m:Elevation.m2.Res.Unburn.fi, ~ unname(quantile(.x, 0.5, na.rm=TRUE)))) %>% #, .names = "mean_{.col}" 
  mutate(param="mean")
lowers.fire <- coeff.fire %>% 
  group_by(Species) %>% 
  summarise(across(.cols = Elevation.m:Elevation.m2.Res.Unburn.fi, ~ unname(quantile(.x, 0.025, na.rm=TRUE)))) %>% #, .names="lower_{.col}" 
  mutate(param="lower")
uppers.fire <- coeff.fire %>% 
  group_by(Species) %>% 
  summarise(across(.cols = Elevation.m:Elevation.m2.Res.Unburn.fi, ~ unname(quantile(.x, 0.975, na.rm=TRUE)))) %>% #, .names="upper_{.col}" 
  mutate(param="upper")

## Step 3: reshape for graphing
all.fire <- rbind(means.fire, lowers.fire, uppers.fire) %>% 
  pivot_longer(!c(Species, param), names_to="Parameter", values_to="Estimate") %>% 
  pivot_wider(names_from=param, values_from="Estimate")

## Step 4: Classic forest plots for fire species

# order of parameters along y axis
order.list.fire <- c("Elevation.m", 
                     "Elevation.m2", 
                     "Resurvey.Burned.fi", 
                     "Resurvey.Unburned.fi",
                     "Elevation.m.Res.Burn.fi", 
                     "Elevation.m.Res.Unburn.fi", 
                     "Elevation.m2.Res.Burn.fi", 
                     "Elevation.m2.Res.Unburn.fi")

all.fire$Parameter <- factor(all.fire$Parameter, levels = rev(order.list.fire))

# tick labels for y axis
vars.fire <- c("Elevation", 
               expression("Elevation" ^ 2), 
               "2015, Burned",
               "2015, Unburned",
               "Elevation * Burned",
               "Elevation * Unburned",
               expression("Elevation" ^ 2 * " * Burned"), 
               expression("Elevation" ^ 2 * " * Unburned"))

all.fire <- all.fire %>% 
  mutate(lowernonzero = ifelse(mean>0 & lower>0, "y", "n")) %>% 
  mutate(uppernonzero = ifelse(mean<0 & upper<0, "y", "n")) %>% 
  mutate(nonzero = ifelse(mean>0, lowernonzero, uppernonzero))

# faceted plot
forestplot.fire <- ggplot(dat=all.fire, aes(y=Parameter, x=mean, xmin=lower, xmax=upper, color=nonzero)) +
  facet_wrap(~Species) +
  geom_point(cex=3) + 
  scale_color_manual(values=c("grey", "black")) +
  geom_errorbarh(height=0.3) + 
  geom_vline(xintercept=0, linetype="dotted") +
  scale_y_discrete(labels=rev(vars.fire)) +
  xlab("Estimate") +
  theme_classic() +
  theme(strip.background = element_rect(fill = "lightgrey", colour = "black", size = 1)) + 
  theme(legend.position = "none")

ggsave("figures/forestplot_coeffs_fire.pdf", forestplot.fire, width=12, height=8)


## Step 5: Repeat steps 2-4 for no-fire species

means.nofire <- coeff.nofire %>% 
  group_by(Species) %>% 
  summarise(across(.cols = Elevation.m:Data.Type.Elevation.m2.nofi, ~ unname(quantile(.x, 0.5, na.rm=TRUE)))) %>%  
  mutate(param="mean")
lowers.nofire <- coeff.nofire %>% 
  group_by(Species) %>% 
  summarise(across(.cols = Elevation.m:Data.Type.Elevation.m2.nofi, ~ unname(quantile(.x, 0.025, na.rm=TRUE)))) %>%  
  mutate(param="lower")
uppers.nofire <- coeff.nofire %>% 
  group_by(Species) %>% 
  summarise(across(.cols = Elevation.m:Data.Type.Elevation.m2.nofi, ~ unname(quantile(.x, 0.975, na.rm=TRUE)))) %>%  
  mutate(param="upper")

all.nofire <- rbind(means.nofire, lowers.nofire, uppers.nofire) %>% 
  pivot_longer(!c(Species, param), names_to="Parameter", values_to="Estimate") %>% 
  pivot_wider(names_from=param, values_from="Estimate")

order.list.nofire <- c("Elevation.m", 
                       "Elevation.m2", 
                       "Data.Type.nofi",
                       "Data.Type.Elevation.m.nofi", 
                       "Data.Type.Elevation.m2.nofi")

all.nofire$Parameter <- factor(all.nofire$Parameter, levels = rev(order.list.nofire))

all.nofire <- all.nofire %>% 
  mutate(lowernonzero = ifelse(mean>0 & lower>0, "y", "n")) %>% 
  mutate(uppernonzero = ifelse(mean<0 & upper<0, "y", "n")) %>% 
  mutate(nonzero = ifelse(mean>0, lowernonzero, uppernonzero))

# tick labels for y axis
vars.nofire <- c("Elevation", 
                 expression("Elevation" ^ 2), 
                 "Time",
                 "Elevation * Time",
                 expression("Elevation" ^ 2 * " * Time"))

# faceted plot
forestplot.nofire <- ggplot(dat=all.nofire, aes(y=Parameter, x=mean, xmin=lower, xmax=upper, color=nonzero)) +
  facet_wrap(~Species) +
  geom_point(cex=2) + 
  scale_color_manual(values=c("grey", "black")) +
  geom_errorbarh(height=0.3) + 
  geom_vline(xintercept=0, linetype="dotted") +
  scale_y_discrete(labels=rev(vars.nofire)) +
  xlab("Estimate") +
  theme_classic() +
  theme(strip.background = element_rect(fill = "lightgrey", colour = "black", size = 1)) +
  theme(legend.position = "none")

ggsave("figures/forestplot_coeffs_nofire.pdf", forestplot.nofire, width=12, height=8)



