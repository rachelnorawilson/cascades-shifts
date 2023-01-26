# Created: Feb. 16, 2021
# Modified: Mar. 10, 2022

#### This script visualizes raw (arithmetic) patterns in rarefied data; see script 6 for visualizing model outputs

### Load libraries
library(tidyverse)
library(Hmisc) #for mean_sdl function to put mean+/-SD on violin plots
library(cowplot) #for multipanel figures
library(viridis) #for color-blind-friendly palettes

### Load data
# rarefied data (as list of dataframes named rare.ALL):
load("data/rare.ALL.Rda")

### Species list

# start with all species
load("data/Species.List.Rda")
species.list <- shifts %>% 
  filter(Species.Code!="MOSS") %>% 
  select(Species=Species.Code)
species.list$Species <- as.character(species.list$Species)

# separate out fire species
species.list.fire <- read.csv("data/3c_new_coefficients.csv", header = TRUE) %>% 
  filter(Type=="Avg")  %>% 
  filter(Fire.Included=="Yes") %>% 
  group_by(Species) %>% 
  summarise(Species=first(Species))

species.list.nofire <- anti_join(species.list, species.list.fire)

### Set up empty matrices to store values
el.mins.raw.leg.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.mins.025.leg.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.maxs.raw.leg.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.maxs.975.leg.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.meds.leg.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.mins.raw.res.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.mins.025.res.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.maxs.raw.res.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.maxs.975.res.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.meds.res.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.highshift.95.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])
el.lowshift.95.nofire <- matrix(nrow=100,ncol=dim(species.list.nofire)[1])

el.mins.raw.leg.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.mins.025.leg.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.maxs.raw.leg.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.maxs.975.leg.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.meds.leg.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.mins.raw.res.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.mins.025.res.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.maxs.raw.res.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.maxs.975.res.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.meds.res.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.highshift.95.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])
el.lowshift.95.fire <- matrix(nrow=100,ncol=dim(species.list.fire)[1])

### Clunky for-loops across rarefied datasets by species

## no-fire species
for(D in 1:100) { # Run time: 40 seconds
  und <- rare.ALL[[D]]
  und$Data.Type <- as.factor(und$Data.Type)
  und$Elevation.m <- as.numeric(und$Elevation.m)
  for(S in 1:dim(species.list.nofire)[1]) {
    species <- as.list(species.list.nofire[S,1])
    und.SPEC <- und %>% 
      filter(Species.Code == species) %>% 
      droplevels()
    und.presence.SPEC.leg <- und.SPEC %>% 
      filter(Data.Type=="Legacy" & Pres.Abs==1) %>% 
      mutate(el.min.raw.leg.nofire = min(Elevation.m),
             el.max.raw.leg.nofire = max(Elevation.m),
             el.med.leg.nofire = median(Elevation.m),
             el.min.025.leg.nofire = quantile(Elevation.m, probs=0.025),
             el.max.975.leg.nofire = quantile(Elevation.m, probs=0.975))
    und.presence.SPEC.res <- und.SPEC %>% 
      filter(Data.Type=="Resurvey" & Pres.Abs==1) %>% 
      mutate(el.min.raw.res.nofire = min(Elevation.m),
             el.max.raw.res.nofire = max(Elevation.m),
             el.med.res.nofire = median(Elevation.m),             
             el.min.025.res.nofire = quantile(Elevation.m, probs=0.025),
             el.max.975.res.nofire = quantile(Elevation.m, probs=0.975))
    el.mins.raw.leg.nofire[D,S] <- und.presence.SPEC.leg$el.min.raw.leg.nofire[1]
    el.mins.025.leg.nofire[D,S] <- und.presence.SPEC.leg$el.min.025.leg.nofire[1]
    el.maxs.raw.leg.nofire[D,S] <- und.presence.SPEC.leg$el.max.raw.leg.nofire[1]
    el.maxs.975.leg.nofire[D,S] <- und.presence.SPEC.leg$el.max.975.leg.nofire[1]
    el.meds.leg.nofire[D,S] <- und.presence.SPEC.leg$el.med.leg.nofire[1]
    el.mins.raw.res.nofire[D,S] <- und.presence.SPEC.res$el.min.raw.res.nofire[1]
    el.mins.025.res.nofire[D,S] <- und.presence.SPEC.res$el.min.025.res.nofire[1]
    el.maxs.raw.res.nofire[D,S] <- und.presence.SPEC.res$el.max.raw.res.nofire[1]
    el.maxs.975.res.nofire[D,S] <- und.presence.SPEC.res$el.max.975.res.nofire[1]
    el.meds.res.nofire[D,S] <- und.presence.SPEC.res$el.med.res.nofire[1]
    el.highshift.95.nofire[D,S] <- und.presence.SPEC.res$el.max.975.res.nofire[1] - und.presence.SPEC.leg$el.max.975.leg.nofire[1]
    el.lowshift.95.nofire[D,S] <- und.presence.SPEC.res$el.min.025.res.nofire[1] - und.presence.SPEC.leg$el.min.025.leg.nofire[1]
  }
}

## fire species
for(D in 1:100) { #Run time: 10 seconds
  und <- rare.ALL[[D]]
  und$Data.Type <- as.factor(und$Data.Type)
  und$Elevation.m <- as.numeric(und$Elevation.m)
  for(S in 1:dim(species.list.fire)[1]) {
    species <- as.list(species.list.fire[S,1])
    und.SPEC <- und %>% 
      filter(Species.Code == species) %>% 
      droplevels()
    und.presence.SPEC.leg <- und.SPEC %>% 
      filter(Data.Type=="Legacy" & Pres.Abs==1) %>% 
      mutate(el.min.raw.leg.fire = min(Elevation.m),
             el.max.raw.leg.fire = max(Elevation.m),
             el.med.leg.fire = median(Elevation.m),
             el.min.025.leg.fire = quantile(Elevation.m, probs=0.025),
             el.max.975.leg.fire = quantile(Elevation.m, probs=0.975))
    und.presence.SPEC.res <- und.SPEC %>% 
      filter(Data.Type=="Resurvey" & Pres.Abs==1) %>% 
      #filter(Fires=="Unburned") %>% #temporary toggle %>% 
      mutate(el.min.raw.res.fire = min(Elevation.m),
             el.max.raw.res.fire = max(Elevation.m),
             el.med.res.fire = median(Elevation.m),
             el.min.025.res.fire = quantile(Elevation.m, probs=0.025),
             el.max.975.res.fire = quantile(Elevation.m, probs=0.975))
    el.mins.raw.leg.fire[D,S] <- und.presence.SPEC.leg$el.min.raw.leg.fire[1]
    el.mins.025.leg.fire[D,S] <- und.presence.SPEC.leg$el.min.025.leg.fire[1]
    el.maxs.raw.leg.fire[D,S] <- und.presence.SPEC.leg$el.max.raw.leg.fire[1]
    el.maxs.975.leg.fire[D,S] <- und.presence.SPEC.leg$el.max.975.leg.fire[1]
    el.meds.leg.fire[D,S] <- und.presence.SPEC.leg$el.med.leg.fire[1]
    el.mins.raw.res.fire[D,S] <- und.presence.SPEC.res$el.min.raw.res.fire[1]
    el.mins.025.res.fire[D,S] <- und.presence.SPEC.res$el.min.025.res.fire[1]
    el.maxs.raw.res.fire[D,S] <- und.presence.SPEC.res$el.max.raw.res.fire[1]
    el.maxs.975.res.fire[D,S] <- und.presence.SPEC.res$el.max.975.res.fire[1]
    el.meds.res.fire[D,S] <- und.presence.SPEC.res$el.med.res.fire[1]
    el.highshift.95.fire[D,S] <- und.presence.SPEC.res$el.max.975.res.fire[1] - und.presence.SPEC.leg$el.max.975.leg.fire[1]
    el.lowshift.95.fire[D,S] <- und.presence.SPEC.res$el.min.025.res.fire[1] - und.presence.SPEC.leg$el.min.025.leg.fire[1]
  }
}

## collapse to averages across rarefied replicates

# no-fire species
el.mins.raw.leg.means.nofire <- as.data.frame(el.mins.raw.leg.nofire) %>% summarise(across(starts_with("V"),mean)) 
el.mins.025.leg.means.nofire <- as.data.frame(el.mins.025.leg.nofire) %>% summarise(across(starts_with("V"),mean))
el.mins.raw.res.means.nofire <- as.data.frame(el.mins.raw.res.nofire) %>% summarise(across(starts_with("V"),mean))
el.mins.025.res.means.nofire <- as.data.frame(el.mins.025.res.nofire) %>% summarise(across(starts_with("V"),mean))
el.maxs.raw.leg.means.nofire <- as.data.frame(el.maxs.raw.leg.nofire) %>% summarise(across(starts_with("V"),mean))
el.maxs.975.leg.means.nofire <- as.data.frame(el.maxs.975.leg.nofire) %>% summarise(across(starts_with("V"),mean))
el.maxs.raw.res.means.nofire <- as.data.frame(el.maxs.raw.res.nofire) %>% summarise(across(starts_with("V"),mean))
el.maxs.975.res.means.nofire <- as.data.frame(el.maxs.975.res.nofire) %>% summarise(across(starts_with("V"),mean))
el.meds.leg.means.nofire <- as.data.frame(el.meds.leg.nofire) %>% summarise(across(starts_with("V"),mean))
el.meds.res.means.nofire <- as.data.frame(el.meds.res.nofire) %>% summarise(across(starts_with("V"),mean))
el.highshift.means.nofire <- as.data.frame(el.highshift.95.nofire) %>% summarise(across(starts_with("V"),mean))
el.highshift.sd.nofire <- as.data.frame(el.highshift.95.nofire) %>% summarise(across(starts_with("V"),sd))
el.lowshift.means.nofire <- as.data.frame(el.lowshift.95.nofire) %>% summarise(across(starts_with("V"),mean))
el.lowshift.sd.nofire <- as.data.frame(el.lowshift.95.nofire) %>% summarise(across(starts_with("V"),sd))

# fire species
el.mins.raw.leg.means.fire <- as.data.frame(el.mins.raw.leg.fire) %>% summarise(across(starts_with("V"),mean))
el.mins.025.leg.means.fire <- as.data.frame(el.mins.025.leg.fire) %>% summarise(across(starts_with("V"),mean))
el.mins.raw.res.means.fire <- as.data.frame(el.mins.raw.res.fire) %>% summarise(across(starts_with("V"),mean))
el.mins.025.res.means.fire <- as.data.frame(el.mins.025.res.fire) %>% summarise(across(starts_with("V"),mean))
el.maxs.raw.leg.means.fire <- as.data.frame(el.maxs.raw.leg.fire) %>% summarise(across(starts_with("V"),mean))
el.maxs.975.leg.means.fire <- as.data.frame(el.maxs.975.leg.fire) %>% summarise(across(starts_with("V"),mean))
el.maxs.raw.res.means.fire <- as.data.frame(el.maxs.raw.res.fire) %>% summarise(across(starts_with("V"),mean))
el.maxs.975.res.means.fire <- as.data.frame(el.maxs.975.res.fire) %>% summarise(across(starts_with("V"),mean))
el.meds.leg.means.fire <- as.data.frame(el.meds.leg.fire) %>% summarise(across(starts_with("V"),mean))
el.meds.res.means.fire <- as.data.frame(el.meds.res.fire) %>% summarise(across(starts_with("V"),mean))
el.highshift.means.fire <- as.data.frame(el.highshift.95.fire) %>% summarise(across(starts_with("V"),mean))
el.highshift.sd.fire <- as.data.frame(el.highshift.95.fire) %>% summarise(across(starts_with("V"),sd))
el.lowshift.means.fire <- as.data.frame(el.lowshift.95.fire) %>% summarise(across(starts_with("V"),mean))
el.lowshift.sd.fire <- as.data.frame(el.lowshift.95.fire) %>% summarise(across(starts_with("V"),sd))

## reshape and join into one frame

# no-fire species
el.mins.raw.leg.tall.nofire <- gather(el.mins.raw.leg.means.nofire, "species", "min.raw.leg", 1:dim(species.list.nofire)[1])
el.mins.025.leg.tall.nofire <- gather(el.mins.025.leg.means.nofire, "species", "min.025.leg", 1:dim(species.list.nofire)[1])
el.mins.raw.res.tall.nofire <- gather(el.mins.raw.res.means.nofire, "species", "min.raw.res", 1:dim(species.list.nofire)[1])
el.mins.025.res.tall.nofire <- gather(el.mins.025.res.means.nofire, "species", "min.025.res", 1:dim(species.list.nofire)[1])
el.maxs.raw.leg.tall.nofire <- gather(el.maxs.raw.leg.means.nofire, "species", "max.raw.leg", 1:dim(species.list.nofire)[1])
el.maxs.975.leg.tall.nofire <- gather(el.maxs.975.leg.means.nofire, "species", "max.975.leg", 1:dim(species.list.nofire)[1])
el.maxs.raw.res.tall.nofire <- gather(el.maxs.raw.res.means.nofire, "species", "max.raw.res", 1:dim(species.list.nofire)[1])
el.maxs.975.res.tall.nofire <- gather(el.maxs.975.res.means.nofire, "species", "max.975.res", 1:dim(species.list.nofire)[1])
el.meds.leg.tall.nofire <- gather(el.meds.leg.means.nofire, "species", "med.leg", 1:dim(species.list.nofire)[1])
el.meds.res.tall.nofire <- gather(el.meds.res.means.nofire, "species", "med.res", 1:dim(species.list.nofire)[1])
el.mean.highshift.95.tall.nofire <- gather(el.highshift.means.nofire, "species", "mean.high", 1:dim(species.list.nofire)[1])
el.sd.highshift.95.tall.nofire <- gather(el.highshift.sd.nofire, "species", "sd.high", 1:dim(species.list.nofire)[1])
el.mean.lowshift.95.tall.nofire <- gather(el.lowshift.means.nofire, "species", "mean.low", 1:dim(species.list.nofire)[1])
el.sd.lowshift.95.tall.nofire <- gather(el.lowshift.sd.nofire, "species", "sd.low", 1:dim(species.list.nofire)[1])

rarefied.change.nofire <- left_join(left_join(left_join(left_join(left_join(left_join(left_join(left_join(left_join(el.mins.raw.leg.tall.nofire, el.mins.raw.res.tall.nofire),el.maxs.raw.leg.tall.nofire), el.maxs.raw.res.tall.nofire), el.meds.leg.tall.nofire), el.meds.res.tall.nofire), el.mins.025.leg.tall.nofire), el.mins.025.res.tall.nofire), el.maxs.975.leg.tall.nofire), el.maxs.975.res.tall.nofire)

rarefied.change.nofire.update <- left_join(left_join(left_join(el.mean.highshift.95.tall.nofire, el.sd.highshift.95.tall.nofire),el.mean.lowshift.95.tall.nofire), el.sd.lowshift.95.tall.nofire)

rarefied.change.nofire <- cbind(rarefied.change.nofire, species.list.nofire)
rarefied.change.nofire$fire <- "no"

rarefied.change.nofire.update <- cbind(rarefied.change.nofire.update, species.list.nofire)
rarefied.change.nofire.update$fire <- "no"

# fire species
el.mins.raw.leg.tall.fire <- gather(el.mins.raw.leg.means.fire, "species", "min.raw.leg", 1:dim(species.list.fire)[1])
el.mins.025.leg.tall.fire <- gather(el.mins.025.leg.means.fire, "species", "min.025.leg", 1:dim(species.list.fire)[1])
el.mins.raw.res.tall.fire <- gather(el.mins.raw.res.means.fire, "species", "min.raw.res", 1:dim(species.list.fire)[1])
el.mins.025.res.tall.fire <- gather(el.mins.025.res.means.fire, "species", "min.025.res", 1:dim(species.list.fire)[1])
el.maxs.raw.leg.tall.fire <- gather(el.maxs.raw.leg.means.fire, "species", "max.raw.leg", 1:dim(species.list.fire)[1])
el.maxs.975.leg.tall.fire <- gather(el.maxs.975.leg.means.fire, "species", "max.975.leg", 1:dim(species.list.fire)[1])
el.maxs.raw.res.tall.fire <- gather(el.maxs.raw.res.means.fire, "species", "max.raw.res", 1:dim(species.list.fire)[1])
el.maxs.975.res.tall.fire <- gather(el.maxs.975.res.means.fire, "species", "max.975.res", 1:dim(species.list.fire)[1])
el.meds.leg.tall.fire <- gather(el.meds.leg.means.fire, "species", "med.leg", 1:dim(species.list.fire)[1])
el.meds.res.tall.fire <- gather(el.meds.res.means.fire, "species", "med.res", 1:dim(species.list.fire)[1])
el.mean.highshift.95.tall.fire <- gather(el.highshift.means.fire, "species", "mean.high", 1:dim(species.list.fire)[1])
el.sd.highshift.95.tall.fire <- gather(el.highshift.sd.fire, "species", "sd.high", 1:dim(species.list.fire)[1])
el.mean.lowshift.95.tall.fire <- gather(el.lowshift.means.fire, "species", "mean.low", 1:dim(species.list.fire)[1])
el.sd.lowshift.95.tall.fire <- gather(el.lowshift.sd.fire, "species", "sd.low", 1:dim(species.list.fire)[1])

rarefied.change.fire <- left_join(left_join(left_join(left_join(left_join(left_join(left_join(left_join(left_join(el.mins.raw.leg.tall.fire, el.mins.raw.res.tall.fire),el.maxs.raw.leg.tall.fire), el.maxs.raw.res.tall.fire), el.meds.leg.tall.fire), el.meds.res.tall.fire), el.mins.025.leg.tall.fire), el.mins.025.res.tall.fire), el.maxs.975.leg.tall.fire), el.maxs.975.res.tall.fire)


rarefied.change.fire.update <- left_join(left_join(left_join(el.mean.highshift.95.tall.fire, el.sd.highshift.95.tall.fire),el.mean.lowshift.95.tall.fire), el.sd.lowshift.95.tall.fire)

rarefied.change.fire <- cbind(rarefied.change.fire, species.list.fire)
rarefied.change.fire$fire <- "yes"

rarefied.change.fire.update <- cbind(rarefied.change.fire.update, species.list.fire)
rarefied.change.fire.update$fire <- "yes"

## calculate range changes
# no-fires species
rarefied.change.nofire <- rarefied.change.nofire %>% 
  mutate(rear.change.raw = min.raw.res - min.raw.leg,
         rear.change.perc = min.025.res - min.025.leg,
         med.change = med.res - med.leg,
         lead.change.raw = max.raw.res - max.raw.leg,
         lead.change.perc = max.975.res - max.975.leg)#,
#span.change = (max.res-min.res) - (max.leg - min.leg))

# fire species
rarefied.change.fire <- rarefied.change.fire %>% 
  mutate(rear.change.raw = min.raw.res - min.raw.leg,
         rear.change.perc = min.025.res - min.025.leg,
         med.change = med.res - med.leg,
         lead.change.raw = max.raw.res - max.raw.leg,
         lead.change.perc = max.975.res - max.975.leg)#,
#span.change = (max.res-min.res) - (max.leg - min.leg))

rarefied.change.calcs <- rbind(rarefied.change.nofire, rarefied.change.fire)
write.csv(rarefied.change.calcs, "data/5_range.change.calcs.csv")

rarefied.change.calcs.update <- rbind(rarefied.change.nofire.update, rarefied.change.fire.update)
write.csv(rarefied.change.calcs.update, "data/5_range.change.calcs.update.csv")

### Statistical tests for differences between fire and no-fire species groups
rarefied.change.calcs <- read_csv("data/5_range.change.calcs.csv")

rear.t <- t.test(rarefied.change.calcs$rear.change.perc[rarefied.change.calcs$fire=="no"], rarefied.change.calcs$rear.change.perc[rarefied.change.calcs$fire=="yes"])
rear.t

rear.t.raw <- t.test(rarefied.change.calcs$rear.change.raw[rarefied.change.calcs$fire=="no"], rarefied.change.calcs$rear.change.raw[rarefied.change.calcs$fire=="yes"])
rear.t.raw

med.t <- t.test(rarefied.change.calcs$med.change[rarefied.change.calcs$fire=="no"], rarefied.change.calcs$med.change[rarefied.change.calcs$fire=="yes"])
med.t

lead.t <- t.test(rarefied.change.calcs$lead.change.perc[rarefied.change.calcs$fire=="no"], rarefied.change.calcs$lead.change.perc[rarefied.change.calcs$fire=="yes"])
lead.t

lead.t.raw <- t.test(rarefied.change.calcs$lead.change.raw[rarefied.change.calcs$fire=="no"], rarefied.change.calcs$lead.change.raw[rarefied.change.calcs$fire=="yes"])
lead.t.raw

### Reshape for violin plotting

# no-fire species
rarefied.change.tall.nofire.raw <- rarefied.change.nofire %>% 
  select(rear.change.raw, med.change, lead.change.raw) %>% 
  gather("edge", "change", 1:3) %>% 
  add_column(fire="no")

rarefied.change.tall.nofire.perc <- rarefied.change.nofire %>% 
  select(rear.change.perc, med.change, lead.change.perc) %>% 
  gather("edge", "change", 1:3) %>% 
  add_column(fire="no")

# fire species
rarefied.change.tall.fire.raw <- rarefied.change.fire %>% 
  select(rear.change.raw, med.change, lead.change.raw) %>% 
  gather("edge", "change", 1:3) %>% 
  add_column(fire="yes")

rarefied.change.tall.fire.perc <- rarefied.change.fire %>% 
  select(rear.change.perc, med.change, lead.change.perc) %>% 
  gather("edge", "change", 1:3) %>% 
  add_column(fire="yes")

rarefied.change.tall.perc <- bind_rows(rarefied.change.tall.nofire.perc, rarefied.change.tall.fire.perc)
rarefied.change.tall.raw <- bind_rows(rarefied.change.tall.nofire.raw, rarefied.change.tall.fire.raw)

level_order.raw = c("rear.change.raw", "med.change", "lead.change.raw")
level_order.perc = c("rear.change.perc", "med.change", "lead.change.perc")

# Violins by range position

violin.plot.raw <- ggplot(rarefied.change.tall.raw, aes(x=factor(edge, level=level_order.raw), y=change, fill=fire)) + 
  geom_violin() +
  scale_fill_viridis(discrete=TRUE, alpha=0.7) +
  #stat_summary(fun=mean, geom="point", cex=2)  +
  theme_classic() +
  geom_hline(yintercept=0, lty="dashed") +
  xlab("") +
  scale_x_discrete(labels=c("Lower\nedge", "Range\ncenter", "Upper\nedge")) +
  ylim(-600,700) +
  ylab(c("Elevational change (m)\n1983-2015")) +
  geom_segment(aes(x = 0.8, xend = 1.2, y = 400, yend = 400)) +
  geom_segment(aes(x = 1.8, xend = 2.2, y = 400, yend = 400)) +
  geom_segment(aes(x = 2.75, xend = 3.15, y = 400, yend = 400)) +
  annotate("text", x=1, y=450, label="ns") + #from t-tests above
  annotate("text", x=2, y=450, label="+") +
  annotate("text", x=2.95, y=440, label="*") 
violin.plot.raw

ggsave("figures/violin_1panel_raw.pdf", violin.plot.raw, device="pdf", width=8, height=5)

violin.plot.perc <- ggplot(rarefied.change.tall.perc, aes(x=factor(edge, level=level_order.perc), y=change, fill=fire)) +#, color=edge, fill=edge)) + 
  geom_violin() +
  geom_point(aes(color=fire), position=position_jitter(seed=1, width=0.2)) +
  scale_fill_viridis(discrete=TRUE, alpha=0.7) +
  #stat_summary(fun=mean, geom="point", cex=2)  +
  theme_classic() +
  geom_hline(yintercept=0, lty="dashed") +
  xlab("") +
  scale_x_discrete(labels=c("Lower\nedge", "Range\ncenter", "Upper\nedge")) +
  theme(axis.text.x=element_text(size=12)) +
  ylim(-600,700) +
  ylab(c("Elevational change (m)\n1983-2015")) +
  geom_segment(aes(x = 0.8, xend = 1.2, y = 400, yend = 400)) +
  geom_segment(aes(x = 1.8, xend = 2.2, y = 400, yend = 400)) +
  geom_segment(aes(x = 2.75, xend = 3.15, y = 400, yend = 400)) +
  annotate("text", x=1, y=450, label="ns") + #from t-tests above
  annotate("text", x=2, y=450, label="+") +
  annotate("text", x=2.95, y=440, label="*") 
violin.plot.perc

ggsave("figures/violin_1panel_perc.pdf", violin.plot.perc, device="pdf", width=8, height=5)



## Freeman-style elevation ranges

rarefied.change.calcs <- read_csv("data/5_range.change.calcs.csv")

# fire and no-fire species separated
rarefied.change.calcs.fire <- rarefied.change.calcs %>% 
  filter(fire=="yes") %>% 
  mutate(species.rank.med = dense_rank(med.leg)+35, #silly workaround to get unique 4-letter species codes on fire species
         #species.rank.min = dense_rank(min.raw.leg), #doesn't work because of ties
         #species.rank.max = dense_rank(max.raw.leg), #doesn't work because of ties
         both.min.raw = pmax(min.raw.leg, min.raw.res),
         both.max.raw = pmin(max.raw.leg, max.raw.res),
         both.min.perc = pmax(min.025.leg, min.025.res),
         both.max.perc = pmin(max.975.leg, max.975.res))

rarefied.change.calcs.nofire <- rarefied.change.calcs %>% 
  filter(fire=="no") %>% 
  mutate(species.rank.med = dense_rank(med.leg),
         #species.rank.min = dense_rank(min.raw.leg), #doesn't work because of ties
         #species.rank.max = dense_rank(max.raw.leg), #doesn't work because of ties
         both.min.raw = pmax(min.raw.leg, min.raw.res),
         both.max.raw = pmin(max.raw.leg, max.raw.res),
         both.min.perc = pmax(min.025.leg, min.025.res),
         both.max.perc = pmin(max.975.leg, max.975.res))

rarefied.change.calcs.facet <- rbind(rarefied.change.calcs.nofire, rarefied.change.calcs.fire)

species.labels <- rarefied.change.calcs.facet$Species[order(rarefied.change.calcs.facet$fire, rarefied.change.calcs.facet$species.rank.med)]

fire.labs <- c("non-fire-experiencing species", "fire-experiencing")
names(fire.labs) <- c("no", "yes")

p.facet <- ggplot(rarefied.change.calcs.facet) + 
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=min.025.leg, ymax=max.975.leg), fill = "#F8766D") + # historic range in red; will show areas of range contractions
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=min.025.res, ymax=max.975.res), fill = "#00BFC4") + # modern range in blue; will show areas of range expansions
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=both.min.perc, ymax=both.max.perc), fill = "#bdbdbd") + # areas common to both in grey
  facet_grid(. ~ fire, scale="free", space="free", labeller=labeller(fire=fire.labs)) +
  scale_x_continuous(breaks=c(1:42), labels=species.labels) +
  ylim(0,2200) +
  xlab("Species") +
  ylab("Elevation (m)") +
  theme_bw() +
  theme(text=element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=90, hjust=0, vjust=0.5))
p.facet

ggsave("figures/5_elevation_ranges_2panel.pdf", p.facet, device="pdf", width=11, height=8)


### Probing what happens on burned vs unburned plots for fire-experiencing species
rarefied.change.calcs.unburned <- read_csv("data/5_range.change.calcs_unburned.csv")
rarefied.change.calcs.burned <- read_csv("data/5_range.change.calcs_burned.csv")
rarefied.change.calcs.sensitivity <- left_join(rarefied.change.calcs.unburned, rarefied.change.calcs.burned, by=c("X1", "species", "min.raw.leg", "max.raw.leg", "med.leg", "min.025.leg", "max.975.leg", "Species", "fire")) #%>% 
#select(-species) %>% 
#rename_with(~ (gsub(".x", ".unburned", .x))) %>% 
#rename_with(~ (gsub(".y", ".burned", .x)))

rarefied.change.calcs.fire.sens <- rarefied.change.calcs.sensitivity %>% 
  filter(fire=="yes") %>% 
  mutate(species.rank.med = dense_rank(med.leg),
         both.min.raw.burn = pmax(min.raw.leg, min.raw.res.y),
         both.max.raw.burn = pmin(max.raw.leg, max.raw.res.y),
         both.min.perc.burn = pmax(min.025.leg, min.025.res.y),
         both.max.perc.burn = pmin(max.975.leg, max.975.res.y),
         both.min.raw.unburn = pmax(min.raw.leg, min.raw.res.x),
         both.max.raw.unburn = pmin(max.raw.leg, max.raw.res.x),
         both.min.perc.unburn = pmax(min.025.leg, min.025.res.x),
         both.max.perc.unburn = pmin(max.975.leg, max.975.res.x))

labs = rarefied.change.calcs.fire.sens$Species
p.fire.unburned <- ggplot(rarefied.change.calcs.fire.sens) + 
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=min.025.leg, ymax=max.975.leg), fill = "#F8766D") + # historic range in red; will show areas of range contractions
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=min.025.res.x, ymax=max.975.res.x), fill = "#00BFC4") + # modern range in blue; will show areas of range expansion
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=both.min.perc.unburn, ymax=both.max.perc.unburn), fill = "#bdbdbd") + # areas common to both in grey
  ylim(0,2200) +
  scale_x_continuous(breaks=c(1:7), labels=labs) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  theme(text=element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor =   element_blank())

p.fire.burned <- ggplot(rarefied.change.calcs.fire.sens) + 
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=min.025.leg, ymax=max.975.leg), fill = "#F8766D") + # historic range in red; will show areas of range contractions
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=min.025.res.y, ymax=max.975.res.y), fill = "#00BFC4") + # modern range in blue; will show areas of range expansion
  geom_rect(aes(xmin=species.rank.med-0.33, xmax=species.rank.med+0.33, ymin=both.min.perc.burn, ymax=both.max.perc.burn), fill = "#bdbdbd") + # areas common to both in grey
  ylim(0,2200) +
  scale_x_continuous(breaks=c(1:7), labels=labs) +
  xlab("") +
  #ylab("Elevation (m)") +
  theme_bw() +
  theme(text=element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor =   element_blank())

range.fig.sens <- plot_grid(p.fire.unburned, p.fire.burned, labels=c("Unburned", "Burned"))
ggsave("figures/elevation_ranges_firespecies_sensitivity.pdf", range.fig, device="pdf", width=11, height=5)




