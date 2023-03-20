# Created: Feb. 16, 2021
# Modified: Feb. 1, 2022

#### This script visualizes model-averaged predictions for each rarefied dataset

library(tidyverse)
library(cowplot)

#### Read in and prepare tables of coefficients

# coefficients from all top models for each rarefied dataset that ran without warnings
# these models use the poly() formulaton for orthogonal elevation^2 terms
coeff.ALLDAT <- read.csv("data/3c_new_coefficients.csv", header = TRUE)

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
    geom_line(linewidth=3.5, linetype="dotted", show.legend=FALSE) +
    theme_classic() +
    scale_color_manual("Time x fire", values=col.pal.fire, labels=c("legacy", "resurvey, burned", "resurvey, unburned")) +
    scale_x_continuous(breaks=poly.ticks, labels=raw.ticks) +
    xlab("") + #Elevation (m)
    ylab("") #Probability of presence
  
  #ggsave(paste("figures/model.preds_ortho_",sp,".pdf",sep=""), gg, width=5, height=5)
  
  assign(paste0("preds_graph_",sp), gg)
} 

# repeat last plot with legend so that legend can be saved for multi-panel fig
gg <- ggplot(graph.dat.means, aes(x = elev.vec.lin, y = preds, color = V2)) + 
  geom_line(data=graph.dat.tall, aes(group=interaction(V2, rep), color=V2), alpha=0.08) +
  geom_line(linewidth=3.5, linetype="dotted") +
  scale_color_manual(values=col.pal.fire, labels=c("1983", "2015, burned", "2015, unburned"), guide = guide_legend(title=NULL)) +
  theme(legend.title=element_blank()) +
  theme_classic() + 
  theme(legend.key.size=unit(1.5, 'cm')) +
  theme(legend.text=element_text(size=14))

legend.fire = get_legend(gg)


#### NO-FIRE SPECIES: big loop to calculate predicted values across each rarefaction for each species

# empty matrices to store predictions
pred.leg.reps = matrix(nrow=length(elev.vec.lin),ncol=100)
pred.res.reps = matrix(nrow=length(elev.vec.lin),ncol=100)

for (i in 1:dim(species.list.nofire)[1]) {
  sp = as.list(species.list.nofire[i,1])
  mods <- coeffs.nofire %>% 
    filter(Species==sp) %>% 
    select(Int=Intercept, 
           Elev=Elevation.m, 
           Elev2=Elevation.m2, 
           Year = Data.Type.nofi,
           YearxElev = Data.Type.Elevation.m.nofi,
           YearxElev2 = Data.Type.Elevation.m2.nofi) %>% 
    mutate_if(is.numeric, ~replace(., is.na(.), 0))
  
  for (j in 1:dim(mods)[1]) {
    pred.leg.reps[,j] = mods$Int[j] + mods$Elev[j]*elev.vec.lin + mods$Elev2[j]*elev.vec.quad
    
    pred.res.reps[,j] = mods$Int[j]  + mods$Elev[j]*elev.vec.lin + mods$Elev2[j]*elev.vec.quad + mods$Year[j] + mods$YearxElev[j]*elev.vec.lin + mods$YearxElev2[j]*elev.vec.quad
  }
  
  t1.reps <- as.data.frame(cbind(elev.vec.lin, 'legacy', pred.leg.reps)) %>% 
    mutate(across(c(3:102), response))
  t1.reps.tall <- gather(t1.reps, "rep", "preds", 3:102) %>% 
    mutate(elev.vec.lin = as.numeric(elev.vec.lin))
  t1.summary <- t1.reps.tall %>% 
    group_by(V2, elev.vec.lin) %>% 
    summarise(mean.resp = mean(preds),
              lower.resp = unname(quantile(preds, c(0.05))),
              upper.resp = unname(quantile(preds, c(0.95))))
  t2.reps <- as.data.frame(cbind(elev.vec.lin, 'resurvey', pred.res.reps)) %>% 
    mutate(across(c(3:102), response))
  t2.reps.tall <- gather(t2.reps, "rep", "preds", 3:102) %>% 
    mutate(elev.vec.lin = as.numeric(elev.vec.lin))
  t2.summary <- t2.reps.tall %>% 
    group_by(V2, elev.vec.lin) %>% 
    summarise(mean.resp = mean(preds),
              lower.resp = unname(quantile(preds, c(0.05))),
              upper.resp = unname(quantile(preds, c(0.95))))
  graph.dat.tall <- bind_rows(t1.reps.tall, t2.reps.tall)  
  
  graph.dat.means <- bind_rows(t1.summary, t2.summary) %>% 
    mutate(preds=mean.resp)
  
  gg <- ggplot(graph.dat.means, aes(x = elev.vec.lin, y = preds, color = V2)) + 
    geom_line(data=graph.dat.tall, aes(group=interaction(V2, rep), color=V2), alpha=0.08, show.legend = FALSE) +
    geom_line(linewidth=3.5, linetype="dotted", show.legend=FALSE) +
    theme_classic() +
    scale_color_manual(name="TIME", values=col.pal.nofire, labels=c("legacy", "resurvey")) + 
    scale_x_continuous(breaks=poly.ticks, labels=raw.ticks) +
    xlab("") + #Elevation (m)
    ylab("") #Probability of presence
  
  #ggsave(paste("figures/model.preds_ortho_",sp,".pdf",sep=""), gg, width=5, height=5)  
  
  assign(paste0("preds_graph_",sp), gg)
}

# repeat last plot with legend so that legend can be saved for supplemental multi-panel fig
gg <- ggplot(graph.dat.means, aes(x = elev.vec.lin, y = preds, color = V2)) + 
  geom_line(data=graph.dat.tall, aes(group=interaction(V2, rep), color=V2), alpha=0.08) +
  geom_line(linewidth=3.5, linetype="dotted") +
  scale_color_manual(values=col.pal.nofire, labels=c("1983", "2015"), guide = guide_legend(title=NULL)) +
  theme(legend.title=element_blank()) +
  theme_classic() + 
  theme(legend.key.size=unit(1.5, 'cm')) +
  theme(legend.text=element_text(size=14))

legend.nofire = get_legend(gg)


#### FIGURE 4: assemble example species into multi-panel figure

# group subplots
multi <- plot_grid(legend.fire,
                   preds_graph_ARUV, #up shift fire 
                   preds_graph_VAME, #expansion fire
                   preds_graph_PAMY, #no shift fire
                   preds_graph_CHUM, #down shift no fire
                   preds_graph_OPHO, #up shift no fire
                   preds_graph_SPBE, #expansion no fire
                   preds_graph_CLUN, #no shift no fire
                   nrow=2, ncol=4,
                   labels=c("","A","B","C","D","E","F","G")) +
  theme(plot.margin = margin(50, 10, 10, 50)) #top, right, bottom, left 

multi.labs <- ggdraw(multi) + 
  draw_label("Downward shift", size=14, x=0.11, y=0.95, hjust=0) +
  draw_label("Upward shift", size=14, x=0.345, y=0.95, hjust=0) +
  draw_label("Overall expansion", size=14, x=0.57, y=0.95, hjust=0) +
  draw_label("No shift", size=14, x=0.8, y=0.95, hjust=0)

multi.x <- ggdraw(add_sub(multi.labs, "Elevation (m)", size=18, x=0.5, y=0.05, hjust=0.5, vjust=0)) 

multi.xy <- ggdraw(add_sub(multi.x, "Probability of presence", size=18, x=0.02, y=2.1, angle=90))

ggsave("figures/6_Figure4_model_preds_multipanel.pdf", multi.xy, width=12, height=8) # too large to sync w Git; placeholder uploaded instead, but can be generated locally



### all other species for supplement

multi.supp.fire <- plot_grid(legend.fire,
                             preds_graph_ACMI,
                             preds_graph_CARU,
                             NULL,
                             preds_graph_CEVE,
                             preds_graph_EPAN,
                             nrow=2, ncol=3,
                             labels=c("","ACMI","CARU","", "CEVE","CHAN"),
                             label_x=0.5,
                             label_y=1) +
  theme(plot.margin = margin(50, 10, 10, 50))

multi.supp.fire.x <- ggdraw(add_sub(multi.supp.fire, "Elevation (m)", size=18, x=0.7, y=0.4, hjust=0.5, vjust=0)) 

multi.supp.fire.xy <- ggdraw(add_sub(multi.supp.fire.x, "Probability of presence", size=18, x=0.37, y=2.1, angle=90))

ggsave("figures/6_supp_model_preds_fire.pdf", multi.supp.fire.xy, width=10, height=8) # too large to sync w Git; placeholder uploaded instead, but can be generated locally



multi.supp.nofire <- plot_grid(legend.nofire,
                               preds_graph_ACCI,
                               preds_graph_ACGL,
                               preds_graph_AMAL,
                               preds_graph_ATFI,
                               preds_graph_CAME,
                               NULL,
                               #preds_graph_CLUN,
                               preds_graph_COCA,
                               preds_graph_GAOV,
                               preds_graph_GASH,
                               preds_graph_GOOB,
                               preds_graph_GYDR,
                               NULL,
                               preds_graph_HIAL,
                               preds_graph_HODI,
                               preds_graph_LIBO,
                               preds_graph_MANE,
                               preds_graph_MEFE,
                               NULL,
                               preds_graph_POMU,
                               preds_graph_PTAQ,
                               preds_graph_RULA,
                               preds_graph_RUPA,
                               preds_graph_RUPE,
                               NULL,
                               preds_graph_RUSP,
                               preds_graph_SOSI,
                               preds_graph_TITR,
                               preds_graph_TRBO,
                               preds_graph_VASI,
                               nrow=5, ncol=6,
                               labels=c("","ACMI","ACGL","AMAL", "ATFI","CAME",
                                        "", "COCA", "GAOV", "GASH", "GOOB","GYDR", 
                                        "", "HIAL", "HODI", "LIBO", "MANE", "MEFE",
                                        "", "POMU", "PTAQ", "RULA", "RUPA", "RUPE",
                                        "", "RUSP", "SOSI", "TITR", "TRBO", "VASI"), #CLUN
                               label_x=0.5,
                               label_y=1) +
  theme(plot.margin = margin(50, 10, 10, 50))

multi.supp.nofire.x <- ggdraw(add_sub(multi.supp.nofire, "Elevation (m)", size=18, x=0.6, y=0.4, hjust=0.5, vjust=0)) 

multi.supp.nofire.xy <- ggdraw(add_sub(multi.supp.nofire.x, "Probability of presence", size=18, x=0.2, y=2.5, angle=90))

ggsave("figures/6_supp_model_preds_nofire.pdf", multi.supp.nofire.xy, width=12, height=10) # too large to sync w Git; placeholder uploaded instead, but can be generated locally


