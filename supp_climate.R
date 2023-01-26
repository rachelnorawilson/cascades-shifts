library(tidyverse)
library(patchwork)
library(lmerTest)

### read in climate files from ClimateWNA
seasonal <- read_csv("data/0_All_Plots_RNW_1970-2015ST.csv")
annual <- read_csv("data/0_All_Plots_RNW_1970-2015YT.csv")

### merge elevation and turn it to numeric
plot.names <- load("data/1.plot.names.Rda") 

seasonal <- left_join(seasonal, plot.names, by=c("ID1" = "Plot.2015"))
annual <- left_join(annual, plot.names, by=c("ID1" = "Plot.2015"))

seasonal$Elevation.m <- as.numeric(seasonal$Elevation.m)
annual$Elevation.m <- as.numeric(annual$Elevation.m)

### filter climate data to select plots that capture the breadth of elevations in the data
# 5 equal intervals for elevation
elev.vec <- seq(min(annual$Elevation.m, na.rm=T), max(annual$Elevation.m, na.rm=T), length.out=5)
# which plots match?

annual$ID1[which(abs(annual$Elevation.m-elev.vec[1])==min(abs(annual$Elevation.m-elev.vec[1]), na.rm=T))] #Bac208
annual$ID1[which(abs(annual$Elevation.m-elev.vec[2])==min(abs(annual$Elevation.m-elev.vec[2]), na.rm=T))] #Ross3020
annual$ID1[which(abs(annual$Elevation.m-elev.vec[3])==min(abs(annual$Elevation.m-elev.vec[3]), na.rm=T))] #Sour4003
annual$ID1[which(abs(annual$Elevation.m-elev.vec[4])==min(abs(annual$Elevation.m-elev.vec[4]), na.rm=T))] #Cari280
annual$ID1[which(abs(annual$Elevation.m-elev.vec[5])==min(abs(annual$Elevation.m-elev.vec[5]), na.rm=T))] #Ste2033
# except Ste is warm for its elevation; replace with another high site like Easy456 or Cas5101

plots=c("Bac208","Ross3020","Sour4003","Cari280","Easy456")
seasonal.sub <- seasonal %>% filter(ID1 %in% plots)
annual.sub <- annual %>% filter(ID1 %in% plots)

### ALTERNATIVE: filter climate data to select plots that capture the max range of differences in temperature among plots
#dat2015 <- annual %>% filter(Year==2015)
#temp.vec <- seq(min(dat2015$MAT, na.rm=T), max(dat2015$MAT, na.rm=T), length.out=5)
# which plots match?
#dat2015$ID1[which(abs(dat2015$MAT-temp.vec[1])==min(abs(dat2015$MAT-temp.vec[1]), na.rm=T))] #Cas5109, Ste2019, Ste2023, Ste2021
#dat2015$ID1[which(abs(dat2015$MAT-temp.vec[2])==min(abs(dat2015$MAT-temp.vec[2]), na.rm=T))] #What6021, What6043, What6005, Ste2031, Ste2032, RaRi5021
#dat2015$ID1[which(abs(dat2015$MAT-temp.vec[3])==min(abs(dat2015$MAT-temp.vec[3]), na.rm=T))]#"HB5158" "STE2024" "Copp6022" "Copp6028" "Copp6024" "Copp6030" "Copp6026" "Ste2078"  "Ste2062"  "Thun7120" "Thun7119" "Thun7118" "Ross3022" "Ross3024" "Ross3018" "Ste2008"  "Ruby1057" "Ruby1059" "Pan1018"  "Pan1021"  "Sour4023" "Sour4033" "Sour4021" "Cari449" "Hozo140"
#dat2015$ID1[which(abs(dat2015$MAT-temp.vec[4])==min(abs(dat2015$MAT-temp.vec[4]), na.rm=T))] #"ROSS4001" "STE2042" "Pyra1019" "Ross3008" "Ross3031" "Ross3010" "Thor217"
#dat2015$ID1[which(abs(dat2015$MAT-temp.vec[5])==min(abs(dat2015$MAT-temp.vec[5]), na.rm=T))] #Bac208

#plots=c("Cas5109","What6021","HB5158","ROSS4001","Bac208")
#seasonal.sub <- seasonal %>% filter(ID1 %in% plots)
#annual.sub <- annual %>% filter(ID1 %in% plots)

### ALTERNATIVE: filter climate data to random subset of plots
#plots <- sample(seasonal$ID1, size=10, replace=FALSE)
#seasonal.sub <- seasonal %>% 
#  filter(ID1 %in% plots) %>% 
#  select(Year, ID1, Tmax_sm, Tmax_wt, PPT_sm)
#annual.sub <- annual %>% 
#  filter(ID1 %in% plots) %>% 
#  select(Year, ID1, MAT, MAP, PAS)


### pull out survey-year data alone
seasonal.sub.survey <- seasonal.sub %>% filter(Year==1983|Year==2015)
annual.sub.survey <- annual.sub %>% filter(Year==1983|Year==2015)

MAT <- ggplot(annual.sub, aes(x=Year, y=MAT, by=ID1)) +
  geom_line(aes(x=Year, y=MAT, color=Elevation.m)) +
  geom_point(aes(x=Year, y=MAT, color=Elevation.m), data=annual.survey, cex=2) +
  scale_color_gradient(low="red", high="blue", name="Elevation (m)") +
  ylab("Mean annual temperature (\u00B0C)") +
  xlab("Year") +
  theme_classic() +
  ggtitle("A")

Twin <- ggplot(seasonal.sub, aes(x=Year, y=Tmax_wt, by=ID1)) +
  geom_line(aes(x=Year, y=Tmax_wt, color=Elevation.m)) +
  geom_point(aes(x=Year, y=Tmax_wt, color=Elevation.m), data=seasonal.survey, cex=2) +
  scale_color_gradient(low="red", high="blue", name="Elevation (m)") +
  ylab("Winter max temperature (\u00B0C)") +
  xlab("Year") +
  theme_classic() +
  ggtitle("B")

Tsum <- ggplot(seasonal.sub, aes(x=Year, y=Tmax_sm, by=ID1)) +
  geom_line(aes(x=Year, y=Tmax_sm, color=Elevation.m)) +
  geom_point(aes(x=Year, y=Tmax_sm, color=Elevation.m), data=seasonal.survey, cex=2) +
  scale_color_gradient(low="red", high="blue", name="Elevation (m)") +
  ylab("Summer max temperature (\u00B0C)") +
  xlab("Year") +
  theme_classic() +
  ggtitle("C")

MAP <- ggplot(annual.sub, aes(x=Year, y=MAP, by=ID1)) +
  geom_line(aes(x=Year, y=MAP, color=Elevation.m)) +
  geom_point(aes(x=Year, y=MAP, color=Elevation.m), data=annual.survey, cex=2) +
  scale_color_gradient(low="red", high="blue", name="Elevation (m)") +
  ylab("Mean annual precipitation (mm)") +
  xlab("Year") +
  theme_classic() +
  ggtitle("D")

Psum <- ggplot(seasonal.sub, aes(x=Year, y=PPT_sm, by=ID1)) +
  geom_line(aes(x=Year, y=PPT_sm, color=Elevation.m)) +
  geom_point(aes(x=Year, y=PPT_sm, color=Elevation.m), data=seasonal.survey, cex=2) +
  scale_color_gradient(low="red", high="blue", name="Elevation (m)") +
  ylab("Summer precipitation (mm)") +
  xlab("Year") +
  theme_classic() +
  ggtitle("E")

Pwin <- ggplot(annual.sub, aes(x=Year, y=PAS, by=ID1)) +
  geom_line(aes(x=Year, y=PAS, color=Elevation.m)) +
  geom_point(aes(x=Year, y=PAS, color=Elevation.m), data=annual.survey, cex=2) +
  scale_color_gradient(low="red", high="blue", name="Elevation (m)") +
  ylab("Precipitation as snow (mm)") +
  xlab("Year") +
  theme_classic() +
  ggtitle("F")

# combine into multi-panel
all <- (MAT | Tsum | Twin) / (MAP | Psum | Pwin) + 
  plot_layout(guides="collect") 
ggsave("figures/suppclimate.pdf", all, width=11, height=8)



### calculate changes over time

MAT.slopes <- annual.sub %>% 
  group_by(ID1) %>%
  do(model = lm(MAT ~ Year, data = .)) %>%
  mutate(MAT.slope=coef(model)["Year"],
         MAT.p=summary(model)$coefficients[2,4]) %>% 
  select(-model)
#0.021C per year --> ~0.7C for 32 years (1983-2015)

MAT.slopes.all <- lmer(MAT ~ Year + (1|ID1), data = annual) 
summary(MAT.slopes.all)
# +0.022 C per year, P<0.0001

MAP.slopes <- annual.sub %>% 
  group_by(ID1) %>%
  do(model = lm(MAP ~ Year, data = .)) %>%
  mutate(MAP.slope=coef(model)["Year"], 
         MAP.p=summary(model)$coefficients[2,4]) %>% 
  select(-model) 

MAP.slopes.all <- lmer(MAP ~ Year + (1|ID1), data = annual) 
summary(MAP.slopes.all)
# -0.2013 mm per year, P=0.249

Tsum.slopes <- seasonal.sub %>% 
  group_by(ID1) %>%
  do(model = lm(Tmax_sm ~ Year, data = .)) %>%
  mutate(Tsum.slope=coef(model)["Year"],
         Tsum.p=summary(model)$coefficients[2,4]) %>% 
  select(-model)

Tsum.slopes.all <- lmer(Tmax_sm ~ Year + (1|ID1), data = seasonal) 
summary(Tsum.slopes.all)
# -0.01687C per year, P<0.0001

Twin.slopes <- seasonal.sub %>% 
  group_by(ID1) %>%
  do(model = lm(Tmax_wt ~ Year, data = .)) %>%
  mutate(Twin.slope=coef(model)["Year"],
         Twin.p=summary(model)$coefficients[2,4]) %>% 
  select(-model)
#0.025C per year --> ~0.8C for 32 years (1983-2015)

Twin.slopes.all <- lmer(Tmax_wt ~ Year + (1|ID1), data = seasonal) 
summary(Twin.slopes.all)
# -0.02562C per year, P<0.0001

Psum.slopes <- seasonal.sub %>% 
  group_by(ID1) %>%
  do(model = lm(PPT_sm ~ Year, data = .)) %>%
  mutate(Psum.slope=coef(model)["Year"],
         Psum.p=summary(model)$coefficients[2,4]) %>% 
  select(-model)

Psum.slopes.all <- lmer(PPT_sm ~ Year + (1|ID1), data = seasonal) 
summary(Psum.slopes.all)
# -0.625mm per year, P<0.0001

Pwin.slopes <- annual.sub %>% 
  group_by(ID1) %>%
  do(model = lm(PAS ~ Year, data = .)) %>%
  mutate(PAS.slope=coef(model)["Year"],
         PAS.p=summary(model)$coefficients[2,4]) %>% 
  select(-model)
#-5 per year --> -160 for 32 years (1983-2015)

Pwin.slopes.all <- lmer(PAS ~ Year + (1|ID1), data = annual) 
summary(Pwin.slopes.all)
# -4.9549mm per year, P<0.0001


Ann.diffs <- annual.sub %>% 
  group_by(ID1) %>%
  #filter(Year==1983 | Year==2015) %>% 
  #pivot_wider(names_from=Year, values_from=c(MAT, MAP, PAS)) %>% 
  mutate(MAT.diff = MAT[Year==2015]-MAT[Year==1983],
         MAP.diff.raw = MAP[Year==2015]-MAP[Year==1983],
         MAP.diff.perc = MAP.diff.raw/MAP[Year==1983]*100,
         PAS.diff.raw = PAS[Year==2015]-PAS[Year==1983],
         PAS.diff.perc = PAS.diff.raw/PAS[Year==1983]*100) %>% 
  select(ID1, MAT.diff, MAP.diff.raw, MAP.diff.perc, PAS.diff.raw, PAS.diff.perc) %>% 
  unique()

Seas.diffs <- seasonal.sub %>% 
  group_by(ID1) %>%
  #filter(Year==1983 | Year==2015) %>% 
  #pivot_wider(names_from=Year, values_from=c(Tmax_sm, Tmax_wt, PPT_sm)) %>% 
  mutate(Tsum.diff = Tmax_sm[Year==2015]-Tmax_sm[Year==1983],
         Twin.diff = Tmax_wt[Year==2015]-Tmax_wt[Year==1983],
         Psum.diff = PPT_sm[Year==2015]-PPT_sm[Year==1983],
         Psum.diff.perc = Psum.diff/PPT_sm[Year==1983]*100) %>% 
  select(ID1, Tsum.diff, Twin.diff, Psum.diff, Psum.diff.perc) %>% 
  unique()


annual.survey <- annual %>% filter(Year==1983 | Year==2015)
seasonal.survey <- seasonal %>% filter(Year==1983 | Year==2015)

MAT.diffs.all <- lmer(MAT ~ as.factor(Year) + (1|ID1), data = annual.survey) 
summary(MAT.diffs.all)
# +2.032 C, P<0.0001

Tsum.diffs.all <- lmer(Tmax_sm ~ as.factor(Year) + (1|ID1), data = seasonal.survey) 
summary(Tsum.diffs.all)
# +4.87373 C, P<0.0001

Twin.diffs.all <- lmer(Tmax_wt ~ as.factor(Year) + (1|ID1), data = seasonal.survey) 
summary(Twin.diffs.all)
# +1.99115 C, P<0.0001

MAP.diffs.all <- lmer(MAP ~ as.factor(Year) + (1|ID1), data = annual.survey) 
summary(MAP.diffs.all)
# -136.721 mm, P<0.0001
MAP.diffs.all.perc <- lmer(MAP/MAP[Year==1983]*100 ~ as.factor(Year) + (1|ID1), data = annual.survey) 
summary(MAP.diffs.all.perc)
# -6.7%, P<0.0001

Psum.diffs.all <- lmer(PPT_sm ~ as.factor(Year) + (1|ID1), data = seasonal.survey) 
summary(Psum.diffs.all)
# -184.048 mm, P<0.0001
Psum.diffs.all.perc <- lmer(PPT_sm/PPT_sm[Year==1983]*100 ~ as.factor(Year) + (1|ID1), data = seasonal.survey) 
summary(Psum.diffs.all.perc)
# -55%, P<0.0001

PAS.diffs.all <- lmer(PAS ~ as.factor(Year) + (1|ID1), data = annual.survey) 
summary(PAS.diffs.all)
# -371.070 mm, P<0.0001
PAS.diffs.all.perc <- lmer(PAS/PAS[Year==1983]*100 ~ as.factor(Year) + (1|ID1), data = annual.survey) 
summary(PAS.diffs.all.perc)
# -48%, P<0.0001

Diffs <- left_join(Ann.diffs, Seas.diffs)

Slopes <- left_join(left_join(left_join(left_join(left_join(MAT.slopes, Tsum.slopes), Twin.slopes), MAP.slopes), Psum.slopes), Pwin.slopes)

write_csv(Diffs, "data/supp_climate_diffs.csv")
write_csv(Slopes, "data/supp_climate_slopes.csv")
