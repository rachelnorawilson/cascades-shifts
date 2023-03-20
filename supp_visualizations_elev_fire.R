# Created: May 10, 2021
# Amended: Mar 20, 2023

# This script will be used to create a violin plot of elevation vs fire history as seen in Figure S2

# IMPORTANT NOTE: unless otherwise indicated, always use Understory_All.csv for these analyses as it is the ONLY file with up-to-date corrections.

# Packages needed:

library(tidyverse)

####### PART 1: All plots regardless of species #########

load("data/1.list.fires.Rda")
burn.data <- list.fires[list.fires$Data.Type == "Resurvey",]
burn.data$Elevation.m <- as.numeric(burn.data$Elevation.m)

violin.plot.burns <- ggplot(burn.data, aes(x = Fires, y = Elevation.m)) +
  geom_violin() +
  theme_classic() +
  xlab("") +
  ylim(100, 2200) +
  ylab("Elevation (m)") +
  theme(legend.position="none", text = element_text(size = 16)) +
  stat_summary(fun=mean, geom="point", cex=2)

ggsave("figures/FigureS2_supp_elev_fire.pdf", violin.plot.burns, width=12, height=8)


