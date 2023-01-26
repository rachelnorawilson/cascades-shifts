# Created: Apr. 15, 2020
# Updated: June 18, 2022

# This script will be used to create 4 heatmaps:
# --> PART 1: The model-averaged coefficients
# --> PART 2: Percent +/- out of total datasets

# IMPORTANT FOR THIS SCRIPT: Heatmaps were manually saved. Export --> specify height based on preview --> choose save directory

# IMPORTANT (GENERAL) NOTE: unless otherwise indicated, always use Understory_All.csv for these analyses as it is the ONLY file with up-to-date corrections.

# Packages needed:

library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(grid)
library(gridExtra)

####### PART 1: AVERAGE AVERAGED COEFFICIENTS #########

## Step 1: Loading and tidying average coefficient data

coeff.ALLDAT <- read.csv("data/3c_new_coefficients.csv", header = TRUE)
coeff.ALLDAT[coeff.ALLDAT$Species == "EPAN", 1] <- paste("CHAN") # Correcting taxonomy issue
coeff.fire <- 
  coeff.ALLDAT[coeff.ALLDAT$Fire.Included == "Yes" & coeff.ALLDAT$Type == "Avg",
               c(1:12, 16:21)]
coeff.nofire <- coeff.ALLDAT[coeff.ALLDAT$Fire.Included == "No" & coeff.ALLDAT$Type == "Avg",
                             c(1:12, 13:15)]

## Step 2: Look at FIRE data

# Step 2a: Tidy data

coeff.fire[is.na(coeff.fire)] <- 0 # Hard code absent coeffs as 0 before averaging
coeff.avg.fire <- aggregate(coeff.fire[c("Elevation.m", 
                                         "Elevation.m2", 
                                         "Resurvey.Burned.fi", 
                                         "Resurvey.Unburned.fi", 
                                         "Elevation.m.Res.Burn.fi", 
                                         "Elevation.m.Res.Unburn.fi",
                                         "Elevation.m2.Res.Burn.fi",
                                         "Elevation.m2.Res.Unburn.fi")], 
                            by = c(coeff.fire["Species"]), mean)

mat.fire <- as.matrix(t(subset(coeff.avg.fire, select = -Species)))
colnames(mat.fire) <- coeff.avg.fire$Species
mat.fire[c("Elevation.m2.Res.Burn.fi", "Elevation.m2.Res.Unburn.fi"), 
         grepl("VAME", colnames(mat.fire))] <- NA

# Step 2b: Creating FIRE vectors for visualization

col.pallette <- rev(brewer.pal(9, "RdBu"))
col.vector.fire <- c("      ACMI", "", "      ARUV", "", 
                     "      CARU", "", "      CEVE", "", 
                     "      CHAN", "", "      PAMY", "", "      VAME", "")
row.vector.fire <- c("Elevation", 
                     expression("Elevation" ^ 2), 
                     "Burned",
                     "Unburned",
                     "Elevation * Burned",
                     "Elevation * Unburned",
                     expression("Elevation" ^ 2 * " * Burned"), 
                     expression("Elevation" ^ 2 * " * Unburned")
)
col.breaks.fire <- c(1:ncol(mat.fire))
range <- max(abs(mat.fire), na.rm = TRUE)
my.breaks <- seq(-range, range, length.out = 10)

# Step 2c: Visualization of FIRE data

pheatmap(mat.fire,
         color = col.pallette,
         cellwidth = 32,
         cellheight = 32,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = "lightgrey",
         legend = TRUE,
         #legend_breaks = leg.label.breaks.fire,
         #legend_labels = leg.label.vector.fire,
         display_numbers = TRUE,
         number_color = "black",
         na_col = "grey",
         #gaps_row = c(1:nrow(mat.fire)),
         gaps_col = col.breaks.fire,
         labels_row = row.vector.fire,
         #labels_col = col.vector.fire,
         fontsize_row = 11,
         fontsize_col = 10,
         angle_col = 0,
         breaks = my.breaks
) # Manually save


## Step 3: Look at NO FIRE data

# Step 3a: Tidy data

coeff.nofire[is.na(coeff.nofire)] <- 0 # Hard code absent coeffs as 0 before averaging
coeff.avg.nofire <- aggregate(coeff.nofire[c("Elevation.m", 
                                             "Elevation.m2", 
                                             "Data.Type.nofi", 
                                             "Data.Type.Elevation.m.nofi", 
                                             "Data.Type.Elevation.m2.nofi")], 
                              by = c(coeff.nofire["Species"]), mean)

mat.nofire <- as.matrix(t(subset(coeff.avg.nofire, select = -Species)))
colnames(mat.nofire) <- coeff.avg.nofire$Species
mat.nofire["Data.Type.Elevation.m2.nofi", 
           grepl("HODI", colnames(mat.nofire))] <- NA
mat.nofire.top <- mat.nofire[, c(1:15)] # Top half
mat.nofire.bot <- mat.nofire[, c(16:29)] # Bottom half

# Step 3b: Creating NO FIRE vectors for visualization

# Top half:

col.pallette <- rev(brewer.pal(9, "RdBu"))
col.vector.nofire.top <- c("      ACCI", "", "      ACGL", "", 
                           "      AMAL", "", "      ATFI", "", 
                           "      CAME", "", "      CHUM", "", 
                           "      CLUN", "", "      COCA", "",
                           "      GAOV", "", "      GASH", "",
                           "      GOOB", "", "      GYDR", "",
                           "      HIAL", "", "      HODI", "")
col.vector.nofire.bot <- c("      LIBO", "", "      MANE", "",
                           "      MEFE", "", "      OPHO", "",
                           "      POMU", "", "      PTAQ", "",
                           "      RULA", "", "      RUPA", "",
                           "      RUPE", "", "      RUSP", "",
                           "      SOSI", "", "      SPBE", "",
                           "      TITR", "", "      TRBO", "",
                           "      VASI", "")
row.vector.nofire <- c("Elevation", 
                       expression("Elevation" ^ 2), 
                       "Year",
                       "Elevation * Year",
                       expression("Elevation" ^ 2 * " * Year"))
col.breaks.nofire.top <- rep(1:ncol(mat.nofire.top))
col.breaks.nofire.bot <- rep(1:ncol(mat.nofire.bot))
range.nof <- max(abs(mat.nofire), na.rm = TRUE)
my.breaks <- seq(-range.nof, range.nof, length.out = 10)

# Step 3c: Visualization of NO FIRE data

no.fire.top <- pheatmap(mat.nofire.top,
                        color = col.pallette,
                        cellwidth = 30,
                        cellheight = 30,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        legend = TRUE,
                        #legend_breaks = leg.label.breaks.fire,
                        #legend_labels = leg.label.vector.fire,
                        display_numbers = TRUE,
                        number_color = "black",
                        na_col = "grey",
                        #gaps_row = c(1:nrow(mat.fire)),
                        gaps_col = col.breaks.nofire.top,
                        labels_row = row.vector.nofire,
                        #labels_col = col.vector.fire,
                        fontsize_row = 11,
                        fontsize_col = 10,
                        angle_col = 0,
                        breaks = my.breaks
)

no.fire.bot <- pheatmap(mat.nofire.bot,
                        color = col.pallette,
                        cellwidth = 30,
                        cellheight = 30,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        legend = FALSE,
                        #legend_breaks = leg.label.breaks.fire,
                        #legend_labels = leg.label.vector.fire,
                        display_numbers = TRUE,
                        number_color = "black",
                        na_col = "grey",
                        #gaps_row = c(1:nrow(mat.fire)),
                        gaps_col = col.breaks.nofire.bot,
                        labels_row = row.vector.nofire,
                        #labels_col = col.vector.fire,
                        fontsize_row = 11,
                        fontsize_col = 10,
                        angle_col = 0,
                        breaks = my.breaks
)

all.heatmaps <- list(no.fire.top[[4]], no.fire.bot[[4]])
for.plotting <- do.call(grid.arrange, all.heatmaps)
for.plotting # Manually save


######### PART 2: Percentage of + vs - coefficients out of total datasets #########

## Step 1: Loading data

perc.fire <- read.csv("data/4_pres_coefficients_percent_fire.csv", header = TRUE)
perc.fire[perc.fire$Species == "EPAN", 1] <- paste("CHAN") # Correcting taxonomy issue
perc.nofire <- read.csv("data/4_pres_coefficients_percent_NOfire.csv", header = TRUE)

## Step 2: Look at FIRE data

# Step 2a: Tidying FIRE data

mat.percfire.large <- t(subset(perc.fire, 
                               select = -c(Species, Fire.Included, Total.Sets, Effect)))
colnames(mat.percfire.large) <- paste(perc.fire$Species, perc.fire$Effect)
mat.percfire.sm <- mat.percfire.large[, !grepl(0, colnames(mat.percfire.large))]
mat.percfire.sm[, grepl("-", colnames(mat.percfire.sm))] <- 
  mat.percfire.sm[, grepl("-", colnames(mat.percfire.sm))] * -1
mat.percfire.sm[is.na(mat.percfire.sm)] <- 0
mat.percfire.sm[c("Elevation.m2.Res.Burn.fi", "Elevation.m2.Res.Unburn.fi"), 
                grepl("VAME", colnames(mat.percfire.sm))] <- NA

# If you need the order reversed
# mat.percfire.rev <- mat.percfire.sm[nrow(mat.percfire.sm):1,]

# Step 2b: Creating FIRE vectors for visualization

col.pallette <- rev(brewer.pal(9, "RdBu"))
col.vector.fire <- c("      ACMI", "", "      ARUV", "", 
                     "      CARU", "", "      CEVE", "", 
                     "      CHAN", "", "      PAMY", "", "      VAME", "")
row.vector.fire <- c("Elevation", 
                     expression("Elevation" ^ 2), 
                     "Burned",
                     "Unburned",
                     "Elevation * Burned",
                     "Elevation * Unburned",
                     expression("Elevation" ^ 2 * " * Burned"), 
                     expression("Elevation" ^ 2 * " * Unburned")
)
leg.label.breaks.fire <- c(-1, -0.5, 0, 0.5, 1)
leg.label.vector.fire <- c("100% Negative", "50% Negative", 
                           "0%", "50% Positive", "100% Positive")
col.breaks.fire <- c(2, 4, 6, 8, 10, 12, 14)

# Step 2c: Visualization of FIRE data

pheatmap(mat.percfire.sm,
         color = col.pallette,
         cellwidth = 30,
         cellheight = 30,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend = TRUE,
         legend_breaks = leg.label.breaks.fire,
         legend_labels = leg.label.vector.fire,
         # display_numbers = TRUE,
         number_color = "black",
         na_col = "grey",
         gaps_col = col.breaks.fire,
         labels_row = row.vector.fire,
         labels_col = col.vector.fire,
         fontsize_row = 12,
         fontsize_col = 15,
         angle_col = 0
) # Manually save


## Step 3: Look at NO FIRE data

# Step 3a: Tidying NO FIRE data

mat.percnofire.large <- t(subset(perc.nofire, 
                                 select = -c(Species, Fire.Included, Total.Sets, Effect)))
colnames(mat.percnofire.large) <- paste(perc.nofire$Species, perc.nofire$Effect)
mat.percnofire.sm <- mat.percnofire.large[, !grepl(0, colnames(mat.percnofire.large))]
mat.percnofire.sm[, grepl("-", colnames(mat.percnofire.sm))] <- 
  mat.percnofire.sm[, grepl("-", colnames(mat.percnofire.sm))] * -1
mat.percnofire.sm[is.na(mat.percnofire.sm)] <- 0
mat.percnofire.sm["Data.Type.Elevation.m2.nofi", 
                  grepl("HODI", colnames(mat.percnofire.sm))] <- NA
mat.percnofire.top <- mat.percnofire.sm[, c(1:30)] # Top half
mat.percnofire.bot <- mat.percnofire.sm[, c(31:58)] # Bottom half


# Step 3b: Creating NO FIRE vectors for visualization

col.pallette <- rev(brewer.pal(9, "RdBu"))
col.vector.nofire.top <- c("      ACCI", "", "      ACGL", "", 
                           "      AMAL", "", "      ATFI", "", 
                           "      CAME", "", "      CHUM", "", 
                           "      CLUN", "", "      COCA", "",
                           "      GAOV", "", "      GASH", "",
                           "      GOOB", "", "      GYDR", "",
                           "      HIAL", "", "      HODI", "",
                           "      LIBO", "")
col.vector.nofire.bot <- c("      MANE", "",
                           "      MEFE", "", "      OPHO", "",
                           "      POMU", "", "      PTAQ", "",
                           "      RULA", "", "      RUPA", "",
                           "      RUPE", "", "      RUSP", "",
                           "      SOSI", "", "      SPBE", "",
                           "      TITR", "", "      TRBO", "",
                           "      VASI", "")
row.vector.nofire <- c("Elevation", 
                       expression("Elevation" ^ 2), 
                       "Year",
                       "Elevation * Year",
                       expression("Elevation" ^ 2 * " * Year"))
leg.label.breaks.nofire <- c(-1, -0.5, 0, 0.5, 1)
leg.label.vector.nofire <- c("100% Negative", "50% Negative", 
                             "0%", "50% Positive", "100% Positive")
col.breaks.nofire.top <- seq(2, ncol(mat.percnofire.top), 2)
col.breaks.nofire.bot <- seq(2, ncol(mat.percnofire.bot), 2)

# Step 3c: Visualization of NO FIRE data

# Top half
nofire.perc.top <- pheatmap(mat.percnofire.top,
                            color = col.pallette,
                            cellwidth = 20,
                            cellheight = 20,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            legend = TRUE,
                            legend_breaks = leg.label.breaks.nofire,
                            legend_labels = leg.label.vector.nofire,
                            # display_numbers = TRUE,
                            number_color = "black",
                            na_col = "grey",
                            gaps_col = col.breaks.nofire.top,
                            labels_row = row.vector.nofire,
                            labels_col = col.vector.nofire.top,
                            fontsize_row = 12,
                            fontsize_col = 12,
                            angle_col = 0
)

# Bottom half
nofire.perc.bot <- pheatmap(mat.percnofire.bot,
                            color = col.pallette,
                            cellwidth = 20,
                            cellheight = 20,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            legend = FALSE,
                            legend_breaks = leg.label.breaks.nofire,
                            legend_labels = leg.label.vector.nofire,
                            # display_numbers = TRUE,
                            number_color = "black",
                            na_col = "grey",
                            gaps_col = col.breaks.nofire.bot,
                            labels_row = row.vector.nofire,
                            labels_col = col.vector.nofire.bot,
                            fontsize_row = 12,
                            fontsize_col = 12,
                            angle_col = 0
)

all.heatmaps.perc <- list(nofire.perc.top[[4]], nofire.perc.bot[[4]])
for.plotting.perc <- do.call(grid.arrange, all.heatmaps.perc)
for.plotting.perc # Manually save







