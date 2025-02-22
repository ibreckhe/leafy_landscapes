## Script to join leaf area data with plot data.
## Author: Ian Breckheimer
## Updated: 1 January 2025

## This script should be run after all field and lab data is QC'd. To see how this script fits into the larger analysis, see the documentation here: https://app.affine.pro/workspace/e094b233-4dce-4599-bc85-a12ab329bda3/lbgAKw_fYkaoE7YpIsfgv?mode=page

## Set working directory
wd <- "/Users/ian/Library/CloudStorage/OneDrive-RMBL/Documents\ -\ Research\ -\ Spatial\ Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/"
setwd(wd)

# Load necessary libraries and authenticate with Google Sheets
library(googlesheets4)
library(dplyr)
library(ggplot2)
gs4_auth(email="ibreckhe@gmail.com")

# Read leaf area data and plot data from CSV and Google Sheets, respectively
area_df_sum <- read.csv("./data/plot_leaf_areas_8_8_2024.csv")
plot_data <- read_sheet("https://docs.google.com/spreadsheets/d/1rSPt5g_-af3EsmYj_nXwJO9Ya9uL4duAkrBlrsVNMSo",
                        range="PL_data!A4:AB317")

# Format plot data for joining
plot_data$PlotID <- toupper(plot_data$PlotID)
plot_data$Num_leaves_scanned <- as.numeric(plot_data$Num_leaves_scanned)

# Join leaf area data with plot data
plot_data_area <- plot_data %>%
  select(PlotID, Species_or_FT, Collected_date, Num_leaves_scanned, Scanned_wet_mass, Scanned_dry_mass, Scanned_leaf_water_pct) %>%
  left_join(area_df_sum, by = c("PlotID" = "plot", "Species_or_FT" = "sample"))

# Calculate LMA for each plot
plot_data_area$LMA <- plot_data_area$Scanned_dry_mass / (plot_data_area$total_leaf_area_sqcm * 0.0001)

# Plot LMA for each functional type
p1 <- ggplot(filter(plot_data_area, Species_or_FT %in% c("forb", "grass", "shrub"))) +
  geom_point(aes(x = Species_or_FT, y = LMA, color = Species_or_FT),
             position = position_jitter(width = 0.2), alpha = 0.1) +
  geom_text(aes(x = Species_or_FT, y = LMA, label = PlotID, color = Species_or_FT),
            position = position_jitter(width = 0.2), size = 2) +
  scale_y_continuous("CWM Leaf Mass Per Area [g/m^2]", limits = c(0, 420)) +
  theme_bw() +
  ggtitle("Functional Group LMA")

# Plot LMA vs. leaf water content
p2 <- ggplot(filter(plot_data_area, Species_or_FT %in% c("forb", "grass", "shrub"))) +
  geom_point(aes(x = Scanned_leaf_water_pct, y = LMA, color = Species_or_FT), alpha = 0.1) +
  geom_smooth(aes(x = Scanned_leaf_water_pct, y = LMA, color = Species_or_FT), method = "lm", se = TRUE) +
  geom_text(aes(x = Scanned_leaf_water_pct, y = LMA, label = PlotID, color = Species_or_FT), size = 2) +
  scale_y_continuous("CWM Leaf Mass Per Area [g/m^2]", limits = c(0, 420)) +
  scale_x_continuous("CWM Leaf Water Content [%]") +
  theme_bw() +
  ggtitle("LMA vs. Leaf Water Content")

# Save plots as PDF
pdf("../LMA_plots.pdf", width = 8.5, height = 11)
gridExtra::grid.arrange(p1, p2, ncol = 1)
dev.off()

# Write joined data to disk
write.csv(plot_data_area, "./data/finalized_plot_data/leafy_landscapes_plot_data_v1.csv", row.names = FALSE)
