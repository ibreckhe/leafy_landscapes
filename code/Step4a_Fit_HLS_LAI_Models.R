## Script to fit LAI models to Harmonized Landsat-Sentinel 2 (HLS) data.
## Author: Ian Breckheimer
## Updated: 1 January 2025

## This script should be run after we are done preparing upscaled drone data and interpolated satellite observations. To see how this script fits into the larger analysis, see the documentation here: https://app.affine.pro/workspace/e094b233-4dce-4599-bc85-a12ab329bda3/aKjjrjpiHUO9Bp0A9mblh?mode=page

# Load libraries
library(dplyr)
library(mgcv)
library(ggplot2)

## Set working directory.
wd <- "/Users/ian/Library/CloudStorage/OneDrive-RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI"
#wd <- "C://Users/IanBreckheimer/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI"
setwd(wd)

## Training data for LAI model
drone_df_join <- read.csv("./output/upscaled_lai_q50_predictions_30m_HLS_interp.csv")

## Exploratory plots.
p1 <- ggplot(drone_df_join) +
  geom_point(aes(x=NDGI, y=pred_lai_q50,color=Site),size=0.2) +
  geom_smooth(aes(x=NDGI, y=pred_lai_q50), color="black", method="gam", formula=y~s(x)) +
  labs(x="HLS NDGI",y="Predicted LAI") +
  theme_bw()

pdf("./figs/HLS_LAI_model_NDGI.pdf",width=6,height=5)
p1
dev.off()

## Preliminary upscaled models
gam1 <- gamm(pred_lai_q50 ~ s(NDGI),random=list(ortho=~1), data=drone_df_join)
summary(gam1$lme)
summary(gam1$gam)
plot(gam1$gam)
saveRDS(gam1,"./output/HLS_LAI_model_gam1.rds")
