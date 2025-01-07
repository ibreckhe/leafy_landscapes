## Script to fit drone-based LAI models and predict LAI using drone imagery.
## Author: Ian Breckheimer
## Updated: 1 January 2025

## This script should be run after drone data has been extracted for field plots. To see how this script fits into the larger analysis, see the documentation here: https://app.affine.pro/workspace/e094b233-4dce-4599-bc85-a12ab329bda3/hakWMNfZ-ZXSIhQycXIYx?mode=page

## Set up workspace.
library(rstanarm)
library(sf)
library(readr)
library(dplyr)
library(terra)
library(ggplot2)
setwd("~/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/")

# Read in the plot data
plot_lab_data_gcc <- read_csv("./data/finalized_plot_data/leafy_landscapes_plot_drone_data_v1.csv")

# Convert collected date to day of year and create a binary variable for early collection
plot_lab_data_gcc$Collected_DOY <- as.numeric(format(as.Date(plot_lab_data_gcc$Collected_date), "%j"))
plot_lab_data_gcc$Collected_early <- plot_lab_data_gcc$Collected_DOY < 185

# Create a variable for the site visit number
plot_lab_data_gcc$Site_visit_num <- as.factor(gsub("[^0-9]", "", plot_lab_data_gcc$visit_code))

# Makes the plots.
ggplot(plot_lab_data_gcc)+
  geom_point(aes(x=Plot_adjusted_LAI,y=Plot_mean_LMA,color=Site_visit_num)) + 
  geom_smooth(aes(x=Plot_adjusted_LAI,y=Plot_mean_LMA,color=Site_visit_num),
              method="lm",se=FALSE) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

ggplot(plot_lab_data_gcc)+
  geom_point(aes(x=plot_gcc_mean,y=Plot_total_canopy_water_content,color=Site_visit_num)) + 
  geom_smooth(aes(x=plot_gcc_mean,y=Plot_total_canopy_water_content,color=Site_visit_num),
              method="gam",se=FALSE) + 
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw()

ggplot(plot_lab_data_gcc)+
  geom_point(aes(x=Plot_adjusted_LAI,y=Plot_total_dry_mass, color=Site_visit_num)) + 
  geom_smooth(aes(x=Plot_adjusted_LAI,y=Plot_total_dry_mass, color=Site_visit_num),
              method="lm",se=FALSE) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()


ggplot(plot_lab_data_gcc)+
  geom_point(aes(x=plot_gcc_mean,y=Plot_mean_LMA,color=Site_visit_num)) + 
  geom_smooth(aes(x=plot_gcc_mean,y=Plot_mean_LMA,color=Site_visit_num),
              method="lm",se=FALSE) + 
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw()

ggplot(plot_lab_data_gcc)+
  geom_point(aes(x=plot_gcc_mean,y=Plot_adjusted_LAI,color=Collected_early)) + 
  geom_smooth(aes(x=plot_gcc_mean,y=Plot_adjusted_LAI,color=Collected_early),
              method="lm",se=FALSE) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

ggplot(plot_lab_data_gcc)+
  geom_point(aes(x=Plot_canopy_water_pct_logit,y=Plot_adjusted_LAI,color=Site_visit_num)) + 
  geom_smooth(aes(x=Plot_canopy_water_pct_logit,y=Plot_adjusted_LAI,color=Site_visit_num),
              method="lm",se=FALSE) + 
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw()

ggplot(plot_lab_data_gcc)+
  geom_point(aes(x=plot_gcc_mean,y=Plot_canopy_water_pct_logit,color=Site_visit_num)) + 
  geom_smooth(aes(x=plot_gcc_mean,y=Plot_canopy_water_pct_logit,color=Site_visit_num),
              method="lm",se=FALSE) + 
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw()

ggplot(plot_lab_data_gcc)+
  geom_point(aes(x=Plot_canopy_water_pct,y=Plot_mean_LMA,color=visit_code)) + 
  geom_text(aes(x=Plot_leaf_water_pct,y=Plot_mean_LMA,label=PlotID),
            method="lm",se=FALSE) + 
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw()

ggplot(plot_lab_data_gcc)+
  geom_point(aes(x=Plot_leaf_water_pct,y=Plot_mean_LMA,color=elev)) + 
  #geom_text(aes(x=Plot_leaf_water_pct,y=Plot_mean_LMA,label=PlotID),
  #          method="lm",se=FALSE) + 
  #scale_x_log10() +
  #scale_y_log10() +
  theme_bw()


## Preliminary models predicting LAI.
plot(plot_lab_data_gcc$crop_gcc_mean,plot_lab_data_gcc$Plot_adjusted_LAI)
plot(plot_lab_data_gcc$plot_gcc_mean,plot_lab_data_gcc$Plot_adjusted_LAI)
plot(plot_lab_data_gcc$crop_ht_sd,plot_lab_data_gcc$Plot_adjusted_LAI)
plot(plot_lab_data_gcc$crop_gcc_ht_sd,plot_lab_data_gcc$Plot_adjusted_LAI)
plot(plot_lab_data_gcc$crop_green_sd,plot_lab_data_gcc$Plot_adjusted_LAI)
plot(plot_lab_data_gcc$plot_bcc_mean,plot_lab_data_gcc$Plot_adjusted_LAI)
plot(plot_lab_data_gcc$plot_bcc_mean,plot_lab_data_gcc$plot_gcc_mean)
plot(plot_lab_data_gcc$crop_green_focal_sd - plot_lab_data_gcc$plot_green_sd ,plot_lab_data_gcc$Plot_adjusted_LAI)

lai_mod1 <- lm(Plot_adjusted_LAI ~ plot_gcc_mean + plot_green_sd + plot_ht_sd + cloudy_yn + cloudy_yn:plot_gcc_mean +
                 Site + Site:plot_gcc_mean,
               data=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E"))))


summary(lai_mod1)
plot_lab_data_gcc$pred_lai_mod1 <- NA
plot_lab_data_gcc$pred_lai_mod1[!is.na(plot_lab_data_gcc$Plot_adjusted_LAI) & !(plot_lab_data_gcc$PlotID %in% c("SN3A","SN1E"))] <- predict(lai_mod1)
plot_lab_data_gcc$pred_lai_resid_mod1 <- plot_lab_data_gcc$Plot_adjusted_LAI - plot_lab_data_gcc$pred_lai_mod1

plot(plot_lab_data_gcc$pred_lai_mod1,plot_lab_data_gcc$Plot_adjusted_LAI,cex=0)
text(plot_lab_data_gcc$pred_lai_mod1,plot_lab_data_gcc$Plot_adjusted_LAI,
     labels=plot_lab_data_gcc$PlotID,cex=0.8,col=as.numeric(factor(plot_lab_data_gcc$visit_code)))
abline(0,1,lty="dotted")


library(loo)
library(shinystan)

lai_mod1 <- stan_glmer((Plot_adjusted_LAI + 0.4) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + crop_ht_sd + cloudy_yn + cloudy_yn:plot_gcc_mean +
                         (1+plot_gcc_mean_sq | Site), family=Gamma(link="log"),
                       data=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN4G"))))
summary(lai_mod1)
#launch_shinystan(lai_mod1)

lai_mod2 <- stan_glmer((Plot_adjusted_LAI + 0.4) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + crop_ht_sd + cloudy_yn + cloudy_yn:plot_gcc_mean +
                         (1+plot_gcc_mean_sq | Site) + (-1 + plot_gcc_mean | Site) + (1 | visit_code), family=Gamma(link="log"),
                       data=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN4G"))))
summary(lai_mod2)
#launch_shinystan(lai_mod2)

lai_mod3 <- stan_glmer((Plot_adjusted_LAI + 0.4) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + plot_ht_sd +
                         (1+plot_gcc_mean_sq | Site) + (1 | visit_code), family=Gamma(link="log"),
                       data=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN4G"))))
summary(lai_mod3)
#launch_shinystan(lai_mod3)

lai_mod4 <- stan_glmer((Plot_adjusted_LAI + 0.4) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + plot_ht_sd +
                         (1+plot_gcc_mean_sq | Site) + (-1 + plot_gcc_mean | Site) + (1 | visit_code), family=Gamma(link="log"),
                       data=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN4G"))))
summary(lai_mod4)
#launch_shinystan(lai_mod4)

# lai_mod5 <- stan_gamm4((Plot_adjusted_LAI + 0.4) ~ s(plot_gcc_mean) + plot_green_sd + plot_ht_sd,
#                           random=~(1 + plot_gcc_mean | Site) + (1 | visit_code), family=Gamma(link="log"),
#                        data=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN4G"))))
# summary(lai_mod5)
# launch_shinystan(lai_mod5)
loo_compare(loo(lai_mod1),loo(lai_mod2),loo(lai_mod3),loo(lai_mod4))


trait_mvmod1 <- stan_mvmer(list((Plot_adjusted_LAI+0.6) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + plot_ht_sd + (1+plot_gcc_mean | Site.x) + (1 | visit_code),
                                (Plot_total_canopy_water_content+50) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + plot_ht_sd + (1+plot_gcc_mean | Site.x) + (1 | visit_code),
                                (Plot_total_dry_mass+50) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + plot_ht_sd + (1+plot_gcc_mean | Site.x) + (1 | visit_code)), 
                           family=list(Gamma(link="log"), Gamma(link="log"), Gamma(link="log")),
                           data=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN4G"))),
                           iter=5000,thin=5)
summary(trait_mvmod1)

trait_mvmod2 <- stan_mvmer(list((Plot_adjusted_LAI+0.6) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + plot_ht_sd + (1+plot_gcc_mean_sq | Site.x) + (1 | visit_code),
                                (Plot_total_canopy_water_content+10) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + plot_ht_sd + (1+plot_gcc_mean_sq | Site.x) + (1 | visit_code),
                                (Plot_total_dry_mass+10) ~ plot_gcc_mean + plot_gcc_mean_sq + plot_green_sd + plot_ht_sd + (1+plot_gcc_mean_sq | Site.x) + (1 | visit_code)), 
                           family=list(Gamma(link="log"), Gamma(link="log"), Gamma(link="log")),
                           data=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN4G"))))
summary(trait_mvmod2)

loo_compare(loo(trait_mvmod1),loo(trait_mvmod2))

## launch_shinystan(trait_mvmod1) # this is not working

lai_pred1 <- posterior_predict(lai_mod1,newdata=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN1A","SN4G"))),
                               re.form=~0)
lai_pred1[is.na(lai_pred1)] <- 0
lai_pred1[is.nan(lai_pred1)] <- 0
lai_quant1 <- t(apply(lai_pred1,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))) - 0.4
summary(lai_quant1[,6] - lai_quant1[,2])

lai_pred2 <- posterior_predict(lai_mod2,newdata=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN1A","SN4G"))),
                               re.form=~(1+plot_gcc_mean_sq | Site) + (1 | visit_code))
lai_pred2[is.na(lai_pred2)] <- 0
lai_pred2[is.nan(lai_pred2)] <- 0
lai_quant2 <- t(apply(lai_pred2,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))) - 0.4
summary(lai_quant2[,6] - lai_quant2[,2])

lai_pred3 <- posterior_predict(lai_mod3,newdata=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN1A","SN4G"))),
                               re.form=~(1+plot_gcc_mean_sq | Site) + (1 | visit_code))
lai_pred3[is.na(lai_pred3)] <- 0
lai_pred3[is.nan(lai_pred3)] <- 0
lai_quant3 <- t(apply(lai_pred3,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))) - 0.4
summary(lai_quant3[,6] - lai_quant3[,2])

lai_pred4 <- posterior_predict(lai_mod4,newdata=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN1A","SN4G"))),
                               re.form=~(1+plot_gcc_mean_sq | Site) + (-1 + plot_gcc_mean | Site) + (1 | visit_code))
lai_pred4[is.na(lai_pred4)] <- 0
lai_pred4[is.nan(lai_pred4)] <- 0
lai_quant4 <- t(apply(lai_pred4,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))) - 0.4
summary(lai_quant4[,6] - lai_quant4[,2])

lai_pred5<- posterior_predict(trait_mvmod1, m=1, newdata=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN1A","SN4G"))))
lai_pred5[is.na(lai_pred5)] <- 0
lai_pred5[is.nan(lai_pred5)] <- 0
lai_quant5 <- t(apply(lai_pred5,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))) - 0.6
summary(lai_quant5[,6] - lai_quant5[,2])

cwc_pred5<- posterior_predict(trait_mvmod1, m=2, newdata=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN1A","SN4G"))))
cwc_pred5[is.na(cwc_pred5)] <- 0
cwc_pred5[is.nan(cwc_pred5)] <- 0
cwc_quant5 <- t(apply(cwc_pred5,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))) - 50
summary(cwc_quant5[,6] - cwc_quant5[,2])

mass_pred5<- posterior_predict(trait_mvmod1, m=3, newdata=filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN1A","SN4G"))))
mass_pred5[is.na(mass_pred5)] <- 0
mass_pred5[is.nan(mass_pred5)] <- 0
mass_quant5 <- t(apply(mass_pred5,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))) - 50
summary(mass_quant5[,6] - mass_quant5[,2])

## Filtered data.
plot_data_filt <- filter(plot_lab_data_gcc,!(PlotID %in% c("SN3A","SN1E","LT2B","SH3C","SB2F","SN1A","SN4G")))
plot_data_filt$pred_lai_mod3_q10 <- lai_quant3[,2]
plot_data_filt$pred_lai_mod3_q25 <- lai_quant3[,3]
plot_data_filt$pred_lai_mod3_q50 <- lai_quant3[,4]
plot_data_filt$pred_lai_mod3_q75 <- lai_quant3[,5]
plot_data_filt$pred_lai_mod3_q90 <- lai_quant3[,6]

library(ggplot2)
ggplot(plot_data_filt)+
  geom_text(aes(x=Plot_adjusted_LAI,y=pred_lai_mod3_q50,label=PlotID),alpha=1,cex=2.9) +
  geom_point(aes(x=Plot_adjusted_LAI,y=pred_lai_mod3_q50),alpha=0.1) +
  #geom_linerange(aes(x=Plot_adjusted_LAI,ymin=pred_lai_mod3_q10,ymax=pred_lai_mod3_q90),lwd=0.2,alpha=0.05) + 
  geom_linerange(aes(x=Plot_adjusted_LAI,ymin=pred_lai_mod3_q25,ymax=pred_lai_mod3_q75),lwd=0.8,alpha=0.1) + 
  geom_abline(intercept=0,slope=1,lty="dotted") + 
  #facet_wrap(~visit_code)+
  theme_bw()

plot_data_filt$pred_lai_mod4_q10 <- lai_quant4[,2]
plot_data_filt$pred_lai_mod4_q25 <- lai_quant4[,3]
plot_data_filt$pred_lai_mod4_q50 <- lai_quant4[,4]
plot_data_filt$pred_lai_mod4_q75 <- lai_quant4[,5]
plot_data_filt$pred_lai_mod4_q90 <- lai_quant4[,6]

library(ggplot2)
ggplot(plot_data_filt)+
  geom_text(aes(x=Plot_adjusted_LAI,y=pred_lai_mod4_q50,label=PlotID),alpha=1,cex=1.9) +
  geom_point(aes(x=Plot_adjusted_LAI,y=pred_lai_mod4_q50),alpha=0.1) +
  #geom_linerange(aes(x=Plot_adjusted_LAI,ymin=pred_lai_mod4_q10,ymax=pred_lai_mod4_q90),lwd=0.2,alpha=0.05) + 
  geom_linerange(aes(x=Plot_adjusted_LAI,ymin=pred_lai_mod4_q25,ymax=pred_lai_mod4_q75),lwd=0.8,alpha=0.1) + 
  geom_abline(intercept=0,slope=1,lty="dotted") + 
  #facet_wrap(~visit_code)+
  theme_bw()

plot_data_filt$pred_lai_mod5_q10 <- lai_quant5[,2]
plot_data_filt$pred_lai_mod5_q25 <- lai_quant5[,3]
plot_data_filt$pred_lai_mod5_q50 <- lai_quant5[,4]
plot_data_filt$pred_lai_mod5_q75 <- lai_quant5[,5]
plot_data_filt$pred_lai_mod5_q90 <- lai_quant5[,6]

p1 <- ggplot(plot_data_filt)+
  geom_text(aes(x=Plot_adjusted_LAI,y=pred_lai_mod5_q50,label=PlotID),alpha=1,cex=1.9) +
  geom_point(aes(x=Plot_adjusted_LAI,y=pred_lai_mod5_q50),alpha=0.1) +
  #geom_linerange(aes(x=Plot_adjusted_LAI,ymin=pred_lai_mod3_q10,ymax=pred_lai_mod3_q90),lwd=0.2,alpha=0.05) + 
  geom_linerange(aes(x=Plot_adjusted_LAI,ymin=pred_lai_mod5_q25,ymax=pred_lai_mod5_q75),lwd=0.8,alpha=0.1) + 
  geom_abline(intercept=0,slope=1,lty="dotted") + 
  scale_x_continuous("Measured Leaf Area Index (Unitless)")+
  scale_y_continuous("UAS Predicted Leaf Area Index (Unitless)")+
  #facet_wrap(~visit_code)+
  theme_bw()

plot_data_filt$pred_cwc_mod5_q10 <- cwc_quant5[,2]
plot_data_filt$pred_cwc_mod5_q25 <- cwc_quant5[,3]
plot_data_filt$pred_cwc_mod5_q50 <- cwc_quant5[,4]
plot_data_filt$pred_cwc_mod5_q75 <- cwc_quant5[,5]
plot_data_filt$pred_cwc_mod5_q90 <- cwc_quant5[,6]

p2 <- ggplot(plot_data_filt)+
  geom_text(aes(x=Plot_total_canopy_water_content,
                y=pred_cwc_mod5_q50, label=PlotID),alpha=1,cex=1.9) +
  geom_point(aes(x=Plot_total_canopy_water_content,
                 y=pred_cwc_mod5_q50), alpha=0.1) +
  #geom_linerange(aes(x=Plot_adjusted_LAI,ymin=pred_lai_mod3_q10,ymax=pred_lai_mod3_q90),lwd=0.2,alpha=0.05) + 
  geom_linerange(aes(x=Plot_total_canopy_water_content,
                     ymin=pred_cwc_mod5_q25,
                     ymax=pred_cwc_mod5_q75),lwd=0.8,alpha=0.1) + 
  geom_abline(intercept=0,slope=1,lty="dotted") + 
  scale_x_continuous("Measured Canopy Water Content (g/m^2)")+
  scale_y_continuous("UAS Predicted Canopy Water Content (g/m^2)")+
  #facet_wrap(~visit_code)+
  theme_bw()

plot_data_filt$pred_mass_mod5_q10 <- mass_quant5[,2]
plot_data_filt$pred_mass_mod5_q25 <- mass_quant5[,3]
plot_data_filt$pred_mass_mod5_q50 <- mass_quant5[,4]
plot_data_filt$pred_mass_mod5_q75 <- mass_quant5[,5]
plot_data_filt$pred_mass_mod5_q90 <- mass_quant5[,6]


p3 <- ggplot(plot_data_filt)+
  geom_text(aes(x=Plot_total_dry_mass, y=pred_mass_mod5_q50,label=PlotID),alpha=1,cex=1.9) +
  geom_point(aes(x=Plot_total_dry_mass,y=pred_mass_mod5_q50),alpha=0.1) +
  #geom_linerange(aes(x=Plot_adjusted_LAI,ymin=pred_lai_mod3_q10,ymax=pred_lai_mod3_q90),lwd=0.2,alpha=0.05) + 
  geom_linerange(aes(x=Plot_total_dry_mass,ymin=pred_mass_mod5_q25,ymax=pred_mass_mod5_q75),lwd=0.8,alpha=0.1) + 
  geom_abline(intercept=0,slope=1,lty="dotted") +
  scale_x_continuous("Measured AGB (g/m^2)")+
  scale_y_continuous("UAS Predicted AGB (g/m^2)")+
  #facet_wrap(~visit_code)+
  theme_bw()

pdf("./figs/lai_cwc_biomass_mv_model5_obs_vs_pred.pdf",width=4,height=10)
gridExtra::grid.arrange(p1,p2,p3)
dev.off()

##3D plot of observed data and model predictions.
library(rgl)
plot3d(plot_data_filt$Plot_adjusted_LAI,plot_data_filt$Plot_total_canopy_water_content,plot_data_filt$Plot_total_dry_mass)
plot3d(plot_data_filt$pred_lai_mod5_q50,plot_data_filt$pred_cwc_mod5_q50,plot_data_filt$pred_mass_mod5_q50)

## Computes residuals for best Bayesian model.
plot_data_filt$pred_lai_mod3_resid <- plot_data_filt$Plot_adjusted_LAI - plot_data_filt$pred_lai_mod3_q50
summary(lm(plot_data_filt$Plot_adjusted_LAI ~ plot_data_filt$pred_lai_mod3_q50))


#### Computes aggregated predictors each mosaic.
# bound <- vect("./data/HLS/HLS_boundary_UTM.gpkg")
# hls_scene <- rast("./data/HLS/HLS.S30.T13SCD.2024186T174911.v2.0.B8A.tif")
# hls_cropped <- crop(hls_scene,ext(bound))
# hls_cropped[] <- 1
# names(hls_cropped) <- "HLS_mask"
# writeRaster(hls_cropped,"./data/HLS/HLS_mask_UER.tif",overwrite=TRUE)
hls_mask_30m <- rast("./data/HLS/HLS_mask_UER.tif")
#hls_mask_1m <- disagg(hls_mask_30m,fact=30,method="near",filename="./data/HLS/HLS_mask_UER_1m.tif",overwrite=TRUE)
hls_mask_1m <- rast("./data/HLS/HLS_mask_UER_1m.tif")

for(i in 1:nrow(image_paths)){
  
  ## Checks to see if output exists.
  outname <- gsub("/ortho/", "/dronedata_predictors/", image_paths[i,]$ortho_path)
  outname <- gsub("_ortho.tif", "_predictors_1m.tif", outname)
  if(file.exists(outname)){
    print(paste(outname, "exists, skipping..."))
    next
  }else{
    print(paste("Preparing", outname))
  }
  
  ## Reads in rasters
  print(paste("Reading inputs..."))
  pred_ortho <- rast(image_paths[i,]$ortho_path)
  pred_dsm <- rast(image_paths[i,]$dsm_path)
  pred_dem <- rast(image_paths[i,]$dem_path)
  pred_dem_1m <- project(pred_dem,hls_mask_1m,method="average",align_only=TRUE)
  
  ## Reprojects to UTM coordinates.
  print(paste("Reprojecting inputs..."))
  pred_dem_proj <- project(pred_dem,pred_dsm,method="bilinear",filename="./scratch/pred_dem_proj.tif",overwrite=TRUE)
  dsm_fun <- function(x){x[1] + 15.55 - x[2]}
  pred_dsm_ht <- app(c(pred_dsm,pred_dem_proj),fun=dsm_fun, filename="./scratch/pred_dsm_ht.tif",
                     overwrite=TRUE,cores=6)
  pred_dsm_ht_sd_proj <- aggregate(pred_dsm_ht,fact=33,fun="sd")
  pred_dsm_ht_sd_1m <- project(pred_dsm_ht_sd_proj,pred_dem_1m,method="bilinear")
  
  ## Computes GCC and squared gcc.
  print(paste("Computing predictors..."))
  pred_ortho_proj <- project(pred_ortho,pred_dem,method="bilinear",filename="./scratch/pred_ortho_proj.tif",
                             overwrite=TRUE)
  gcc_fun <- function(x){x[2]/(x[1]+x[2]+x[3])}
  pred_gcc_proj <- app(pred_ortho_proj,fun=gcc_fun,filename="./scratch/pred_gcc_proj.tif",
                       overwrite=TRUE, cores=6)
  pred_gcc_mean <- aggregate(pred_gcc_proj,fact=33,fun="mean")
  pred_gcc_mean_1m <- project(pred_gcc_mean,pred_dem_1m,method="bilinear")
  pred_gcc_mean_sq_1m <- pred_gcc_mean_1m ^ 2
  
  ## Computes greenness sd.
  pred_green_sd <- aggregate(pred_ortho_proj[[2]],fact=33,fun="sd")
  pred_green_sd_1m <- project(pred_green_sd,pred_dem_1m,method="bilinear")
  
  ## Creates a mask for height sd far outside of the range of the data.
  ht_mask <- pred_dsm_ht_sd_1m < 1.0
  
  ## Writes out predictors.
  print(paste("Writing output..."))
  preds_1m <- mask(c(pred_dsm_ht_sd_1m, pred_gcc_mean_1m, 
                     pred_gcc_mean_sq_1m, pred_green_sd_1m),ht_mask,maskvalues=0)
  names(preds_1m) <- c("plot_ht_sd","plot_gcc_mean", "plot_gcc_mean_sq", "plot_green_sd")
  writeRaster(preds_1m, file=outname, overwrite=TRUE)
  rm(pred_ortho_proj,pred_gcc_proj,pred_dem_proj,pred_dsm_ht)
}


## Predicts LAI using the best Bayesian model.
site_preds <- list.files("./data/dronedata_predictors/",pattern="1m.tif$",full.names=TRUE)
#for(i in 1:4){
for(i in 1:length(site_preds)){
  
  outname <- gsub("/data/ortho/", "/output/lai_predictions/", image_paths[i,]$ortho_path)
  outname <- gsub("_ortho.tif", "_predicted_LAI_1m_mod5.tif", outname)
  
  #if(file.exists(outname)){
  #  print(paste(outname, "exists, skipping..."))
  #  next
  #}else{
  #  print(paste("Creating model predictions for", outname))
  #}
  
  test_preds_1m <- rast(site_preds[i])
  
  test_preds_1m_df <- as.data.frame(test_preds_1m,xy=TRUE)
  test_preds_1m_df$cloudy_yn <- image_paths[i,]$cloudy_yn
  test_preds_1m_df$Site.x <- image_paths[i,]$Site
  test_preds_1m_df$visit_code <- image_paths[i,]$visit_code
  test_preds_1m_df_comp <- test_preds_1m_df[complete.cases(test_preds_1m_df),]
  
  test_lai_pred <- posterior_predict(trait_mvmod1,m=1, newdata=test_preds_1m_df_comp, draws=100,
                                     re.form=~(1+plot_gcc_mean | Site.x)) - 0.6
  test_lai_quants <- t(apply(test_lai_pred,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975),na.rm=TRUE))
  
  ## Sets predictions for very low gcc to zero.
  test_lai_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,4] <- 0 
  test_lai_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,3] <- 0 
  test_lai_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,2] <- 0 
  test_lai_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,1] <- 0 
  
  out_df_q50 <- data.frame(x=test_preds_1m_df_comp$x, 
                           y=test_preds_1m_df_comp$y,
                           pred_lai_q50=test_lai_quants[,4])
  
  out_rast_q50 <- rast(out_df_q50)
  
  out_df_q10 <- data.frame(x=test_preds_1m_df_comp$x, 
                           y=test_preds_1m_df_comp$y,
                           pred_lai_q10=test_lai_quants[,2])
  out_rast_q10 <- rast(out_df_q10)
  
  out_df_q90 <- data.frame(x=test_preds_1m_df_comp$x, 
                           y=test_preds_1m_df_comp$y,
                           pred_lai_q90=test_lai_quants[,6])
  out_rast_q90 <- rast(out_df_q90)
  out_rast_int80 <- out_rast_q90 - out_rast_q10
  names(out_rast_int80) <- "pred_lai_int80"
  
  out_rast <- c(out_rast_q50,out_rast_q10,out_rast_q90,out_rast_int80)
  names(out_rast) <- c("pred_lai_q50","pred_lai_q10","pred_lai_q90","pred_lai_int80")
  crs(out_rast) <- crs(test_preds_1m)
  plot(out_rast[[1]],range=c(0,8),main=outname)
  writeRaster(out_rast, file=outname, overwrite=TRUE)
}

## Predicts CWC using the best Bayesian model.
site_preds <- list.files("./data/dronedata_predictors/",pattern="1m.tif$",full.names=TRUE)
#for(i in 1:4){
for(i in 1:length(site_preds)){
  
  outname <- gsub("/data/ortho/", "/output/cwc_predictions/", image_paths[i,]$ortho_path)
  outname <- gsub("_ortho.tif", "_predicted_CWC_1m_mod5.tif", outname)
  
  #if(file.exists(outname)){
  #  print(paste(outname, "exists, skipping..."))
  #  next
  #}else{
  #  print(paste("Creating model predictions for", outname))
  #}
  
  test_preds_1m <- rast(site_preds[i])
  
  test_preds_1m_df <- as.data.frame(test_preds_1m,xy=TRUE)
  test_preds_1m_df$cloudy_yn <- image_paths[i,]$cloudy_yn
  test_preds_1m_df$Site.x <- image_paths[i,]$Site
  test_preds_1m_df$visit_code <- image_paths[i,]$visit_code
  test_preds_1m_df_comp <- test_preds_1m_df[complete.cases(test_preds_1m_df),]
  
  test_cwc_pred <- posterior_predict(trait_mvmod1,m=2, newdata=test_preds_1m_df_comp, draws=100,
                                     re.form=~(1+plot_gcc_mean | Site.x)) - 50
  test_cwc_quants <- t(apply(test_cwc_pred,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975),na.rm=TRUE))
  
  ## Sets predictions for very low gcc to zero.
  test_cwc_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,4] <- 0 
  test_cwc_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,3] <- 0 
  test_cwc_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,2] <- 0 
  test_cwc_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,1] <- 0 
  
  out_df_q50 <- data.frame(x=test_preds_1m_df_comp$x, 
                           y=test_preds_1m_df_comp$y,
                           pred_cwc_q50=test_cwc_quants[,4])
  
  out_rast_q50 <- rast(out_df_q50)
  
  out_df_q10 <- data.frame(x=test_preds_1m_df_comp$x, 
                           y=test_preds_1m_df_comp$y,
                           pred_cwc_q10=test_cwc_quants[,2])
  out_rast_q10 <- rast(out_df_q10)
  
  out_df_q90 <- data.frame(x=test_preds_1m_df_comp$x, 
                           y=test_preds_1m_df_comp$y,
                           pred_cwc_q90=test_cwc_quants[,6])
  out_rast_q90 <- rast(out_df_q90)
  out_rast_int80 <- out_rast_q90 - out_rast_q10
  names(out_rast_int80) <- "pred_cwc_int80"
  
  out_rast <- c(out_rast_q50,out_rast_q10,out_rast_q90,out_rast_int80)
  names(out_rast) <- c("pred_cwc_q50","pred_cwc_q10","pred_cwc_q90","pred_cwc_int80")
  crs(out_rast) <- crs(test_preds_1m)
  plot(out_rast[[1]],range=c(0,5000),main=outname)
  writeRaster(out_rast, file=outname, overwrite=TRUE)
}

site_preds <- list.files("./data/dronedata_predictors/",pattern="1m.tif$",full.names=TRUE)
#for(i in 1:4){
for(i in 1:length(site_preds)){
  
  outname <- gsub("/data/ortho/", "/output/biomass_predictions/", image_paths[i,]$ortho_path)
  outname <- gsub("_ortho.tif", "_predicted_biomass_1m_mod5.tif", outname)
  
  #if(file.exists(outname)){
  #  print(paste(outname, "exists, skipping..."))
  #  next
  #}else{
  #  print(paste("Creating model predictions for", outname))
  #}
  
  test_preds_1m <- rast(site_preds[i])
  
  test_preds_1m_df <- as.data.frame(test_preds_1m,xy=TRUE)
  test_preds_1m_df$cloudy_yn <- image_paths[i,]$cloudy_yn
  test_preds_1m_df$Site.x <- image_paths[i,]$Site
  test_preds_1m_df$visit_code <- image_paths[i,]$visit_code
  test_preds_1m_df_comp <- test_preds_1m_df[complete.cases(test_preds_1m_df),]
  
  test_mass_pred <- posterior_predict(trait_mvmod1,m=3, newdata=test_preds_1m_df_comp, draws=100,
                                      re.form=~(1+plot_gcc_mean | Site.x)) - 50
  test_mass_quants <- t(apply(test_mass_pred,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975),na.rm=TRUE))
  
  ## Sets predictions for very low gcc to zero.
  test_mass_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,4] <- 0 
  test_mass_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,3] <- 0 
  test_mass_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,2] <- 0 
  test_mass_quants[test_preds_1m_df_comp$plot_gcc_mean < 0.345,1] <- 0 
  
  out_df_q50 <- data.frame(x=test_preds_1m_df_comp$x, 
                           y=test_preds_1m_df_comp$y,
                           pred_mass_q50=test_mass_quants[,4])
  
  out_rast_q50 <- rast(out_df_q50)
  
  out_df_q10 <- data.frame(x=test_preds_1m_df_comp$x, 
                           y=test_preds_1m_df_comp$y,
                           pred_mass_q10=test_mass_quants[,2])
  out_rast_q10 <- rast(out_df_q10)
  
  out_df_q90 <- data.frame(x=test_preds_1m_df_comp$x, 
                           y=test_preds_1m_df_comp$y,
                           pred_mass_q90=test_mass_quants[,6])
  out_rast_q90 <- rast(out_df_q90)
  out_rast_int80 <- out_rast_q90 - out_rast_q10
  names(out_rast_int80) <- "pred_mass_int80"
  
  out_rast <- c(out_rast_q50,out_rast_q10,out_rast_q90,out_rast_int80)
  names(out_rast) <- c("pred_mass_q50","pred_mass_q10","pred_mass_q90","pred_mass_int80")
  crs(out_rast) <- crs(test_preds_1m)
  plot(out_rast[[1]],range=c(0,800),main=outname)
  writeRaster(out_rast, file=outname, overwrite=TRUE)
}

