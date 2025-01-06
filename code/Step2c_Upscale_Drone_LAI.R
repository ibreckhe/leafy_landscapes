## Script to filter and upscale drone LAI predictions to the HLS grid.
## Author: Ian Breckheimer
## Updated: 1 January 2025

## This script should be run after drone-based LAI models are finalized. To see how this script fits into the larger analysis, see the documentation here: https://app.affine.pro/workspace/e094b233-4dce-4599-bc85-a12ab329bda3/hakWMNfZ-ZXSIhQycXIYx?mode=page

## Set up workspace.
library(rstanarm)
library(sf)
library(readr)
library(dplyr)
library(terra)
library(ggplot2)
#setwd("C:\\Users/IanBreckheimer/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/")
setwd("~/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/")

## Computes site-level masks for the LAI predictions.
lai_preds <- list.files("./output/lai_predictions/",pattern="mod5.tif$",full.names=TRUE)

hls_mask_30m <- rast("./data/HLS/HLS_drone_mask_30m.tif")

unc_mask_fun <- function(mean_rasts,unc_rasts,template_rast,
                         unc_thresh=3.0,mean_thresh=15){
  mean_rasts <- rast(lapply(mean_rasts,crop,template_rast,extend=TRUE))
  zero_rasts <- !(mean_rasts==0)
  unc_rasts <- rast(lapply(unc_rasts,crop,template_rast,extend=TRUE))
  unc_rel_rasts <- unc_rasts / mean_rasts
  rel_unc_masks <- unc_rel_rasts > unc_thresh | mean_rasts > mean_thresh
  rel_unc_masks <- rel_unc_masks * zero_rasts
  rel_unc_mask <- sum(rel_unc_masks)
  rel_unc_mask_out <- app(rel_unc_mask,fun=function(x){ifelse(x==0,1,NA)})
  return(rel_unc_mask_out)
}

jf_pred_rasts <- lapply(lai_preds[1:4],rast,lyrs=1)
jf_unc_rasts <- lapply(lai_preds[1:4],rast,lyrs=4)
jf_rel_unc_mask <- unc_mask_fun(jf_pred_rasts,jf_unc_rasts,jf_pred_rasts[[1]],
                                unc_thresh=2.0,mean_thresh=10)
plot(jf_rel_unc_mask)

lt_pred_rasts <- lapply(lai_preds[5:8],rast,lyrs=1)
lt_unc_rasts <- lapply(lai_preds[5:8],rast,lyrs=4)
lt_rel_unc_mask <- unc_mask_fun(lt_pred_rasts,lt_unc_rasts,lt_pred_rasts[[1]],
                                unc_thresh=3.0,mean_thresh=15)
plot(lt_rel_unc_mask)

pb_pred_rasts <- lapply(lai_preds[9:12],rast,lyrs=1)
pb_unc_rasts <- lapply(lai_preds[9:12],rast,lyrs=4)
pb_rel_unc_mask <- unc_mask_fun(pb_pred_rasts,pb_unc_rasts,pb_pred_rasts[[1]],
                                unc_thresh=3.0,mean_thresh=15)
plot(pb_rel_unc_mask)

rg_pred_rasts <- lapply(lai_preds[13:16],rast,lyrs=1)
rg_unc_rasts <- lapply(lai_preds[13:16],rast,lyrs=4)
rg_rel_unc_mask <- unc_mask_fun(rg_pred_rasts,rg_unc_rasts,rg_pred_rasts[[1]],
                                unc_thresh=3.0,mean_thresh=15)
plot(rg_rel_unc_mask)

sr_pred_rasts <- lapply(lai_preds[17:20],rast,lyrs=1)
sr_unc_rasts <- lapply(lai_preds[17:20],rast,lyrs=4)
sr_rel_unc_mask <- unc_mask_fun(sr_pred_rasts,sr_unc_rasts,sr_pred_rasts[[1]],
                                unc_thresh=2.0,mean_thresh=10)
plot(sr_rel_unc_mask)

sn_pred_rasts <- lapply(lai_preds[21:24],rast,lyrs=1)
sn_unc_rasts <- lapply(lai_preds[21:24],rast,lyrs=4)
sn_rel_unc_mask <- unc_mask_fun(sn_pred_rasts,sn_unc_rasts,sn_pred_rasts[[1]],
                                unc_thresh=3.0,mean_thresh=15)
plot(sn_rel_unc_mask)

sn_pred_rasts <- lapply(lai_preds[21:24],rast,lyrs=1)
sn_unc_rasts <- lapply(lai_preds[21:24],rast,lyrs=4)
sn_rel_unc_mask <- unc_mask_fun(sn_pred_rasts,sn_unc_rasts,sn_pred_rasts[[1]],
                                unc_thresh=3.0,mean_thresh=15)
plot(sn_rel_unc_mask)

sb_pred_rasts <- lapply(lai_preds[25:28],rast,lyrs=1)
sb_unc_rasts <- lapply(lai_preds[25:28],rast,lyrs=4)
sb_rel_unc_mask <- unc_mask_fun(sb_pred_rasts,sb_unc_rasts,sb_pred_rasts[[1]],
                                unc_thresh=3.0,mean_thresh=15)
plot(sb_rel_unc_mask)

sh_pred_rasts <- lapply(lai_preds[29:32],rast,lyrs=1)
sh_unc_rasts <- lapply(lai_preds[29:32],rast,lyrs=4)
sh_rel_unc_mask <- unc_mask_fun(sh_pred_rasts,sh_unc_rasts,sh_pred_rasts[[1]],
                                unc_thresh=3.0,mean_thresh=15)
plot(sh_rel_unc_mask)

mask_list <- list(jf_rel_unc_mask,lt_rel_unc_mask,pb_rel_unc_mask,rg_rel_unc_mask,
                  sr_rel_unc_mask,sn_rel_unc_mask,sb_rel_unc_mask,sh_rel_unc_mask)
mask_zero <- lapply(mask_list,function(x){subst(x,NA,0)})

hls_mask_fun <- function(drone_mask,hls_mask){
  hls_crop <- crop(hls_mask,drone_mask)
  hls_proj <- project(drone_mask,hls_crop,method="average")
  hls_proj_mask <- hls_proj >= 0.95
  return(hls_proj_mask)
}

mask_zero_hls <- lapply(mask_zero,hls_mask_fun,hls_mask=hls_mask_30m)
mask_na_hls <- lapply(mask_zero_hls,function(x){subst(x,0,NA)})
mask_hls_full <- sum(rast(lapply(mask_na_hls,function(x){extend(x,hls_mask_30m)})),na.rm=TRUE)
writeRaster(mask_hls_full,"./data/HLS/HLS_drone_mask_30m.tif",overwrite=TRUE)
mask_hls_df <- as.data.frame(mask_hls_full,xy=TRUE)

## Looks at distribution of mean predictions and uncertainty.
hls_mean_fun <- function(drone_pred,hls_mask){
  hls_crop <- crop(hls_mask,drone_pred)
  hls_proj <- project(drone_pred,hls_crop,method="average")
  hls_out <- mask(hls_proj,hls_crop)
  return(hls_out)
}

jf_hls_means <- lapply(jf_pred_rasts,hls_mean_fun,hls_mask=mask_hls_full)
names(jf_hls_means) <- lai_preds[1:4]
jf_hls_df <- lapply(jf_hls_means,function(x){as.data.frame(x,xy=TRUE)})

lt_hls_means <- lapply(lt_pred_rasts,hls_mean_fun,hls_mask=mask_hls_full)
lt_hls_df <- lapply(lt_hls_means,function(x){as.data.frame(x,xy=TRUE)})

pb_hls_means <- lapply(pb_pred_rasts,hls_mean_fun,hls_mask=mask_hls_full)
pb_hls_df <- lapply(pb_hls_means,function(x){as.data.frame(x,xy=TRUE)})

rg_hls_means <- lapply(rg_pred_rasts,hls_mean_fun,hls_mask=mask_hls_full)
rg_hls_df <- lapply(rg_hls_means,function(x){as.data.frame(x,xy=TRUE)})

sr_hls_means <- lapply(sr_pred_rasts,hls_mean_fun,hls_mask=mask_hls_full)
sr_hls_df <- lapply(sr_hls_means,function(x){as.data.frame(x,xy=TRUE)})

sn_hls_means <- lapply(sn_pred_rasts,hls_mean_fun,hls_mask=mask_hls_full)
sn_hls_df <- lapply(sn_hls_means,function(x){as.data.frame(x,xy=TRUE)})

sb_hls_means <- lapply(sb_pred_rasts,hls_mean_fun,hls_mask=mask_hls_full)
sb_hls_df <- lapply(sb_hls_means,function(x){as.data.frame(x,xy=TRUE)})

sh_hls_means <- lapply(sh_pred_rasts,hls_mean_fun,hls_mask=mask_hls_full)
sh_hls_df <- lapply(sh_hls_means,function(x){as.data.frame(x,xy=TRUE)})

hls_pred_df <- c(jf_hls_df,lt_hls_df,pb_hls_df,rg_hls_df,sr_hls_df,sn_hls_df,sb_hls_df,sh_hls_df)
names(hls_pred_df) <- lai_preds

hls_pred_df_comb <- bind_rows(hls_pred_df,.id = "id")
hls_pred_df_comb$ortho <- gsub("./output/lai_predictions//","",hls_pred_df_comb$id)
hls_pred_df_comb$Site <- stringr::str_split_fixed(hls_pred_df_comb$ortho,"_",2)[,1]
ortho_year <- stringr::str_split_fixed(hls_pred_df_comb$ortho,"_",3)[,2]
ortho_month<- stringr::str_split_fixed(hls_pred_df_comb$ortho,"_",4)[,3]
ortho_day <- stringr::str_split_fixed(hls_pred_df_comb$ortho,"_",5)[,4]
hls_pred_df_comb$Date <- as.Date(paste(ortho_year,ortho_month,ortho_day,sep="-"))
hls_pred_df_comb$DOY <- as.numeric(format(hls_pred_df_comb$Date,"%j"))
hls_pred_df_comb$PointID <- paste0(hls_pred_df_comb$x,"E_", hls_pred_df_comb$y,"N")

## Preliminary plots.
library(ggplot2)
p1 <- ggplot(hls_pred_df_comb)+
  ggtitle("Predicted LAI (posterior median) at 30m Resolution") +
  geom_point(aes(x=Date,y=pred_lai_q50,color=Site),
             size=0.4,position=position_jitter(width=2)) +
  geom_smooth(aes(x=Date,y=pred_lai_q50,color=Site,group=PointID),se=FALSE,lwd=0.1,alpha=0.05) +
  #scale_x_continuous("Day of Year") +
  scale_y_continuous("Predicted Leaf Area Index (unitless)") +
  facet_wrap(~Site,scales="fixed",ncol=4) +
  theme_bw()+
  theme(legend.position="none")

pdf("./figs/upscaled_lai_q50_predictions_30m_v5.pdf",width=10,height=7)
print(p1)
dev.off()

## Writes predictions to disk.
write_csv(hls_pred_df_comb,"./output/upscaled_lai_q50_predictions_30m_mod5.csv")
