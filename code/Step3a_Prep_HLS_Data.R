## Script to prepare 2024 HLS satellite data for analysis
## Author: Ian Breckheimer
## Updated: 1 January 2025

## This script should be run after we are done preparing upscaled drone data. To see how this script fits into the larger analysis, see the documentation here: https://app.affine.pro/workspace/e094b233-4dce-4599-bc85-a12ab329bda3/SRJhHNf7zqPPvk6edPJZ4?mode=page

## Sets up workspace
library(tidyr)
library(dplyr)
library(sf)
library(terra)

## Set working directory
setwd("~/Library/CloudStorage/OneDrive-RMBL/Documents\ -\ Research\ -\ Spatial\ Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI")
#setwd("C:/Users/IanBreckheimer/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/")

## Reads in template extent.
hls_ext <- rast("./data/HLS/HLS_mask_UER.tif")

## Gets lat and lon bounds from extent.
hls_ext <- project(ext(hls_ext), from="EPSG:32613", to="EPSG:4326")

## 
# define a region of interest
roi <- c(
  lon_min = hls_ext[1], lat_min = hls_ext[3],
  lon_max = hls_ext[2], lat_max = hls_ext[4]
)
names(roi) <- c("lon_min", "lat_min", "lon_max", "lat_max")

### Gets metadata stats for scenes.
hls_l30_meta <- read.csv("./data/HLS/HLS_L30_2024/metadata/HLSL30-020-Fmask-Statistics-QA.csv")

### Lists tif files for HLS data.
hls_files <- list.files("./data/HLS/HLS_L30_2024/", pattern = ".tif$", full.names = TRUE)
b01_rast <- rast(hls_files[grep("B01", hls_files)])
b02_rast <- rast(hls_files[grep("B02", hls_files)])
b03_rast <- rast(hls_files[grep("B03", hls_files)])
b04_rast <- rast(hls_files[grep("B04", hls_files)])
b05_rast <- rast(hls_files[grep("B05", hls_files)])
b06_rast <- rast(hls_files[grep("B06", hls_files)])
b07_rast <- rast(hls_files[grep("B07", hls_files)])
b09_rast <- rast(hls_files[grep("B09", hls_files)])
b10_rast <- rast(hls_files[grep("B10", hls_files)])
b11_rast <- rast(hls_files[grep("B11", hls_files)])
fmask_rast <- rast(hls_files[grep("Fmask", hls_files)])

## combines clear scenes into a list.
#ls_clear_scenes <- c(5,10,14,19,23,26,27,31,32,33,37,47,49,54)
# ls_list <- list(b01_rast, b02_rast, 
#                 b03_rast, b04_rast, 
#                 b05_rast, b06_rast, 
#                 b07_rast, b10_rast, 
#                 b11_rast, fmask_rast)
ls_list <- list(b01_rast, b02_rast, 
                b03_rast, b04_rast, 
                b05_rast, b06_rast, 
                b07_rast, b10_rast, 
                b11_rast, fmask_rast)

## Extracts data from filtered locations of drone data pixels.
drone_locs <- rast("./data/HLS/HLS_drone_mask_30m.tif")
drone_loc_df <- as.data.frame(drone_locs,xy=TRUE)
drone_loc_mat <- as.matrix(drone_loc_df[,1:2])

extract_fun <- function(x){
  ext_data <- terra::extract(x,y=drone_loc_mat)
  out_wide <- data.frame(x=drone_loc_df[,"x"],
                         y=drone_loc_df[,"y"],
                         ext_data)
  out_long <- pivot_longer(out_wide, cols=contains("doy"), names_to="source", values_to="value")
  return(out_long)
}

ls_extracted <- lapply(ls_list, extract_fun)
ls_bind <- bind_rows(ls_extracted)
ls_bind$Band <- stringr::str_split_fixed(ls_bind$source, "_", 4)[,2]
ls_bind$Date <- as.Date(gsub("doy","",stringr::str_split_fixed(ls_bind$source, "_", 4)[,3]),
                        format="%Y%j")
ls_bind_wide <- pivot_wider(ls_bind[,-3], names_from=c("Band"),values_from="value")
ls_bind_wide$qc_clear <- apply(matrix(as.numeric(sapply(ls_bind_wide$Fmask, intToBits,simplify=TRUE)),ncol=32,byrow=TRUE)[,1:5],
                               FUN=sum, MARGIN=1) == 0
ls_bind_wide$PointID <- paste0(ls_bind_wide$x,"E_",ls_bind_wide$y,"N")
ls_bind_filt <- filter(ls_bind_wide,qc_clear %in% c(TRUE))

LS_bands <- c("Aerosol"="B01",
              "Blue"="B02",
              "Green"="B03",
              "Red"="B04",
              "NIR"="B05",
              "SWIR1"="B06",
              "SWIR2"="B07")
LS_keep <- c("x","y","Date","Fmask","PointID",LS_bands)
ls_band_filt <- ls_bind_filt[,LS_keep]
colnames(ls_band_filt) <- c("x","y","Date","Fmask","PointID","Aerosol","Blue","Green","Red","NIR","SWIR1","SWIR2")

ls_band_filt$NDVI <- (ls_band_filt$NIR - ls_band_filt$Red) / (ls_band_filt$NIR + ls_band_filt$Red)
ls_band_filt$EVI2 <- 2.4 * (ls_band_filt$NIR - ls_band_filt$Red) / (ls_band_filt$NIR + ls_band_filt$Red + 1)
ls_band_filt$NDWI1 <- (ls_band_filt$SWIR1 - ls_band_filt$NIR) / (ls_band_filt$SWIR1 + ls_band_filt$NIR)
ls_band_filt$ANDWI <- (ls_band_filt$Blue + ls_band_filt$Green + ls_band_filt$Red - ls_band_filt$NIR - ls_band_filt$SWIR1 - ls_band_filt$SWIR2) / 
  (ls_band_filt$Blue + ls_band_filt$Green + ls_band_filt$Red + ls_band_filt$NIR + ls_band_filt$SWIR1 + ls_band_filt$SWIR2)
ls_band_filt$NDGI <- (0.65 * ls_band_filt$Green + (1 - 0.65) * ls_band_filt$NIR - ls_band_filt$Red) / 
  (0.65 * ls_band_filt$Green + (1 - 0.65) * ls_band_filt$NIR + ls_band_filt$Red)

ls_bind_long <- pivot_longer(ls_band_filt,cols=c("Aerosol","Blue","Green","Red","NIR","SWIR1","SWIR2",
                                                 "NDVI","EVI2","NDWI1","ANDWI","NDGI"),
                             names_to="Band",values_to="value")

library(ggplot2)
ggplot(ls_bind_long,aes(x=Date,y=value,color=Band))+
  geom_point(position=position_jitter(width=2,height=0.1),pch=".")+
  #scale_y_continuous(limits=c(0,1))+
  facet_wrap(facets=~Band,scales="free_y")+
  theme_bw()

## Sentinel 2 data.
### Lists tif files for HLS data.
extend_fun <- function(x,ext){
  x <- rast(x)
  x <- extend(x,ext)
  return(x)
}

hls_s30_files <- list.files("./data/HLS/HLS_S30_2024/", pattern = ".tif$", full.names = TRUE)
template_ext <- ext(rast(hls_s30_files[1]))

b01_rast <- rast(lapply(hls_s30_files[grep("B01", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b02_rast <- rast(lapply(hls_s30_files[grep("B02", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b03_rast <- rast(lapply(hls_s30_files[grep("B03", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b04_rast <- rast(lapply(hls_s30_files[grep("B04", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b05_rast <- rast(lapply(hls_s30_files[grep("B05", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b06_rast <- rast(lapply(hls_s30_files[grep("B06", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b07_rast <- rast(lapply(hls_s30_files[grep("B07", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b08_rast <- rast(lapply(hls_s30_files[grep("B08", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b08a_rast <- rast(lapply(hls_s30_files[grep("B8A", hls_s30_files)], 
                         FUN=extend_fun, ext=template_ext))
b09_rast <- rast(lapply(hls_s30_files[grep("B09", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b10_rast <- rast(lapply(hls_s30_files[grep("B10", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b11_rast <- rast(lapply(hls_s30_files[grep("B11", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))
b12_rast <- rast(lapply(hls_s30_files[grep("B12", hls_s30_files)], 
                        FUN=extend_fun, ext=template_ext))

fmask_rast <- rast(lapply(hls_s30_files[grep("Fmask", hls_s30_files)], 
                          FUN=extend_fun, ext=template_ext))

#s2_clear_scenes <- c(8,14,24,32,36,40,47,49,53,59,76)
# s2_list <- list(b01_rast[[s2_clear_scenes]], b02_rast[[s2_clear_scenes]], 
#                 b03_rast[[s2_clear_scenes]], b04_rast[[s2_clear_scenes]], 
#                 b05_rast[[s2_clear_scenes]], b06_rast[[s2_clear_scenes]], 
#                 b07_rast[[s2_clear_scenes]], b08_rast[[s2_clear_scenes]], 
#                 b08a_rast[[s2_clear_scenes]], b10_rast[[s2_clear_scenes]], 
#                 b11_rast[[s2_clear_scenes]], b12_rast[[s2_clear_scenes]], 
#                 fmask_rast[[s2_clear_scenes]])

s2_list <- list(b01_rast, b02_rast, 
                b03_rast, b04_rast, 
                b05_rast, b06_rast, 
                b07_rast, b08_rast, 
                b08a_rast, b10_rast, 
                b11_rast, b12_rast, 
                fmask_rast)

extract_fun <- function(x){
  ext_data <- terra::extract(x,y=drone_loc_mat)
  out_wide <- data.frame(x=drone_loc_df[,"x"],
                         y=drone_loc_df[,"y"],
                         ext_data)
  out_long <- pivot_longer(out_wide, cols=contains("doy"), names_to="source", values_to="value")
  return(out_long)
}

s2_extracted <- lapply(s2_list, extract_fun)
s2_bind <- bind_rows(s2_extracted)
s2_bind$Band <- stringr::str_split_fixed(s2_bind$source, "_", 4)[,2]
s2_bind$Date <- as.Date(gsub("doy","",stringr::str_split_fixed(s2_bind$source, "_", 4)[,3]),
                        format="%Y%j")
s2_bind_wide <- pivot_wider(s2_bind[,-3], names_from=c("Band"),values_from="value")
s2_bind_wide$qc_clear <- apply(matrix(as.numeric(sapply(s2_bind_wide$Fmask, intToBits,simplify=TRUE)),ncol=32,byrow=TRUE)[,1:5],
                               FUN=sum, MARGIN=1) == 0
s2_bind_wide$PointID <- paste0(s2_bind_wide$x,"E_",s2_bind_wide$y,"N")
s2_bind_filt <- filter(s2_bind_wide,qc_clear %in% c(TRUE))

s2_bands <- c("Aerosol"="B01",
              "Blue"="B02",
              "Green"="B03",
              "Red"="B04",
              "NIR"="B08",
              "SWIR1"="B11",
              "SWIR2"="B12")
s2_keep <- c("x","y","Date","Fmask","PointID",s2_bands)
s2_band_filt <- s2_bind_filt[,s2_keep]
colnames(s2_band_filt) <- c("x","y","Date","Fmask","PointID","Aerosol","Blue","Green","Red","NIR","SWIR1","SWIR2")

s2_band_filt$NDVI <- (s2_band_filt$NIR - s2_band_filt$Red) / (s2_band_filt$NIR + s2_band_filt$Red)
s2_band_filt$EVI2 <- 2.4 * (s2_band_filt$NIR - s2_band_filt$Red) / (s2_band_filt$NIR + s2_band_filt$Red + 1)
s2_band_filt$NDWI1 <- (s2_band_filt$SWIR1 - s2_band_filt$NIR) / (s2_band_filt$SWIR1 + s2_band_filt$NIR)
s2_band_filt$ANDWI <- (s2_band_filt$Blue + s2_band_filt$Green + s2_band_filt$Red - s2_band_filt$NIR - s2_band_filt$SWIR1 - s2_band_filt$SWIR2) / 
  (s2_band_filt$Blue + s2_band_filt$Green + s2_band_filt$Red + s2_band_filt$NIR + s2_band_filt$SWIR1 + s2_band_filt$SWIR2)
s2_band_filt$NDGI <- (0.65 * s2_band_filt$Green + (1 - 0.65) * s2_band_filt$NIR - s2_band_filt$Red) / 
  (0.65 * s2_band_filt$Green + (1 - 0.65) * s2_band_filt$NIR + s2_band_filt$Red)

s2_bind_long <- pivot_longer(s2_band_filt,cols=c("Aerosol","Blue","Green","Red","NIR","SWIR1","SWIR2",
                                                 "NDVI","EVI2","NDWI1","ANDWI","NDGI"),
                             names_to="Band",values_to="value")

ggplot(s2_bind_long,aes(x=Date,y=value,color=Band))+
  geom_point(position=position_jitter(width=1,height=0.1),pch=".")+
  #scale_y_continuous(limits=c(0,1))+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-09-01")))+
  facet_wrap(facets=~Band,scales="free_y")+
  theme_bw()

## Plots both Landsat and Sentinel 2 NDVI.
ls_ndvi <- ls_band_filt[,c("x", "y", "PointID", "Date","Fmask",
                           "Aerosol","Blue","Green","Red","NIR","SWIR1","SWIR2",
                           "NDVI", "EVI2", "NDWI1","ANDWI","NDGI")]
s2_ndvi <- s2_band_filt[,c("x", "y", "PointID", "Date", "Fmask",
                           "Aerosol","Blue","Green","Red","NIR","SWIR1","SWIR2",
                           "NDVI", "EVI2", "NDWI1", "ANDWI","NDGI")]

ls_ndvi$Source <- "HLS-L30"
s2_ndvi$Source <- "HLS-S30"

ndvi_bind <- bind_rows(ls_ndvi, s2_ndvi)

ndvi_sp <- st_as_sf(ndvi_bind, coords=c("x","y"), crs=st_crs(32613))
polys <- st_read("./data/Leafy_landscapes_drone_boundaries.geojson")
points_vect <- vect(ndvi_sp)
polys_vect <- vect(polys)
ndvi_bind$Site <- gsub("_Clip","",terra::extract(polys_vect,points_vect)$layer)

sample_plots <- unique(ndvi_bind$PointID)[sample(1:length(unique(ndvi_bind$PointID)), 150)]
ndvi_samp <- ndvi_bind[ndvi_bind$PointID %in% sample_plots,]

ggplot(filter(ndvi_bind))+
  geom_point(aes(x=Date,y=NDGI,color=Site),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=NDGI,color=Site),method="loess",se=FALSE,span=0.6)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  #scale_y_continuous(limits=c(0,1))+
  #facet_wrap(facets=~PointID)+
  theme_bw()

ggplot(filter(ndvi_bind))+
         geom_point(aes(x=Date,y=EVI2,color=Site),
                    position=position_jitter(width=1),size=0.8)+
         geom_smooth(aes(x=Date,y=EVI2,color=Site),method="loess",se=FALSE,span=0.6)+
         scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
         scale_y_continuous(limits=c(0,1))+
         #facet_wrap(facets=~PointID)+
         theme_bw()

ggplot(filter(ndvi_bind))+
  geom_point(aes(x=Date,y=ANDWI,color=Site),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=ANDWI,color=Site),method="loess",se=FALSE,span=0.6)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  scale_y_continuous(limits=c(-1,1))+
  #facet_wrap(facets=~PointID)+
  theme_bw()

ggplot(filter(ndvi_samp,Site=="StrandHill"))+
  geom_point(aes(x=Date,y=NDGI,color=Source),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=NDGI,group=PointID),method="loess",se=FALSE,span=0.5,lwd=0.5)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  scale_y_continuous(limits=c(0,1))+
  facet_wrap(facets=~PointID)+
  theme_bw()+
  theme(legend.position="none")

ggplot(filter(ndvi_samp,Site=="LupineTrail"))+
  geom_point(aes(x=Date,y=NDGI,color=Source),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=NDGI),method="loess",se=FALSE,span=0.5,lwd=0.5)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  scale_y_continuous(limits=c(0,1))+
  facet_wrap(facets=~PointID)+
  theme_bw()+
  theme(legend.position="none")

ggplot(filter(ndvi_samp,Site=="SlideRockRoad"))+
  geom_point(aes(x=Date,y=NDGI,color=Source),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=NDGI,color=PointID),method="loess",se=FALSE,span=0.5,lwd=0.5)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  scale_y_continuous(limits=c(0,1))+
  facet_wrap(facets=~PointID)+
  theme_bw()+
  theme(legend.position="none")

ggplot(filter(ndvi_samp,Site=="SnodgrassTrailhead"))+
  geom_point(aes(x=Date,y=NDGI,color=Source),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=NDGI,color=PointID),method="loess",se=FALSE,span=0.5,lwd=0.5)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  scale_y_continuous(limits=c(0,1))+
  facet_wrap(facets=~PointID)+
  theme_bw()+
  theme(legend.position="none")

ggplot(filter(ndvi_samp,Site=="JuddFalls"))+
  geom_point(aes(x=Date,y=NDGI,color=Source),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=NDGI,color=PointID),method="loess",se=FALSE,span=0.5,lwd=0.5)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  scale_y_continuous(limits=c(0,1))+
  facet_wrap(facets=~PointID)+
  theme_bw()+
  theme(legend.position="none")

ggplot(filter(ndvi_samp,Site=="RustlersGulch"))+
  geom_point(aes(x=Date,y=NDGI,color=Source),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=NDGI,color=PointID),method="loess",se=FALSE,span=0.5,lwd=0.5)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  scale_y_continuous(limits=c(0,1))+
  facet_wrap(facets=~PointID)+
  theme_bw()+
  theme(legend.position="none")

ggplot(filter(ndvi_samp,Site=="SouthBaldy"))+
  geom_point(aes(x=Date,y=NDGI,color=Source),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=NDGI,color=PointID),method="loess",se=FALSE,span=0.5,lwd=0.5)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  scale_y_continuous(limits=c(0,1))+
  facet_wrap(facets=~PointID)+
  theme_bw()+
  theme(legend.position="none")

ggplot(filter(ndvi_samp,Site=="ParadiseBasin"))+
  geom_point(aes(x=Date,y=NDGI,color=Source),
             position=position_jitter(width=1),size=0.8)+
  geom_smooth(aes(x=Date,y=NDGI),method="loess",se=FALSE,span=0.8,span=0.5,lwd=0.5)+
  scale_x_date(limits=c(as.Date("2024-05-01"),as.Date("2024-10-01")))+
  scale_y_continuous(limits=c(0,1))+
  facet_wrap(facets=~PointID)+
  theme_bw()+
  theme(legend.position="none")

## Brings in sample dates.
drone_df <- read.csv("./output/upscaled_lai_q50_predictions_30m_mod5.csv")
drone_df$Date <- as.Date(drone_df$Date, format="%Y-%m-%d")
drone_dates <- unique(as.Date(drone_df$Date, format="%Y-%m-%d"))

## Function for interpolating landsat observations.

ndvi_long <- pivot_longer(ndvi_bind,cols=c("Aerosol","Blue","Green","Red","NIR","SWIR1","SWIR2",
                                              "NDVI","EVI2","NDWI1","ANDWI","NDGI"),
                          names_to="Band",values_to="value")
  
test_df <- filter(ndvi_long,PointID==sample_plots[1])
test_dates <- unique(test_df$Date)
out_dates <- seq(min(test_dates),max(test_dates),by=7)

loess_fun <- function(x){
  loess_mod <- loess(x$value ~ as.numeric(x$Date), span=0.5)
  loess_pred <- predict(loess_mod, newdata=as.numeric(drone_dates))
  loess_out <- data.frame(Date=drone_dates, value=loess_pred)
  return(loess_out)
}
test_loess <- test_df %>% group_by(PointID,Band) %>%
                          group_modify(~loess_fun(.x))

ndvi_long_interp <- ndvi_long %>% group_by(PointID,Band) %>%
                                  group_modify(~loess_fun(.x))
ndvi_wide_interp <- pivot_wider(ndvi_long_interp,names_from="Band",values_from="value")

drone_df_join <- left_join(drone_df,ndvi_wide_interp,by=c("Date"="Date","PointID"="PointID"))
write.csv(drone_df_join,"./output/upscaled_lai_q50_predictions_30m_HLS_interp.csv",row.names=FALSE)
#drone_df_join <- read.csv("./output/upscaled_lai_q50_predictions_30m_HLS_interp.csv")

## Plotting the interpolated values.
ggplot(drone_df_join)+
  geom_point(aes(x=NDGI, y=pred_lai_q50, color=Site),size=0.3)+
  geom_smooth(aes(x=NDGI, y=pred_lai_q50), method="gam",se=FALSE,span=0.8)+
  #facet_wrap(facets=~as.factor(Site)+as.factor(Date))+
  theme_bw()

## Testing the upscaling using 2024 data.

# ## Preliminary upscaled models
# library(mgcv)
# gam1 <- gamm(pred_lai_q50 ~ s(NDGI),random=list(ortho=~1), data=drone_df_join)
# summary(gam1$lme)
# summary(gam1$gam)
# plot(gam1$gam)
# saveRDS(gam1,"./output/HLS_LAI_model_gam1.rds")
# 
# #cloud_codes <- c(66, 70, 74 ,78, 130, 134, 138, 142, 194, 198, 202, 206)
# 
# ndgi_fun <- function(x){
#   ndgi <- ifelse(x[4] %in% c(64, 68, 80, 84, 96, 100, 112, 116, 128, 
#                                132, 144, 148, 160, 164, 176, 180,
#                                192, 196, 208, 212, 224, 228, 240, 244), 
#          (0.65 * x[1] + (1 - 0.65) * x[3] - x[2]) / 
#          (0.65 * x[1] + (1 - 0.65) * x[3] + x[2]), 
#          NA)
#   ndgi[ndgi < -1] <- -1
#   ndgi[ndgi > 1.2] <- 1.2
#   return(ndgi)
# }
# 
# 
# for(i in 1:nlyr(ls_list[[1]])){
#   print(paste("Predicting LAI for",i,"of",nlyr(ls_list[[1]])))
#   test_HLS <- c(ls_list[[3]][[i]],
#                 ls_list[[4]][[i]],
#                 ls_list[[5]][[i]],
#                 ls_list[[10]][[i]])
#   test_ndgi <- app(test_HLS,ndgi_fun,cores=4)
#   #plot(test_ndgi,range=c(-1,1))
#   names(test_ndgi) <- "NDGI"
#   test_lai <- terra::predict(test_ndgi, model=gam1$gam)
#   names(test_lai) <- gsub("B03","LAI",names(test_HLS)[1])
#   outname <- paste0("./output/lai_predictions_HLS/",gsub("B03","LAI",names(test_lai)),".tif")
#   plot(test_lai,range=c(0,8),main=outname)
#   writeRaster(test_lai,filename=outname,overwrite=TRUE)
# }
# 
# for(i in 1:nlyr(s2_list[[1]])){
#   print(paste("Predicting LAI for",i,"of",nlyr(s2_list[[1]])))
#   test_HLS <- c(s2_list[[3]][[i]],
#                 s2_list[[4]][[i]],
#                 s2_list[[8]][[i]],
#                 s2_list[[13]][[i]])
#   test_ndgi <- app(test_HLS,ndgi_fun,cores=4)
#   #plot(test_ndgi,range=c(-1,1))
#   names(test_ndgi) <- "NDGI"
#   test_lai <- terra::predict(test_ndgi, model=gam1$gam)
#   names(test_lai) <- gsub("B03","LAI",names(test_HLS)[1])
#   outname <- paste0("./output/lai_predictions_HLS/",gsub("B03","LAI",names(test_lai)),".tif")
#   plot(test_lai,range=c(0,8),main=outname)
#   writeRaster(test_lai,filename=outname,overwrite=TRUE)
# }
# 
# ## Interpolates the LAI predictions to generate continuous time-series.
# lai_df <- data.frame(path=list.files(path="./output/lai_predictions_HLS",pattern=".tif$",full.names=TRUE))
# lai_df$Date <- as.Date(gsub("doy","",stringr::str_split_fixed(lai_df$path,"_|\\.",8)[,7]),format="%Y%j")
# lai_df <- lai_df[order(lai_df$Date),]
# 
# lai_rast <- rast(lai_df$path)
# 
# out_dates <- seq(as.Date("2024-04-01"),as.Date("2024-10-01"),by=1)
# 
# loess_fun <- function(x,dates=as.numeric(lai_df$Date)){
#   ifelse(all(is.na(x)),return(rep(NA,length(out_dates))),
#          try({
#            loess_mod <- loess(x ~ dates, span=0.4,
#                               control = loess.control(surface = "direct"))
#            loess_pred <- predict(loess_mod, newdata=as.numeric(out_dates))
#            loess_pred[loess_pred < 0] <- 0
#            loess_pred[loess_pred > 10] <- 10
#            return(round(loess_pred * 1000))
#          }))
#   
# }
# 
# #plot(lai_rast[[49]])
# #test_ext <- draw()
# #test_lai <- crop(lai_rast,test_ext)
# 
# lai_interp <- app(lai_rast,fun=loess_fun, filename="./output/lai_predictions_HLS_interp_2024.tif",
#                    overwrite=TRUE,wopt=list(datatype="INT2U",
#                    names=paste0("LAI-",format(out_dates,"%Y-doy%j"))))
# #names(laiinterp) <- out_dates
# #plot(lai_interp,range=c(0,6000))
# #lai_mem <- lai_rast + 0
# 
# ## Computes max LAI and integrated LAI.
# pred_interval <- as.numeric(out_dates[2] - out_dates[1])
# integral_fun <- function(x,interval=pred_interval){
#   ifelse(all(is.na(x)),return(NA),
#          try({
#            integral <- sum(x/1000 * interval)
#            return(integral)
#          }))
# }
# lai_integral <- app(lai_interp,integral_fun,filename="./output/lai_predictions_HLS_integral_2024.tif",
#                      overwrite=TRUE)
# 
# ## Computes date of max LAI.
# out_doy <- as.numeric(format(out_dates,"%j"))
# peak_lai_fun <- function(x){
#   ifelse(all(is.na(x)),return(NA),
#          try({
#            max_date <- out_doy[which.max(x)]
#            return(max_date)
#          }))
# }
# lai_peak <- app(lai_interp,peak_lai_fun,filename="./output/lai_predictions_HLS_peakDOY_2024.tif",
#                  overwrite=TRUE)
# 
# amp_lai_fun <- function(x){
#   ifelse(all(is.na(x)),return(NA),
#          try({
#            amp <- quantile(x,0.95,na.rm=TRUE) - quantile(x,0.05,na.rm=TRUE)
#            return(amp)
#          }))
# }
# lai_amp <- app(lai_interp,amp_lai_fun,filename="./output/lai_predictions_HLS_LAIamp_2024.tif",
#                  overwrite=TRUE)
# 
# 
# ## Computes dates reaching quantiles of max LAI.
# gup_lai_fun <- function(x,quantile_thresh=0.25){
#   ifelse(all(is.na(x)),return(NA),
#          try({
#            x_quants <- quantile(x,probs=c(0.01,0.99),na.rm=TRUE)
#            x_thresh <- x_quants[1] + (x_quants[2] - x_quants[1]) * quantile_thresh
#            quant_date <- out_doy[which.max(x >= x_thresh)]
#            return(quant_date)
#          }))
# }
# lai_greenup <- app(lai_interp,gup_lai_fun,quantile_thresh=0.5, 
#                     filename="./output/lai_predictions_HLS_midGUP_2024.tif",
#                     overwrite=TRUE)
# 
# bdn_lai_fun <- function(x,quantile_thresh=0.5){
#   ifelse(all(is.na(x)),return(NA),
#          try({
#            x_quants <- quantile(x,probs=c(0.01,0.99),na.rm=TRUE)
#            x_thresh <- x_quants[1] + (x_quants[2] - x_quants[1]) * quantile_thresh
#            quant_date <- rev(out_doy)[which.max(rev(x) >= x_thresh)]
#            return(quant_date)
#          }))
# }
# lai_browndown <- app(lai_interp,bdn_lai_fun,quantile_thresh=0.5, 
#                     filename="./output/lai_predictions_HLS_midBDN_2024.tif",
#                     overwrite=TRUE)
