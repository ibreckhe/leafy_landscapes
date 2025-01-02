## Script to extract drone data at leafy landscapes plot locations.

## Set up workspace.
library(rstanarm)
library(sf)
library(terra)
library(readr)
library(dplyr)
setwd("C:\\Users/IanBreckheimer/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/")
#setwd("~/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/")

## writes list list of orthos and DSMs to code by site visit
# ortho_files <- list.files("./data/ortho",full.names = TRUE)
# dsm_files <- list.files("./data/dsm/",full.names = TRUE)
# image_paths <- data.frame(ortho_path=ortho_files,
#                           dsm_path=dsm_files,
#                           visit_code=NA)
# write_csv(image_paths, "./data/drone_ortho_visit_codes.csv")

## Reads in data frame for image paths and visit codes.
image_paths <- read_csv("./data/drone_ortho_visit_codes_2023_2024.csv")
plot_polygons <- vect("./data/landscape_poly_nudge_wgs84.gpkg")
plot_polygons$visit_code <- substr(plot_polygons$Label_of_F,1,3)

## Joins polygons and image paths by visit code.
plot_polygons_join <- merge(plot_polygons,image_paths,
                            by="visit_code",all.x=TRUE)
plot_polygons_join$plot_area_sqm <- st_area(st_as_sf(plot_polygons_join))

# # Creates raster templates for each site.
# templates_ll <- image_paths$dsm_path[image_paths$template=="Y"]
# for(i in 1:length(templates_ll)){
#   outname <- gsub("/dsm/","/dronedata_templates/",templates_ll[i])
#   if(file.exists(outname)){
#     print("Template exists, skipping...")
#     next
#   }else{
#     print(paste("Projecting template",templates_ll[i]))
#     template_ll <- rast(templates_ll[i])
#     template_bin <- app(template_ll,fun=function(x){!is.na(x)})
#     template_proj <- terra::project(template_bin,y="EPSG:32613",res=0.03,method="near",
#                                     filename=outname,
#                                     overwrite=TRUE)
#   }
#  }

## Creates resampled LiDAR DEMs for each site.
 template_paths <- list.files("./data/dronedata_templates",pattern=".tif$",full.names=TRUE)
 dem_1m <- rast("Y://data_products/released/release3/UG_dem_1m_v1.tif")
 
 # for(i in 1:length(template_paths)){
 #   dem_filename <- gsub("/dronedata_templates/","/site_DEMs/",template_paths[i])
 #   if(file.exists(dem_filename)){
 #     print("DEM exists, skipping...")
 #     next()
 #   }else{
 #     print(paste("Resampling DEM for",template_paths[i]))
 #     site_template <- rast(template_paths[i])
 #     dem_crop <- crop(dem_1m,site_template,snap="out")
 #     #dem_filename <- gsub("dsm.tif$","dem.tif$",dem_filename)
 #     dem_proj <- project(dem_crop,site_template,method="bilinear",filename=dem_filename)
 #   }
 # }


## Function to crop the orthos and dsms for each plot.
crop_resample_field_plots <- function(plot_vect,ortho,dsm,template,dem,offset=15.55,
                                      buffer_width=0.96,write_output=TRUE){
    plot_cent <- centroids(plot_vect)
    plot_outer <- buffer(plot_cent,width=buffer_width * 1.5)
    ortho_crop <- crop(ortho,plot_outer)
    ortho_proj <- project(ortho_crop,template,align=TRUE,method="bilinear")
    ortho_crop2 <- crop(ortho_proj,buffer(project(plot_cent,crs(template)),width=buffer_width))
    ortho_mask <- mask(ortho_crop2,project(plot_vect,crs(template)))
    dem_crop <- crop(dem,ortho_crop2)
    ht_proj <- project(dsm,ortho_crop2,method="bilinear") + offset - dem_crop
    ht_mask <- mask(ht_proj,project(plot_vect,crs(template)))

    out <- c(ortho_crop2,ht_proj,ortho_mask,ht_mask,dem_crop)
    names(out) <- c("r_crop","g_crop","b_crop","ht_crop",
                    "r_mask","g_mask","b_mask","ht_mask",
                    "elev_DEM")
    if(write_output==TRUE){
      out_name <- paste0("./data/dronedata_crops/",plot_vect$Label_of_F,"-dronedata.tif")
      out <- writeRaster(out,out_name,overwrite=TRUE)
    }
    return(out)
}

## Tests function.
res_test <- crop_resample_field_plots(plot_vect=plot_polygons_join[20,],
                                      ortho=rast(plot_polygons_join[20,]$ortho_path),
                                      dsm=rast(plot_polygons_join[20,]$dsm_path),
                                      template=rast(plot_polygons_join[20,]$template_path),
                                      dem=rast(plot_polygons_join[20,]$dem_path),
                                      buffer_width=0.96, offset=15.55,
                                      write_output = FALSE)

## Loops through all plots and extracts raster chips for each plot.
chip_list <- list()
chip_paths <- rep(NA,nrow(plot_polygons_join))

elev <- rep(NA,nrow(plot_polygons_join))
gcc_mean <- rep(NA,nrow(plot_polygons_join))
gcc_prop <- rep(NA,nrow(plot_polygons_join))
ht_mean <- rep(NA,nrow(plot_polygons_join))
ht_max <- rep(NA,nrow(plot_polygons_join))
ht_sd <- rep(NA,nrow(plot_polygons_join))
ht_gcc_sd <- rep(NA,nrow(plot_polygons_join))
green_sd <- rep(NA,nrow(plot_polygons_join))
green_focal_sd <- rep(NA,nrow(plot_polygons_join))
bcc_mean <- rep(NA,nrow(plot_polygons_join))

plot_gcc_mean <- rep(NA,nrow(plot_polygons_join))
plot_gcc_prop <- rep(NA,nrow(plot_polygons_join))
plot_ht_mean <- rep(NA,nrow(plot_polygons_join))
plot_ht_max <- rep(NA,nrow(plot_polygons_join))
plot_ht_sd <- rep(NA,nrow(plot_polygons_join))
plot_ht_gcc_sd <- rep(NA,nrow(plot_polygons_join))
plot_green_sd <- rep(NA,nrow(plot_polygons_join))
plot_green_focal_sd <- rep(NA,nrow(plot_polygons_join))
plot_bcc_mean <- rep(NA,nrow(plot_polygons_join))

gcc_fun <- function(x){x[2]/(x[1]+x[2]+x[3])}
bcc_fun <- function(x){x[3]/(x[1]+x[2]+x[3])}

#ht_fun <- function(x){x[1] - x[2]}

for(i in 1:nrow(plot_polygons_join)){
  print(paste("Extracting data for plot",plot_polygons_join$Label_of_F[i]))
  chip_list[[i]] <- crop_resample_field_plots(plot_vect=plot_polygons_join[i,],
                            ortho=rast(plot_polygons_join[i,]$ortho_path),
                            dsm=rast(plot_polygons_join[i,]$dsm_path),
                            template=rast(plot_polygons_join[i,]$template_path),
                            dem=rast(plot_polygons_join[i,]$dem_path),
                            buffer_width=0.96,write_output = TRUE)
  names(chip_list)[i] <- plot_polygons_join$Label_of_F[i]
  
  ## Summary stats for each chip.
  chip_paths[i] <- sources(chip_list[[i]])[1]
  
  elev[i] <- global(chip_list[[i]][[9]],fun="mean",na.rm=TRUE)
  chip_gcc <- app(chip_list[[i]][[1:3]],fun=gcc_fun)
  chip_bcc <- app(chip_list[[i]][[1:3]],fun=bcc_fun)
  gcc_mean[i] <- global(chip_gcc,fun="mean",na.rm=TRUE)
  gcc_prop[i] <- global(chip_gcc > 0.39,fun="mean",na.rm=TRUE)
  ht_mean[i] <- global(chip_list[[i]][[4]],fun="mean",na.rm=TRUE)
  ht_max[i] <- global(chip_list[[i]][[4]],fun="max",na.rm=TRUE)
  ht_sd[i] <- global(chip_list[[i]][[4]],fun="sd",na.rm=TRUE)
  ht_gcc_sd[i] <- global(chip_list[[i]][[4]] * (chip_gcc > 0.39),fun="sd",na.rm=TRUE)
  green_sd[i] <- global(chip_list[[i]][[2]],fun="sd",na.rm=TRUE)
  green_focal_sd[i] <- global(focal(chip_list[[i]][[2]],w=matrix(1,3,3),fun="sd"),fun="mean",na.rm=TRUE)
  bcc_mean[i] <- global(chip_bcc,fun="mean",na.rm=TRUE)
  
  plot_gcc <- app(chip_list[[i]][[5:7]],fun=gcc_fun)
  plot_bcc <- app(chip_list[[i]][[5:7]],fun=bcc_fun)
  plot_gcc_mean[i] <- global(plot_gcc,fun="mean",na.rm=TRUE)
  plot_gcc_prop[i] <- global(plot_gcc > 0.39,fun="mean",na.rm=TRUE)
  plot_ht_mean[i] <- global(chip_list[[i]][[8]],fun="mean",na.rm=TRUE)
  plot_ht_max[i] <- global(chip_list[[i]][[8]],fun="max",na.rm=TRUE)
  plot_ht_sd[i] <- global(chip_list[[i]][[8]],fun="sd",na.rm=TRUE)
  plot_ht_gcc_sd[i] <- global(chip_list[[i]][[8]] * (plot_gcc > 0.39),fun="sd",na.rm=TRUE)
  plot_green_sd[i] <- global(chip_list[[i]][[6]],fun="sd",na.rm=TRUE)
  plot_green_focal_sd[i] <- global(focal(chip_list[[i]][[6]],w=matrix(1,3,3),fun="sd"),fun="mean",na.rm=TRUE)
  plot_bcc_mean[i] <- global(plot_bcc,fun="mean",na.rm=TRUE)
}
plot(chip_list[[150]])

## Adds summary stats to polygon data.
plot_polygons_join$elev <- elev
plot_polygons_join$crop_gcc_mean <- gcc_mean
plot_polygons_join$crop_gcc_prop <- gcc_prop
plot_polygons_join$crop_ht_mean <- ht_mean
plot_polygons_join$crop_ht_max <- ht_max
plot_polygons_join$crop_ht_sd <- ht_sd
plot_polygons_join$crop_gcc_ht_sd <- ht_gcc_sd
plot_polygons_join$crop_green_sd <- green_sd
plot_polygons_join$crop_green_focal_sd <- green_focal_sd
plot_polygons_join$crop_bcc_mean <- bcc_mean


plot_polygons_join$plot_gcc_mean <- plot_gcc_mean
plot_polygons_join$plot_gcc_prop <- plot_gcc_prop
plot_polygons_join$plot_ht_mean <- plot_ht_mean
plot_polygons_join$plot_ht_max <- plot_ht_max
plot_polygons_join$plot_ht_sd <- plot_ht_sd
plot_polygons_join$plot_gcc_ht_sd <- plot_ht_gcc_sd
plot_polygons_join$plot_green_sd <- plot_green_sd
plot_polygons_join$plot_green_focal_sd <- plot_green_focal_sd
plot_polygons_join$plot_bcc_mean <- plot_bcc_mean


chip_paths_sub <- gsub("/Users/ian/Library/CloudStorage/OneDrive-RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/",
                      ".",chip_paths,fixed=TRUE) 

plot_polygons_join$chip_path <- chip_paths_sub

## Quick visualizations.
plot(plot_polygons_join$plot_ht_sd,plot_polygons_join$plot_gcc_prop,cex=0)
text(plot_polygons_join$plot_ht_sd,plot_polygons_join$plot_gcc_prop,
     labels=plot_polygons_join$Label_of_F,cex=0.8,col=as.numeric(factor(plot_polygons_join$visit_code)))

plot(plot_polygons_join$crop_ht_sd,plot_polygons_join$plot_gcc_mean,cex=0)
text(plot_polygons_join$crop_ht_sd,plot_polygons_join$plot_gcc_mean,
     labels=plot_polygons_join$Label_of_F,cex=0.8,col=as.numeric(factor(plot_polygons_join$visit_code)))

plot(plot_polygons_join$crop_gcc_mean,plot_polygons_join$plot_gcc_mean,cex=0)
text(plot_polygons_join$crop_gcc_mean,plot_polygons_join$plot_gcc_mean,
     labels=plot_polygons_join$Label_of_F,cex=0.8,col=as.numeric(factor(plot_polygons_join$visit_code)))

plot(plot_polygons_join$crop_green_sd,plot_polygons_join$plot_green_sd,cex=0)
text(plot_polygons_join$crop_green_sd,plot_polygons_join$plot_green_sd,
     labels=plot_polygons_join$Label_of_F,cex=0.8,col=as.numeric(factor(plot_polygons_join$visit_code)))

plot(plot_polygons_join$plot_green_sd,plot_polygons_join$plot_green_focal_sd,cex=0)
text(plot_polygons_join$plot_green_sd,plot_polygons_join$plot_green_focal_sd,
     labels=plot_polygons_join$Label_of_F,cex=0.8,col=as.numeric(factor(plot_polygons_join$visit_code)))


## Joins with laboratory data.
plot_lab_data <- read_csv("./data/finalized_plot_data/leafy_landscapes_plot_data_v1.csv")
plot_lab_data_gcc <- left_join(plot_lab_data,as.data.frame(plot_polygons_join),by=c("PlotID"="Label_of_F"))

## Calculates slope-adjusted LAI with plot areas calculated from polygons.
plot_lab_data_gcc$Plot_adjusted_LAI <- plot_lab_data_gcc$Plot_total_leaf_area / as.numeric(plot_lab_data_gcc$plot_area_sqm)
plot_lab_data_gcc$Plot_adjusted_LAI[plot_lab_data_gcc$Plot_adjusted_LAI == 0] <- 0.1
plot_lab_data_gcc$plot_gcc_mean_sq <- plot_lab_data_gcc$plot_gcc_mean^2
plot_lab_data_gcc$Plot_leaf_water_pct <- plot_lab_data_gcc$Plot_total_canopy_water_content / (plot_lab_data_gcc$Plot_total_canopy_water_content + plot_lab_data_gcc$Plot_total_dry_mass)
plot_lab_data_gcc$Plot_leaf_water_pct_logit <- logit(plot_lab_data_gcc$Plot_leaf_water_pct)


write_csv(plot_lab_data_gcc,file="./data/finalized_plot_data/leafy_landscapes_plot_drone_data_v1.csv")
plot_lab_data_gcc <- read_csv("./data/finalized_plot_data/leafy_landscapes_plot_drone_data_v1.csv")

## Plots of covariation between LAI, LMA, and Canopy Water Content.
library(ggplot2)
plot_lab_data_gcc$Collected_DOY <- as.numeric(format(as.Date(plot_lab_data_gcc$Collected_date), "%j"))
plot_lab_data_gcc$Collected_early <- plot_lab_data_gcc$Collected_DOY < 185
plot_lab_data_gcc$Site_visit_num <- NA
plot_lab_data_gcc$Site_visit_num[grepl("1", plot_lab_data_gcc$visit_code)] <- "1"
plot_lab_data_gcc$Site_visit_num[grepl("2", plot_lab_data_gcc$visit_code)] <- "2"
plot_lab_data_gcc$Site_visit_num[grepl("3", plot_lab_data_gcc$visit_code)] <- "3"
plot_lab_data_gcc$Site_visit_num[grepl("4", plot_lab_data_gcc$visit_code)] <- "4"


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


## Computes site-level masks for the LAI predictions.
lai_preds <- list.files("./output/lai_predictions/",pattern="mod5.tif$",full.names=TRUE)

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

