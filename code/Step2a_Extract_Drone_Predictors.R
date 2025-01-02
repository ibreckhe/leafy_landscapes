## Script to extract drone data at leafy landscapes plot locations.

## Set up workspace.
library(rstanarm)
library(sf)
library(terra)
library(readr)
library(dplyr)
#setwd("C:\\Users/IanBreckheimer/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/")
setwd("~/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/")

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
# template_paths <- list.files("./data/dronedata_templates",pattern=".tif$",full.names=TRUE)
# dem_1m <- rast("Y://data_products/released/release3/UG_dem_1m_v1.tif")
 
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
