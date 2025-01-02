## Script to apply the HLS LAI model to Gunnison Basin HLS data from 2015-2024.
## Author: Ian Breckheimer
## Updated: 1 January 2025

## This script should be run after we are done fitting models relating HLS data to upscaled observations. To see how this script fits into the larger analysis, see the documentation here: https://app.affine.pro/workspace/e094b233-4dce-4599-bc85-a12ab329bda3/aKjjrjpiHUO9Bp0A9mblh?mode=page

# Load libraries
library(dplyr)
library(terra)
library(readr)
library(tidyr)
library(stringr)
library(earthdatalogin)

## Set working directory.
wd <- "/Users/ian/Library/CloudStorage/OneDrive-RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI"
#wd <- "C://Users/IanBreckheimer/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI"
setwd(wd)

## Parameters for data access via NASA Earth Data login.
edl_netrc()

## Training data for LAI model
drone_df_join <- read.csv("./output/upscaled_lai_q50_predictions_30m_HLS_interp.csv")

## Reads in fit model.
gam1 <- readRDS("./output/HLS_LAI_model_gam1.rds")

## Reads in study extent.
#study_extent <- ext(rast("./data/HLS/HLS_mask_UER.tif"))
study_extent <- ext(305073,388098,4256064,4328667)

## Gets the granule list from AppEEARS.
granule_list <- c(as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2015-granule-list.txt",header=FALSE)[,1]),
                  as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2016-granule-list.txt",header=FALSE)[,1]),
                  as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2017-granule-list.txt",header=FALSE)[,1]),
                  as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2018-granule-list.txt",header=FALSE)[,1]),
                  as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2019-granule-list.txt",header=FALSE)[,1]),
                  as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2020-granule-list.txt",header=FALSE)[,1]),
                  as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2021-granule-list.txt",header=FALSE)[,1]),
                  as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2022-granule-list.txt",header=FALSE)[,1]),
                  as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2023-granule-list.txt",header=FALSE)[,1]),
                  as.vector(read.csv("./data/HLS/test_granules/HLS-L30-2024-granule-list.txt",header=FALSE)[,1]))

products <- str_split_fixed(as.character(granule_list),pattern=fixed("."),n=20)[,7]
dates <- str_split_fixed(as.character(granule_list),pattern=fixed("."),n=20)[,14]
bands <- str_split_fixed(as.character(granule_list),pattern=fixed("."),n=20)[,17]
tiles <- str_split_fixed(as.character(granule_list),pattern=fixed("."),n=20)[,13]
imagery_df <- data.frame(url=granule_list,
                         product=products,
                         band=bands,
                         date=dates,
                         tiles=tiles)

study_tiles <- c("T13SCD", "T13SCC")
fmask_codes <- c(64, 68, 80, 84, 96, 100, 112, 116, 128, 
                 132, 144, 148, 160, 164, 176, 180,
                 192, 196, 208, 212, 224, 228, 240, 244)
imagery_df <- filter(imagery_df, tiles %in% study_tiles, band %in% c("B03", "B04", "B05", "Fmask"))
imagery_df$date_scene <- paste0(imagery_df$date, "-", imagery_df$tiles)
imagery_dates_scenes <- unique(imagery_df$date_scene)
imagery_dates <- stringr::str_split_fixed(imagery_dates_scenes,pattern="-",n=2)[,1]
imagery_scenes <- stringr::str_split_fixed(imagery_dates_scenes,pattern="-",n=2)[,2]

##
#test_urls <- imagery_df$url[imagery_df$date == unique(imagery_df$date)[1]]
#test_rast <- crop(rast(test_urls), study_extent)

## Function for calculating LAI predictors.
ndgi_fun <- function(x){
  ndgi <- ifelse(x[4] %in% c(64, 68, 80, 84, 96, 100, 112, 116, 128, 
                             132, 144, 148, 160, 164, 176, 180,
                             192, 196, 208, 212, 224, 228, 240, 244), 
                 (0.65 * x[1] + (1 - 0.65) * x[3] - x[2]) / 
                   (0.65 * x[1] + (1 - 0.65) * x[3] + x[2] + 0.001), 
                 NA)
  ndgi[ndgi < -1] <- -1
  ndgi[ndgi > 1.2] <- 1.2
  return(ndgi)
}

nirv_fun <- function(x){
  nirv <- ifelse(x[4] %in% c(64, 68, 80, 84, 96, 100, 112, 116, 128, 
                             132, 144, 148, 160, 164, 176, 180,
                             192, 196, 208, 212, 224, 228, 240, 244), 
                 ((x[3]-x[2]) / (x[3]+x[2] + 0.001)) * (((x[3]-x[2]) / (x[3]+x[2]) + 0.001) > 0) * x[3] , 
                 NA)
  nirv[nirv < 0] <- 0
  nirv[nirv > 1.1] <- 1.1
  return(nirv)
}

## Testing function
#test_ndgi <- app(test_rast,fun=ndgi_fun)
#names(test_ndgi) <- "NDGI"
#test_lai <- predict(test_ndgi, model=gam1$gam, type="response")
#plot(test_lai)

## Loops through all scenes and predicts LAI.
calc_veg_index_fun <- function(imagery_df,pred_date_scene,pred_extent,overwrite=FALSE){
  setwd(wd)
  edl_netrc()
  urls <- imagery_df$url[imagery_df$date_scene == pred_date_scene]
  prod <- imagery_df$product[imagery_df$date_scene == pred_date_scene][1]
  tile <- imagery_df$tiles[imagery_df$date_scene == pred_date_scene][1]
  pred_date <- stringr::str_split_fixed(pred_date_scene,pattern="-",n=2)[,1]
  
  print(paste0("Predicting Veg indices for ", prod,"at tile",tile," on ", pred_date))
  
  outname_ndgi <- paste0("./scratch/HLS_NDGI_scenes/HLS_", 
                         tile, "_",prod,"_", pred_date, "_NDGI.tif")
  outname_nirv <- gsub("NDGI","NIRv",outname_ndgi)
  
  if(file.exists(outname_ndgi) & file.exists(outname_nirv) & overwrite==FALSE){
    print(paste("File",outname_ndgi,"or",outname_nirv,"already exists. Skipping..."))
    return(NULL)
  }
  
  print("Preparing Inputs...")
  input_rast <- terra::crop(rast(urls), pred_extent)
  
  print("Calculating Outputs...")
  out_ndgi <- terra::app(input_rast,fun=ndgi_fun,cores=1)
  names(out_ndgi) <- "NDGI"
  out_nirv <- terra::app(input_rast,fun=nirv_fun,cores=1)
  names(out_nirv) <- "NIRv"
  
  #print("Predicting LAI...")
  #out_lai <- terra::round(terra::predict(input_pred, gam1$gam, type="response") * 1000)
  #out_lai[out_lai < 0] <- 0
  #out_lai[out_lai > 10000] <- 10000
  #crs(out_lai) <- "EPSG:32613"
  
  terra::writeRaster(out_ndgi, 
                     filename=outname_ndgi,
                     filetype="COG",
                     overwrite=TRUE)
  
  terra::writeRaster(out_nirv, 
                     filename=outname_nirv,
                     filetype="COG",
                     overwrite=TRUE)
  #return(out_ndgi)
}

library(foreach)
library(doParallel)

## Set up parallel processing.
#cl <- makeCluster(8)
#registerDoParallel(cl)

## HLS L30 2016.
out_test_lai <- foreach(i=1:length(imagery_dates), .combine=list, 
                        .packages=c("terra","mgcv","earthdatalogin"),
                        .errorhandling = "pass") %do% {
                          print(paste("Processing scene", i, "of", length(imagery_dates), ": " ,imagery_dates_scenes[i]))
                          edl_netrc()
                          calc_veg_index_fun(imagery_df,imagery_dates_scenes[i],pred_extent=ext(305073,388098,4256064,4328667),
                                             overwrite=FALSE)
                        }
#stopCluster(cl)

s2_granule_list <- c(as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2015-granule-list.txt",header=FALSE)[,1]),
                     as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2016-granule-list.txt",header=FALSE)[,1]),
                     as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2017-granule-list.txt",header=FALSE)[,1]),
                     as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2018-granule-list.txt",header=FALSE)[,1]),
                     as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2019-granule-list.txt",header=FALSE)[,1]),
                     as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2020-granule-list.txt",header=FALSE)[,1]),
                     as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2021-granule-list.txt",header=FALSE)[,1]),
                     as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2022-granule-list.txt",header=FALSE)[,1]),
                     as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2023-granule-list.txt",header=FALSE)[,1]),
                     as.vector(read.csv("./data/HLS/test_granules/HLS-S30-2024-granule-list.txt",header=FALSE)[,1]))


s2_products <- str_split_fixed(as.character(s2_granule_list),pattern=fixed("."),n=20)[,7]
s2_dates <- str_split_fixed(as.character(s2_granule_list),pattern=fixed("."),n=20)[,14]
s2_bands <- str_split_fixed(as.character(s2_granule_list),pattern=fixed("."),n=20)[,17]
s2_tiles <- str_split_fixed(as.character(s2_granule_list),pattern=fixed("."),n=20)[,13]
s2_imagery_df <- data.frame(url=s2_granule_list,
                            product=s2_products,
                            band=s2_bands,
                            date=s2_dates,
                            tiles=s2_tiles)

study_tile <- c("T13SCD","T13SCC")
fmask_codes <- c(64, 68, 80, 84, 96, 100, 112, 116, 128, 
                 132, 144, 148, 160, 164, 176, 180,
                 192, 196, 208, 212, 224, 228, 240, 244)
s2_imagery_df <- filter(s2_imagery_df, tiles %in% study_tiles, band %in% c("B03", "B04", "B08", "Fmask"))
s2_imagery_df$date_scene <- paste0(s2_imagery_df$date, "-", s2_imagery_df$tiles)
s2_imagery_dates_scenes <- unique(s2_imagery_df$date_scene)
s2_imagery_dates <- stringr::str_split_fixed(s2_imagery_dates_scenes,pattern="-",n=2)[,1]
s2_imagery_scenes <- stringr::str_split_fixed(s2_imagery_dates_scenes,pattern="-",n=2)[,2]

## Set up parallel processing.
cl <- makeCluster(8)
registerDoParallel(cl)

## HLS L30.
out_test_lai <- foreach(i=1:length(s2_imagery_dates), .combine=list, 
                        .packages=c("terra","mgcv","earthdatalogin"),
                        .errorhandling = "pass") %do% {
                          print(paste("Processing scene", i, "of", length(s2_imagery_dates), ": " ,s2_imagery_dates_scenes[i]))
                          edl_netrc()
                          calc_veg_index_fun(s2_imagery_df,s2_imagery_dates_scenes[i],pred_extent=ext(305073,388098,4256064,4328667),
                                             overwrite=FALSE)
                        }
stopCluster(cl)


### Predicts LAI from the veg index data.
pred_lai_fun <- function(image_path,lai_model,tile,pred_date,pred_extent,overwrite=FALSE){
  setwd(wd)
  
  pred_doy <- format(pred_date,format="%j")
  outname_lai <- gsub("NDGI","LAI",image_path)
  outname_lai <- gsub("_scenes","_scenes_UG",outname_lai)
  
  if(file.exists(outname_lai) & overwrite==FALSE){
    print(paste("File",outname_lai, "already exists. Skipping..."))
    return(NULL)
  }
  
  print("Preparing Inputs...")
  input_rast <- terra::crop(rast(image_path), pred_extent)
  
  print("Predicting LAI...")
  out_lai <- terra::round(terra::predict(input_rast, lai_model$gam, type="response") * 1000)
  out_lai[out_lai < 0] <- 0
  out_lai[out_lai > 10000] <- 10000
  crs(out_lai) <- "EPSG:32613"
  
  terra::writeRaster(out_lai, 
                     filename=outname_lai,
                     datatype="INT2U",
                     overwrite=TRUE)
  #return(out_ndgi)
}

## Reads in paths of all LAI predictions and parses dates from filenames
lai_df_all <- data.frame(path=list.files(path="./scratch/HLS_NDGI_scenes/",pattern=".tif$",full.names=TRUE))
lai_df_all$Date <- as.Date(substr(stringr::str_split_fixed(lai_df_all$path,"_|\\.",9)[,7],1,7),format="%Y%j")
lai_df_all$Tile <- substr(stringr::str_split_fixed(lai_df_all$path,"_|\\.",9)[,5],1,6)
lai_df_all$Year <- format(lai_df_all$Date,format="%Y")
lai_df_all <- lai_df_all[order(lai_df_all$Date),]

lai_df_year <- filter(lai_df_all, Year=="2023", Tile %in% c("T13SCC","T13SCD"))
#lai_pred_test <- pred_lai_fun(image_path=lai_df_year_tile$path[70], 
#                              tile=lai_df_year_tile$Tile[70], 
#                              pred_date=lai_df_year_tile$Date[70],
#                              pred_extent=ext(305073,388098,4256064,4328667),
#                              lai_model=gam1, overwrite=TRUE)

## Predicts LAI for all scenes.
library(foreach)
library(doParallel)

## Set up parallel processing.
cl <- makeCluster(8)
registerDoParallel(cl)

foreach(i=1:nrow(lai_df_all), .combine=list, 
        .packages=c("terra","mgcv","earthdatalogin"),
        .errorhandling = "pass") %dopar% {
          outlist <- pred_lai_fun(image_path=lai_df_all$path[i], 
                                  tile=lai_df_all$Tile[i], 
                                  pred_date=lai_df_all$Date[i],
                                  pred_extent=ext(305073,388098,4256064,4328667),
                                  lai_model=gam1, overwrite=FALSE)
        }

stopCluster(cl)
