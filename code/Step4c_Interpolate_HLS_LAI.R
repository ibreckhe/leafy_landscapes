## Script to interpolate and gap-fill HLS LAI predictions before computing seasonal summaries.
## Author: Ian Breckheimer
## Updated: 1 January 2025

## This script should be run after generating LAI predictions for each HLS satellite scene. To see how this script fits into the larger analysis, see the documentation here: https://app.affine.pro/workspace/e094b233-4dce-4599-bc85-a12ab329bda3/aKjjrjpiHUO9Bp0A9mblh?mode=page

## Load libraries
library(dplyr)
library(terra)
library(readr)
library(tidyr)
library(stringr)

## Set working directory.
wd <- "/Users/ian/Library/CloudStorage/OneDrive-RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI"
#wd <- "C://Users/IanBreckheimer/OneDrive - RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI"
setwd(wd)

## Reads in paths of all LAI predictions and parses dates from filenames
lai_df_all <- data.frame(path=list.files(path="./scratch/HLS_NDGI_scenes/",pattern=".tif$",full.names=TRUE))
lai_df_all$Date <- as.Date(substr(stringr::str_split_fixed(lai_df_all$path,"_|\\.",9)[,7],1,7),format="%Y%j")
lai_df_all$Tile <- substr(stringr::str_split_fixed(lai_df_all$path,"_|\\.",9)[,5],1,6)
lai_df_all$Year <- format(lai_df_all$Date,format="%Y")
lai_df_all <- lai_df_all[order(lai_df_all$Date),]

lai_df_year <- filter(lai_df_all, Year=="2023", Tile %in% c("T13SCC","T13SCD"))

#### Interpolates LAI for each year.
loess_fun <- function(x,dates=as.numeric(lai_df_year$Date)){
  ifelse(all(is.na(x)),return(rep(NA,length(out_dates))),
         try({
           loess_mod <- loess(x ~ dates, span=0.75,
                              control = loess.control(surface = "direct"))
           loess_pred <- predict(loess_mod, newdata=as.numeric(out_dates))
           loess_pred[loess_pred < 0] <- 0
           loess_pred[loess_pred > 10000] <- 10000
           return(round(loess_pred))
         }))
  
}
years <- c(2015:2024)

for(i in 1:length(years)){
  
  print(paste("Interpolating and summarizing for tile T13SCD in year",years[i]))
  lai_df_year <- filter(lai_df_all,Year==years[i], Tile=="T13SCD")
  
  lai_df_year$outpath <- gsub("NDGI","LAI",lai_df_year$path)
  lai_df_year$outpath <- gsub("_scenes","_scenes_UG",lai_df_year$outpath)
  
  
  out_dates <- seq(as.Date(paste0(years[i],"-90"),format="%Y-%j"),
                   as.Date(paste0(years[i],"-305"),format="%Y-%j"), by=5)
  
  interp_filename <- paste0("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_",years[i],"_T13SCD.tif")
  
  #if(!file.exists(interp_filename)){
  lai_rast_year <- rast(lai_df_year$outpath)
  
  #plot(lai_rast_year[[25]])
  #test_ext <- ext(325089, 331939, 4310449, 4315310)
  #lai_rast_year_test <- crop(lai_rast_year,test_ext)
  
  lai_interp <- app(lai_rast_year,fun=loess_fun, filename=interp_filename,
                    overwrite=TRUE, wopt=list(datatype="INT2U",
                                              names=paste0("LAI-",format(out_dates,"%Y-doy%j"))))
  names(lai_interp) <- out_dates
  #}else{
  #  print(paste("Loading interpolation from existing files..."))
  #  lai_interp <- rast(interp_filename)
  #}
  
  plot(lai_interp[[seq(1,nlyr(lai_interp),by=floor(nlyr(lai_interp)/3))]],range=c(0,6000),main=years[i])
  #lai_mem <- lai_rast + 0
  
  pred_interval <- as.numeric(out_dates[2] - out_dates[1])
  integral_fun <- function(x,interval=pred_interval){
    ifelse(all(is.na(x)),return(NA),
           try({
             integral <- sum((x/1000) * interval)
             return(integral)
           }))
  }
  integral_outname <- paste0("./output/lai_predictions_HLS_integral_",years[i],"_T13SCD.tif")
  lai_integral <- app(lai_interp,integral_fun,filename=integral_outname,
                      overwrite=TRUE)
  
  ## Computes date of max LAI.
  out_doy <- as.numeric(format(out_dates,"%j"))
  peak_lai_fun <- function(x){
    ifelse(all(is.na(x)),return(NA),
           try({
             max_date <- out_doy[which.max(x)]
             return(max_date)
           }))
  }
  lai_peak <- app(lai_interp, peak_lai_fun,filename=gsub("integral","peakDOY",integral_outname),
                  overwrite=TRUE)
  
  ## Computes the date halfway through the cumulative LAI distribution.
  half_lai_fun <- function(x){
    period <- out_doy[2] - out_doy[1]
    ifelse(all(is.na(x)),return(NA),
           try({
             lai_cum <- cumsum(x*period)
             lai_half <- max(lai_cum) / 2
             half_date <- out_doy[which.max(lai_cum > lai_half)] - period/2
             return(half_date)
           }))
  }
  lai_half <- app(lai_interp,half_lai_fun,filename=gsub("integral","halfDOY",integral_outname),
                  overwrite=TRUE)
  
  ## Computes the date 1/4 of the way through the cumulative LAI distribution.
  oneforth_lai_fun <- function(x){
    period <- out_doy[2] - out_doy[1]
    ifelse(all(is.na(x)),return(NA),
           try({
             lai_cum <- cumsum(x*period)
             lai_quarter <- max(lai_cum) / 4
             quart_date <- out_doy[which.max(lai_cum > lai_quarter)] - period/2
             return(quart_date)
           }))
  }
  lai_oneforth <- app(lai_interp,oneforth_lai_fun,filename=gsub("integral","oneforthDOY",integral_outname),
                      overwrite=TRUE)
  
  threeforth_lai_fun <- function(x){
    period <- out_doy[2] - out_doy[1]
    ifelse(all(is.na(x)),return(NA),
           try({
             lai_cum <- cumsum(x*period)
             lai_threequart <- max(lai_cum) * 0.75
             threequart_date <- out_doy[which.max(lai_cum > lai_threequart)] - period/2
             return(threequart_date)
           }))
  }
  lai_threeforth<- app(lai_interp,threeforth_lai_fun,filename=gsub("integral","threeforthDOY",integral_outname),
                       overwrite=TRUE)
  
  amp_lai_fun <- function(x){
    ifelse(all(is.na(x)),return(NA),
           try({
             amp <- quantile(x,0.95,na.rm=TRUE) - quantile(x,0.05,na.rm=TRUE)
             return(amp)
           }))
  }
  lai_amp <- app(lai_interp,amp_lai_fun,filename=gsub("integral","LAIamp",integral_outname),
                 overwrite=TRUE)
  
  
  ## Computes dates reaching quantiles of max LAI.
  gup_lai_fun <- function(x,quantile_thresh=0.25){
    ifelse(all(is.na(x)),return(NA),
           try({
             x_quants <- quantile(x,probs=c(0.01,0.99),na.rm=TRUE)
             x_thresh <- x_quants[1] + (x_quants[2] - x_quants[1]) * quantile_thresh
             quant_date <- out_doy[which.max(x >= x_thresh)]
             return(quant_date)
           }))
  }
  lai_greenup <- app(lai_interp,gup_lai_fun,quantile_thresh=0.5, 
                     filename=gsub("integral","midGUP",integral_outname),
                     overwrite=TRUE)
  
  bdn_lai_fun <- function(x,quantile_thresh=0.5){
    ifelse(all(is.na(x)),return(NA),
           try({
             x_quants <- quantile(x,probs=c(0.01,0.99),na.rm=TRUE)
             x_thresh <- x_quants[1] + (x_quants[2] - x_quants[1]) * quantile_thresh
             quant_date <- rev(out_doy)[which.max(rev(x) >= x_thresh)]
             return(quant_date)
           }))
  }
  lai_browndown <- app(lai_interp,bdn_lai_fun,quantile_thresh=0.5, 
                       filename=gsub("integral","midBDN",integral_outname),
                       overwrite=TRUE)
}

years <- c(2015:2024)

for(i in 1:length(years)){
  
  print(paste("Interpolating and summarizing for tile T13SCC in year",years[i]))
  lai_df_year <- filter(lai_df_all,Year==years[i], Tile=="T13SCC")
  
  lai_df_year$outpath <- gsub("NDGI","LAI",lai_df_year$path)
  lai_df_year$outpath <- gsub("_scenes","_scenes_UG",lai_df_year$outpath)
  
  
  out_dates <- seq(as.Date(paste0(years[i],"-90"),format="%Y-%j"),
                   as.Date(paste0(years[i],"-305"),format="%Y-%j"), by=5)
  
  interp_filename <- paste0("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_",years[i],"_T13SCC.tif")
  
  #if(!file.exists(interp_filename)){
  lai_rast_year <- rast(lai_df_year$outpath)
  
  #plot(lai_rast_year[[32]])
  #test_ext <- draw()
  #test_ext <- ext(333707, 343430, 4281911, 4292076)
  #lai_rast_year_test <- crop(lai_rast_year,test_ext)
  
  lai_interp <- app(lai_rast_year,fun=loess_fun, filename=interp_filename,
                    overwrite=TRUE, wopt=list(datatype="INT2U",
                                              names=paste0("LAI-",format(out_dates,"%Y-doy%j"))))
  names(lai_interp) <- out_dates
  #}else{
  #  print(paste("Loading interpolation from existing files..."))
  #  lai_interp <- rast(interp_filename)
  #}
  
  plot(lai_interp[[seq(1,nlyr(lai_interp),by=floor(nlyr(lai_interp)/3))]],range=c(0,6000),main=years[i])
  #lai_mem <- lai_rast + 0
  
  pred_interval <- as.numeric(out_dates[2] - out_dates[1])
  integral_fun <- function(x,interval=pred_interval){
    ifelse(all(is.na(x)),return(NA),
           try({
             integral <- sum((x/1000) * interval)
             return(integral)
           }))
  }
  integral_outname <- paste0("./output/lai_predictions_HLS_integral_",years[i],"_T13SCC.tif")
  lai_integral <- app(lai_interp,integral_fun,filename=integral_outname,
                      overwrite=TRUE)
  
  ## Computes date of max LAI.
  out_doy <- as.numeric(format(out_dates,"%j"))
  peak_lai_fun <- function(x){
    ifelse(all(is.na(x)),return(NA),
           try({
             max_date <- out_doy[which.max(x)]
             return(max_date)
           }))
  }
  lai_peak <- app(lai_interp, peak_lai_fun,filename=gsub("integral","peakDOY",integral_outname),
                  overwrite=TRUE)
  
  ## Computes the date halfway through the cumulative LAI distribution.
  half_lai_fun <- function(x){
    period <- out_doy[2] - out_doy[1]
    ifelse(all(is.na(x)),return(NA),
           try({
             lai_cum <- cumsum(x*period)
             lai_half <- max(lai_cum) / 2
             half_date <- out_doy[which.max(lai_cum > lai_half)] - period/2
             return(half_date)
           }))
  }
  lai_half <- app(lai_interp,half_lai_fun,filename=gsub("integral","halfDOY",integral_outname),
                  overwrite=TRUE)
  
  ## Computes the date 1/4 of the way through the cumulative LAI distribution.
  oneforth_lai_fun <- function(x){
    period <- out_doy[2] - out_doy[1]
    ifelse(all(is.na(x)),return(NA),
           try({
             lai_cum <- cumsum(x*period)
             lai_quarter <- max(lai_cum) / 4
             quart_date <- out_doy[which.max(lai_cum > lai_quarter)] - period/2
             return(quart_date)
           }))
  }
  lai_oneforth <- app(lai_interp,oneforth_lai_fun,filename=gsub("integral","oneforthDOY",integral_outname),
                      overwrite=TRUE)
  
  threeforth_lai_fun <- function(x){
    period <- out_doy[2] - out_doy[1]
    ifelse(all(is.na(x)),return(NA),
           try({
             lai_cum <- cumsum(x*period)
             lai_threequart <- max(lai_cum) * 0.75
             threequart_date <- out_doy[which.max(lai_cum > lai_threequart)] - period/2
             return(threequart_date)
           }))
  }
  lai_threeforth<- app(lai_interp,threeforth_lai_fun,filename=gsub("integral","threeforthDOY",integral_outname),
                       overwrite=TRUE)
  
  amp_lai_fun <- function(x){
    ifelse(all(is.na(x)),return(NA),
           try({
             amp <- quantile(x,0.95,na.rm=TRUE) - quantile(x,0.05,na.rm=TRUE)
             return(amp)
           }))
  }
  lai_amp <- app(lai_interp,amp_lai_fun,filename=gsub("integral","LAIamp",integral_outname),
                 overwrite=TRUE)
  
  
  ## Computes dates reaching quantiles of max LAI.
  gup_lai_fun <- function(x,quantile_thresh=0.25){
    ifelse(all(is.na(x)),return(NA),
           try({
             x_quants <- quantile(x,probs=c(0.01,0.99),na.rm=TRUE)
             x_thresh <- x_quants[1] + (x_quants[2] - x_quants[1]) * quantile_thresh
             quant_date <- out_doy[which.max(x >= x_thresh)]
             return(quant_date)
           }))
  }
  lai_greenup <- app(lai_interp,gup_lai_fun,quantile_thresh=0.5, 
                     filename=gsub("integral","midGUP",integral_outname),
                     overwrite=TRUE)
  
  bdn_lai_fun <- function(x,quantile_thresh=0.5){
    ifelse(all(is.na(x)),return(NA),
           try({
             x_quants <- quantile(x,probs=c(0.01,0.99),na.rm=TRUE)
             x_thresh <- x_quants[1] + (x_quants[2] - x_quants[1]) * quantile_thresh
             quant_date <- rev(out_doy)[which.max(rev(x) >= x_thresh)]
             return(quant_date)
           }))
  }
  lai_browndown <- app(lai_interp,bdn_lai_fun,quantile_thresh=0.5, 
                       filename=gsub("integral","midBDN",integral_outname),
                       overwrite=TRUE)
}

## Joins data from the different tiles.
interps_t13scc <- list.files("./scratch/HLS_LAI_interps/",pattern="T13SCC.tif$",
                             full.names=TRUE)
rasts_t13scc <- lapply(interps_t13scc,rast)
interps_t13scd <- list.files("./scratch/HLS_LAI_interps/",pattern="T13SCD.tif$",
                             full.names=TRUE)
rasts_t13scd <- lapply(interps_t13scd,rast)

interp_years <- 2015:2024
interp_mosaics <- list()

for(i in 1:length(interp_years)){
  print(paste("Merging files for year",interp_years[i]))
  interp_mosaics[[i]] <- merge(rasts_t13scc[[i]],rasts_t13scd[[i]],
                               filename=paste0("./scratch/HLS_LAI_interps/UG_LAI_Interp_",interp_years[i],".tif"),
                               overwrite=TRUE)
  
}

t13scc_files <- list.files("./output/",pattern="T13SCC.tif$",full.names=TRUE)
t13scd_files <- list.files("./output/",pattern="T13SCD.tif$",full.names=TRUE)


prods <- c("midGUP","midBDN","LAIamp","halfDOY","threeforthDOY","oneforthDOY","peakDOY","integral")
prod_merges <- list()
prod_means <- list()
prod_medians <- list()
prod_anoms <- list()

for(i in 1:length(prods)){
  print(paste("Merging tiles for product",prods[i]))
  t13scc_prod <- rast(t13scc_files[grepl(prods[i],t13scc_files)])
  t13scd_prod <- rast(t13scd_files[grepl(prods[i],t13scd_files)])
  prod_filename <- paste0("./output/UG_",prods[i],"_HLS_2015_2024.tif")
  
  prod_merge <- merge(t13scc_prod,t13scd_prod)
  names(prod_merge) <- paste0(prods[i],"_",interp_years)
  #writeRaster(prod_merge,filename=prod_filename,overwrite=TRUE)
  #prod_merges[[i]] <- prod_merge
  
  #prod_mean <- mean(prod_merge)
  #names(prod_mean) <- paste0(prods[i],"_2015_2024_mean")
  #writeRaster(prod_mean,filename=gsub("2024","2024_mean",prod_filename),overwrite=TRUE)
  #prod_means[[i]] <- prod_mean
  
  prod_median <- median(prod_merge)
  names(prod_median) <- paste0(prods[i],"_2015_2024_median")
  writeRaster(prod_median,filename=gsub("2024","2024_median",prod_filename),overwrite=TRUE)
  prod_means[[i]] <- prod_mean
  
  #prod_anom <- prod_merge - prod_mean
  #names(prod_anom) <- paste0(prods[i],"_anom_",interp_years) 
  #writeRaster(prod_anom,filename=gsub("HLS","HLS_anom",prod_filename),overwrite=TRUE)
  #prod_anoms[[i]] <- prod_anom
}

## Visualizes interpolated LAI.
library(ggplot2)
library(dplyr)
library(tidyterra)
library(terra)

plot_years <- 2015:2024

for(i in 1:length(plot_years)){
  print(paste("Creating plot for year",plot_years[i]))
  lai_interp_year <- rast(paste0("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_",
                                 plot_years[i],".tif"))
  
  lai_plot <- ggplot()+
    geom_spatraster(data=lai_interp_year[[seq(16,nlyr(lai_interp_year),by=3)]])+
    scale_fill_whitebox_c(limits=c(0,5000),palette="viridi")+
    facet_wrap(facets=~lyr)+
    theme_minimal()
  
  plot_path <- paste0("./figs/HLS_lai_interp_",plot_years[i],".pdf")
  
  pdf(plot_path,width=10,height=10)
  print(lai_plot)
  dev.off()
  
}

## Gets lists of files for each attribute.
integral_files <- list.files("./output/",pattern="integral",full.names=TRUE)
integral_rast <- rast(integral_files)
integral_mean <- mean(integral_rast)
integral_anom <- integral_rast - integral_mean
writeRaster(integral_anom,filename="./output/lai_int_anom_2015_2024.tif",overwrite=TRUE)

oneforth_files <- list.files("./output/",pattern="oneforth",full.names=TRUE)
oneforth_rast <- rast(oneforth_files)
oneforth_mean <- mean(oneforth_rast)
oneforth_anom <- oneforth_rast - oneforth_mean
writeRaster(oneforth_anom,filename="./output/lai_quarter_anom_2015_2024.tif",overwrite=TRUE)

half_files <- list.files("./output/",pattern="halfDOY",full.names=TRUE)
half_rast <- rast(half_files)
half_mean <- mean(half_rast)
half_anom <- half_rast - half_mean
writeRaster(half_anom,filename="./output/lai_half_anom_2015_2024.tif",overwrite=TRUE)

threeq_files <- list.files("./output/",pattern="threeforthDOY",full.names=TRUE)
threeq_rast <- rast(threeq_files)
threeq_mean <- mean(threeq_rast)
threeq_anom <- threeq_rast - threeq_mean
writeRaster(threeq_anom,filename="./output/lai_threeq_anom_2015_2024.tif",overwrite=TRUE)

## Calculates annual medians and anomalies.
interp_year_2015 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2015.tif")
interp_year_2016 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2016.tif")
interp_year_2017 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2017.tif")
interp_year_2018 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2018.tif")
interp_year_2019 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2019.tif")
interp_year_2020 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2020.tif")
interp_year_2021 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2021.tif")
interp_year_2022 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2022.tif")
interp_year_2023 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2023.tif")
interp_year_2024 <- rast("./scratch/HLS_LAI_interps/lai_predictions_HLS_interp_mod5_2024.tif")

interp_mean <- mean(interp_year_2015,
                    interp_year_2016,
                    interp_year_2017,
                    interp_year_2018,
                    interp_year_2019,
                    interp_year_2020,
                    interp_year_2021,
                    interp_year_2022,
                    interp_year_2023,
                    interp_year_2024,
                    par=TRUE)
names(interp_mean) <- gsub("LAI-2016-","LAI-",names(interp_year_2016))
writeRaster(interp_mean,filename="./output/lai_predictions_HLS_interp_mean_2015_2024.tif",
            overwrite=TRUE)

winter_mean <- mean(interp_mean[[c(1:20)]],
                    filename="./output/lai_winter_mean_2015_2023.tif",
                    overwrite=TRUE)

peak_mean <- mean(interp_mean[[c(91:151)]],
                  filename="./output/lai_peak_mean_2015_2023.tif",
                  overwrite=TRUE)

gs_amplitude_mean <- max(interp_mean[[26:66]]) - min(interp_mean[[26:66]])
writeRaster(gs_amplitude_mean,filename="./output/lai_gs_amplitude_mean_2015_2023.tif",
            overwrite=TRUE)

lai_mask <- winter_mean < 2000 & gs_amplitude_mean > 200

crop_ext <- ext(315310, 337181, 4300236, 4322160)

anom_fun <- function(x){
  out <- mask(x - interp_mean,lai_mask,maskvalues=0,updatevalue=NA)
}

interp_list <- list(interp_year_2015,
                    interp_year_2016,
                    interp_year_2017,
                    interp_year_2018,
                    interp_year_2019,
                    interp_year_2020,
                    interp_year_2021,
                    interp_year_2022,
                    interp_year_2023,
                    interp_year_2024)

anom_list <- list()
anom_years <- 2015:2024
for(i in 1:length(interp_list)){
  anom_list[[i]] <- anom_fun(interp_list[[i]])
  names(anom_list[[i]]) <- gsub("2016",anom_years[i],names(interp_year_2016))
}

integral_maps <- list.files("./output/",pattern="integral",full.names=TRUE)
integral_rast <- rast(integral_maps)
integral_mean <- mask(mean(integral_rast),lai_mask,
                      maskvalues=0,updatevalue=NA)
integral_anom <- integral_rast - integral_mean
integral_prop_anom <- integral_rast / integral_mean
names(integral_anom) <- anom_years
names(integral_prop_anom) <- anom_years
writeRaster(integral_prop_anom,filename="./output/lai_int_proportional_anoms_2015_2024.tif",
            overwrite=TRUE)
writeRaster(integral_anom,filename="./output/lai_int_anoms_2015_2024.tif",
            overwrite=TRUE)

## Extracts time-series for leafy landscapes pixels.
drone_df <- read.csv("./output/upscaled_lai_q50_predictions_30m_mod5.csv")
drone_df$Date <- as.Date(drone_df$Date, format="%Y-%m-%d")
drone_sites <- drone_df[!(duplicated(drone_df$PointID)),]
drone_sites_samp <- drone_sites %>% group_by(Site) %>% slice_sample(n=100)

site_extract_fun <- function(x){
  ext_wide <- cbind(drone_sites_samp$PointID,
                    drone_sites_samp$Site,
                    drone_sites_samp$x,
                    drone_sites_samp$y,
                    terra::extract(x,as.matrix(cbind(drone_sites_samp$x,drone_sites_samp$y)),
                                   method="simple"))
  ext_long <- pivot_longer(ext_wide,cols=contains("LAI"))
  colnames(ext_long) <- c("PointID","Site","x","y","year_DOY","HLS_LAI")
  return(ext_long)
}
site_extract_allyears <- lapply(interp_list,site_extract_fun)
site_extract_means <- site_extract_fun(interp_mean)
site_extract_sitemeans <- site_extract_means %>% group_by(Site,year_DOY) %>% summarise(HLS_LAI_mean=mean(HLS_LAI))
site_extract_sitemeans$DOY <- as.numeric(gsub("LAI-doy","",site_extract_sitemeans$year_DOY))


site_extract_df <- bind_rows(site_extract_allyears)
site_extract_df$Year <- stringr::str_split_fixed(site_extract_df$year_DOY,pattern="-",n=3)[,2]
site_extract_df$DOY <- as.numeric(gsub("doy","",stringr::str_split_fixed(site_extract_df$year_DOY,pattern="-",n=3)[,3]))
#site_extract_df$Year <- factor(site_extract_df$Year,levels=c(2018,2015,2020,2021,2022,2016,2017,2023,2019))
site_extract_df$Site <- factor(site_extract_df$Site,
                               levels=c("JuddFalls","LupineTrail","ParadiseBasin","StrandHill","SouthBaldy",
                                        "SnodgrassTrailhead","SlideRockRoad", "RustlersGulch"))
point_mean_df <- site_extract_df %>% group_by(Site,Year,DOY) %>% summarise(HLS_LAI_year_mean=mean(HLS_LAI))
site_extract_df <- left_join(site_extract_df,site_extract_sitemeans,by=c("Site","DOY"))
site_extract_df <- left_join(site_extract_df,point_mean_df, by=c("Site","Year","DOY"))


drone_site_lai_plot <- ggplot(filter(site_extract_df,DOY > 120 & DOY < 300))+
  geom_line(aes(x=DOY,y=HLS_LAI/1000,group=interaction(Year,PointID)),lwd=0.2,alpha=0.1)+
  geom_line(aes(x=DOY,y=HLS_LAI_mean/1000,group=interaction(Year,PointID),linetype="2015-2024 Mean"))+
  geom_line(aes(x=DOY,y=HLS_LAI_year_mean/1000,group=interaction(Year,PointID),linetype="Site Yearly Mean"))+
  scale_linetype_manual("", values=c("Site Yearly Mean"=1,
                                     "2015-2024 Mean"=3))+
  scale_y_continuous("HLS Estimated LAI",limits=c(0,5.5))+
  facet_grid(facets=Year~Site)+
  theme_bw()+
  theme(legend.position="bottom")

pdf("./figs/drone_sites_HLS_LAI_2015_2024.pdf",width=12,height=10)
print(drone_site_lai_plot)
dev.off()

## Plots integrated LAI anomalies.
integral_anom_plot <- ggplot()+
  geom_spatraster(data=integral_anom)+
  scale_fill_distiller(name="Integated\nLAI Anomaly",
                       limits=c(-150,150),palette="RdBu",direction = 1)+
  facet_wrap(facets="lyr",ncol=5)+
  theme_minimal()

pdf("./figs/lai_integral_anoms_2015_2024.pdf",width=11,height=6)
print(integral_anom_plot)
dev.off()

prop_anom_plot <- ggplot()+
  geom_spatraster(data=integral_prop_anom)+
  scale_fill_distiller(name="Proportional\nInt. LAI \nAnomaly",
                       limits=c(0.5,1.5),palette="RdBu",direction = 1)+
  facet_wrap(facets="lyr",ncol=5)+
  theme_minimal()

pdf("./figs/lai_integral_prop_anoms_2015_2024.pdf",width=11,height=6)
print(prop_anom_plot)
dev.off()


plot_doy <- as.character(c(151,181,211,241,271))
plot_layers <- which(substr(names(interp_year_2016),13,15) %in% plot_doy) 

lai_anom_plots <- list()
plot_years <- 2015:2024
first_anoms <- anom_list[[1]][[plot_layers]]
label_dates <- format(as.Date(paste0("2015-",plot_doy),format="%Y-%j"),
                      format="%B %d")
names(first_anoms) <- label_dates
lai_anom_plots[[1]] <- ggplot()+
  geom_spatraster(data=first_anoms)+
  scale_fill_distiller(limits=c(-1500,1500),palette="RdBu",direction = 1)+
  scale_y_continuous(as.character(plot_years[1]),
                     breaks=c(38.8,38.9,39.0,39.1))+
  scale_x_continuous(breaks=c(-107.15,-107.0,-106.85))+
  facet_wrap(facets=~lyr,ncol=5)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

for(i in 2:(length(plot_years)-1)){
  print(paste("Creating plot for year",plot_years[i]))
  lai_anom_year <- anom_list[[i]][[plot_layers]]
  lai_anom_plots[[i]] <- ggplot()+
    geom_spatraster(data=lai_anom_year)+
    scale_fill_distiller(limits=c(-1500,1500),palette="RdBu",direction = 1)+
    scale_y_continuous(as.character(plot_years[i]),
                       breaks=c(38.8,38.9,39.0,39.1))+
    scale_x_continuous(breaks=c(-107.15,-107.0,-106.85))+
    facet_wrap(facets=~lyr,ncol=5)+
    theme_minimal()+
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())
}

last_anoms <- anom_list[[length(plot_years)]][[plot_layers]]

lai_anom_plots[[length(plot_years)]] <- ggplot()+
  geom_spatraster(data=last_anoms)+
  scale_fill_distiller(name="LAI Anomaly", limits=c(-1500,1500),
                       palette="RdBu", direction = 1,
                       breaks=c(-1500,-1000,-500,0,500,1000,1500),
                       labels=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5))+
  scale_y_continuous(as.character(plot_years[length(plot_years)]),
                     breaks=c(38.8,38.9,39.0,39.1))+
  scale_x_continuous(breaks=c(-107.15,-107.0,-106.85))+
  facet_wrap(facets=~lyr,ncol=5)+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.key.width = unit(2,"cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank())

pdf("./figs/lai_anoms_2015_2024.pdf",width=12,height=22)
gridExtra::grid.arrange(lai_anom_plots[[1]],
                        lai_anom_plots[[2]],
                        lai_anom_plots[[3]],
                        lai_anom_plots[[4]],
                        lai_anom_plots[[5]],
                        lai_anom_plots[[6]],
                        lai_anom_plots[[7]],
                        lai_anom_plots[[8]],
                        lai_anom_plots[[9]],
                        lai_anom_plots[[10]],
                        heights=c(1.12,rep(1,8),1.42),
                        ncol=1)
dev.off()

## Creates an animated GIF map of LAI anomalies across the season
for(i in 61:65){
  
  print(paste("Creating animation frame ", i))
  p1 <- ggplot()+
    labs(subtitle=names(interp_list[[4]][[i]]))+
    geom_spatraster(data=interp_list[[4]][[i]])+
    scale_fill_whitebox_c(limits=c(0,5000),palette="viridi")+
    #facet_wrap(facets=~lyr)+
    theme_minimal()+
    theme(legend.position = "none")
  
  p2 <- ggplot()+
    labs(subtitle=names(interp_list[[5]][[i]]))+
    geom_spatraster(data=interp_list[[5]][[i]])+
    scale_fill_whitebox_c(limits=c(0,5000),palette="viridi")+
    #facet_wrap(facets=~lyr)+
    theme_minimal()+
    theme(legend.position = "none")
  
  p3 <- ggplot()+
    labs(subtitle=names(interp_list[[10]][[i]]))+
    geom_spatraster(data=interp_list[[10]][[i]])+
    scale_fill_whitebox_c(name="LAI",limits=c(0,5000),breaks=c(0,2000,4000),
                          labels=c(0,2,4),palette="viridi")+
    #facet_wrap(facets=~lyr)+
    theme_minimal()
  
  png(paste0("./figs/animation_frames/lai_frame_",i,".png"),width=11,height=6,units="in",res=300)
  gridExtra::grid.arrange(p1, p2, p3,ncol=3,widths=c(1,1,1.2))
  dev.off()
}
