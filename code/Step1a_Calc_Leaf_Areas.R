## Script to calculate leaf area from scanned images of leaves.

## Load necessary libraries
library(terra)

## Set working directory
setwd("/Users/ian/Library/CloudStorage/OneDrive-RMBL/Documents\ -\ Research\ -\ Spatial\ Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/data/leaf_scans/")

## List image files.
img_files <- list.files(".",pattern = ".jpg$")

## Tests function logic
# img1 <- rast(img_files[1])
# hsv_fun <- function(x) {
#   as.numeric(rgb2hsv(matrix(c(x[1], x[2], x[3]), nrow=3, ncol=1) ))
# }
# rhsv <- app(img1, hsv_fun)
# 
# thresh_fun <- function(x){
#   ifelse(x[1] > 0.01 & x[1] < 0.4 &
#          x[2] > 0.001 & x[2] < 1 &
#          x[3] > 0 & x[3] < 0.7, 1, 0)
# }
# img_thresh <- app(rhsv, thresh_fun)
# 
# plotRGB(img1)
# plot(img_thresh, add=TRUE,alpha=0.5, col=c(1,2))

## Function to calculate leaf area
calc_leaf_area <- function(img_path,
                           out_path="../leaf_masks/",
                           plot_path="../leaf_plots/",
                           dpi=300,
                           h_thresh=c(0.01,0.4), 
                           s_thresh=c(0,1), 
                           v_thresh=c(0,0.7),
                           write_output=TRUE,
                           overwrite=TRUE){
  
  print(paste("Now processing: ", img_path))
  out_name <- paste(out_path,gsub(".jpg$","_thresh.tif",img_path),sep="")
  
  if(file.exists(out_name) & overwrite==FALSE){
    print(paste("Output file already exists. Skipping..."))
    return(NULL)
  }
  
  img <- rast(img_path)
  
  print(paste("Converting to HSV color space..."))
  
  hsv_fun <- function(x) {
    as.numeric(rgb2hsv(matrix(c(x[1], x[2], x[3]), nrow=3, ncol=1) ))
  }
  img_hsv <- app(img, hsv_fun)
  
  print(paste("Thresholding image..."))
  
  thresh_fun <- function(x){
    ifelse(x[1] >= h_thresh[1] & x[1] <= h_thresh[2] &
           x[2] >= s_thresh[1] & x[2] <= s_thresh[2] &
           x[3] >= v_thresh[1] & x[3] <= v_thresh[2], 1, 0)
  }
  img_thresh <- app(img_hsv, thresh_fun)
  
  img_area <- global(img_thresh,"sum") * ((1/dpi)*(1/dpi)) * 6.4516
  print(paste("Estimated leaf area: ", img_area, "cm^2"))
  area_df <- data.frame(file=img_path,leaf_area_sqcm=img_area)
  colnames(area_df) <- c("file","leaf_area_sqcm")
  
  plot_name <- paste(plot_path,gsub(".jpg$","_overlay.pdf",img_path),sep="")
  
  par(oma=c(2,2,2,2),mar=c(1,1,1,1))
  pdf(plot_name,width=8.5,height=11)
  plotRGB(img,axes=TRUE,main=img_path)
  plot(img_thresh, add=TRUE,alpha=0.5, col=c(1,2))
  dev.off()
  
  ## Writing outputs to file.
  if(write_output==TRUE){
    writeRaster(img_thresh, out_name,
                overwrite=overwrite)
  }
  return(list(img_thresh,area_df))
}

## Tests function.
img_test <- calc_leaf_area(img_files[450],
                           out_path="../leaf_masks/",
                           plot_path="../leaf_plots/",
                           dpi=300,
                           h_thresh=c(0.0,1), 
                           s_thresh=c(0.05,1), 
                           v_thresh=c(0,0.75),
                           write_output=FALSE,
                           overwrite=TRUE)

## Loop through all images.
out_list <- list()
area_list <- list()
for(i in 1:length(img_files)){
  out <- calc_leaf_area(img_files[i],
                             out_path="../leaf_masks/",
                             plot_path="../leaf_plots/",
                             dpi=300,
                             h_thresh=c(0,1), 
                             s_thresh=c(0.05,1), 
                             v_thresh=c(0,0.75),
                             write_output=TRUE,
                             overwrite=TRUE)
  if(!is.null(out)){
    out_list[[i]] <- out
    area_list[[i]] <- out_list[[i]][[2]]
  }
}


## Combines outputs into a single data frame and writes to disk.
library(dplyr)
area_df <- bind_rows(area_list)

area_df$h_thresh_1 <- 1
area_df$h_thresh_2 <- 1
area_df$s_thresh_1 <- 0.05
area_df$s_thresh_2 <- 1
area_df$v_thresh_1 <- 0
area_df$v_thresh_2 <- 0.75

write.csv(area_df,file="../scan_leaf_areas_8_8_2024.csv")

## Sums leaf area for all images in a sample.
library(stringr)
area_df$plot <- toupper(str_split_fixed(area_df$file,"-",n=3)[,1])
area_df$sample <- tolower(str_split_fixed(area_df$file,"-",n=3)[,2])
area_df$scan <- gsub(".jpg","",str_split_fixed(area_df$file,"-",n=3)[,3])

area_df_sum <- area_df %>% 
  group_by(plot,sample) %>%
  summarise(total_leaf_area_sqcm=sum(leaf_area_sqcm),
            n_scans=n())

write.csv(area_df_sum,file="../plot_leaf_areas_8_8_2024.csv",row.names=FALSE)

## Joins with rest of plot data in Google Sheets.
library(googlesheets4)
gs4_auth(email="ibreckhe@gmail.com")

plot_data <- read_sheet("https://docs.google.com/spreadsheets/d/1rSPt5g_-af3EsmYj_nXwJO9Ya9uL4duAkrBlrsVNMSo",
                        range="PL_data!A4:AB317")
plot_data$PlotID <- toupper(plot_data$PlotID)
plot_data$Num_leaves_scanned <- as.numeric(plot_data$Num_leaves_scanned)

plot_data_area <- plot_data %>% 
  select(PlotID,Species_or_FT,Collected_date,Num_leaves_scanned,Scanned_wet_mass,Scanned_dry_mass,Scanned_leaf_water_pct) %>%
  left_join(area_df_sum,by=c("PlotID"="plot",
                             "Species_or_FT"="sample"))
plot_data_area$LMA <- plot_data_area$Scanned_dry_mass / (plot_data_area$total_leaf_area_sqcm * 0.0001)

## Summarises LMA for each functional type.
library(ggplot2)

p1 <- ggplot(filter(plot_data_area,Species_or_FT %in% c("forb","grass","shrub")))+
  geom_point(aes(x=Species_or_FT,y=LMA,color=Species_or_FT), 
             position=position_jitter(width=0.2),alpha=0.1)+
  geom_text(aes(x=Species_or_FT,y=LMA,label=PlotID, color=Species_or_FT), 
            position=position_jitter(width=0.2),size=2)+
  scale_y_continuous("CWM Leaf Mass Per Area [g/m^2]",limits=c(0,420))+
  theme_bw()+
  ggtitle("Functional Group LMA")

p2 <- ggplot(filter(plot_data_area,Species_or_FT %in% c("forb","grass","shrub")))+
  geom_point(aes(x=Scanned_leaf_water_pct,y=LMA,color=Species_or_FT),alpha=0.1)+
  geom_smooth(aes(x=Scanned_leaf_water_pct,y=LMA,color=Species_or_FT),method="lm",se=TRUE)+
  geom_text(aes(x=Scanned_leaf_water_pct,y=LMA,label=PlotID, color=Species_or_FT),size=2)+
  scale_y_continuous("CWM Leaf Mass Per Area [g/m^2]",limits=c(0,420))+
  scale_x_continuous("CWM Leaf Water Content [%]")+
  theme_bw()+
  ggtitle("LMA vs. Leaf Water Content")

pdf("../LMA_plots.pdf",width=8.5,height=11)
gridExtra::grid.arrange(p1,p2,ncol=1)
dev.off()
