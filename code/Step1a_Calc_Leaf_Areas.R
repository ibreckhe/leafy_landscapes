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

