## Script to calculate leaf area from scanned images of leaves.
## Author: Ian Breckheimer
## Updated: 1 January 2025

## This script should be run after all samples are scanned and saved as .jpg files. To see how this script fits into the larger analysis, see the documentation here: https://app.affine.pro/workspace/e094b233-4dce-4599-bc85-a12ab329bda3/lbgAKw_fYkaoE7YpIsfgv?mode=page

## Load necessary libraries
# 'terra' for raster data handling, 'dplyr' for data manipulation, and 'stringr' for string manipulation
library(terra)
library(dplyr)
library(stringr)

# Set working directory to where the leaf scan images are located
setwd("/Users/ian/Library/CloudStorage/OneDrive-RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/data/leaf_scans/")

# List all image files in the current directory with the .jpg extension
img_files <- list.files(".", pattern = ".jpg$")

# Define a function to calculate leaf area from an image file
# This function takes the image path as input, along with various parameters for image thresholding and output settings
calc_leaf_area <- function(img_path, out_path = "../leaf_masks/", plot_path = "../leaf_plots/", dpi = 300, h_thresh = c(0.01, 0.4), s_thresh = c(0, 1), v_thresh = c(0, 0.7), write_output = TRUE, overwrite = TRUE) {
  # Define the output file name
  out_name <- paste0(out_path, gsub(".jpg$", "_thresh.tif", img_path))

  # If the output file already exists and overwrite is set to FALSE, skip processing this image
  if (file.exists(out_name) & !overwrite) {
    message("Output file already exists. Skipping...")
    return(NULL)
  }

  # Read the image file into a raster object
  img <- rast(img_path)

  # Convert the image from RGB to HSV color space
  hsv_fun <- function(x) {
    as.numeric(rgb2hsv(matrix(x, nrow = 3, ncol = 1)))
  }
  img_hsv <- app(img, hsv_fun)

  # Threshold the image based on the HSV values and the specified thresholds
  thresh_fun <- function(x) {
    ifelse(x[1] >= h_thresh[1] & x[1] <= h_thresh[2] & x[2] >= s_thresh[1] & x[2] <= s_thresh[2] & x[3] >= v_thresh[1] & x[3] <= v_thresh[2], 1, 0)
  }
  img_thresh <- app(img_hsv, thresh_fun)

  # Calculate the total leaf area in square centimeters
  img_area <- global(img_thresh, "sum") * ((1 / dpi) ^ 2) * 6.4516
  message(paste("Estimated leaf area:", img_area, "cm^2"))

  # Create a data frame with the file name and leaf area
  area_df <- data.frame(file = img_path, leaf_area_sqcm = img_area)

  # Define the output file name for the plot
  plot_name <- paste0(plot_path, gsub(".jpg$", "_overlay.pdf", img_path))

  # Create a plot of the original image with the thresholded image overlaid
  pdf(plot_name, width = 8.5, height = 11)
  plotRGB(img, axes = TRUE, main = img_path)
  plot(img_thresh, add = TRUE, alpha = 0.5, col = c(1, 2))
  dev.off()

  # Write the thresholded image to a file
  if (write_output) {
    writeRaster(img_thresh, out_name, overwrite = overwrite)
  }

  # Return the thresholded image and the data frame with the leaf area
  return(list(img_thresh = img_thresh, area_df = area_df))
}

## Tests function with the 450th image in the list
img_test <- calc_leaf_area(img_files[450],
                           out_path="../leaf_masks/",
                           plot_path="../leaf_plots/",
                           dpi=300,
                           h_thresh=c(0.0,1),
                           s_thresh=c(0.05,1),
                           v_thresh=c(0,0.75),
                           write_output=FALSE,
                           overwrite=TRUE)

# Loop through all images, calculate leaf area for each image, and store the results in a list
out_list <- list()
area_list <- list()
for (i in seq_along(img_files)) {
  out <- calc_leaf_area(img_files[i])
  if (!is.null(out)) {
    out_list[[i]] <- out
    area_list[[i]] <- out$area_df
  }
}

# Combine the results for all images into a single data frame, and add the threshold values as columns
area_df <- bind_rows(area_list)
area_df$h_thresh_1 <- h_thresh[1]
area_df$h_thresh_2 <- h_thresh[2]
area_df$s_thresh_1 <- s_thresh[1]
area_df$s_thresh_2 <- s_thresh[2]
area_df$v_thresh_1 <- v_thresh[1]
area_df$v_thresh_2 <- v_thresh[2]

# Write the data frame to a CSV file
write.csv(area_df, "../scan_leaf_areas_8_8_2024.csv", row.names = FALSE)

# Summarize the leaf area for all images in a sample by grouping the data frame by plot and sample, and calculating the sum and count of leaf area
area_df$plot <- toupper(str_split_fixed(area_df$file, "-", n = 3)[, 1])
area_df$sample <- tolower(str_split_fixed(area_df$file, "-", n = 3)[, 2])
area_df$scan <- gsub(".jpg", "", str_split_fixed(area_df$file, "-", n = 3)[, 3])

area_df_sum <- area_df %>%
  group_by(plot, sample) %>%
  summarise(total_leaf_area_sqcm = sum(leaf_area_sqcm),
            n_scans = n())

# Write the summarized data frame to a CSV file
write.csv(area_df_sum, "../plot_leaf_areas_8_8_2024.csv", row.names = FALSE)
  summarise(total_leaf_area_sqcm = sum(leaf_area_sqcm), n_scans = n())
write.csv(area_df_sum, "../plot_leaf_areas_8_8_2024.csv", row.names = FALSE)