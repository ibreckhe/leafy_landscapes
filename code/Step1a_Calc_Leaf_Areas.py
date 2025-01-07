# Import necessary libraries
import os
import cv2
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from skimage import color, io

# Set working directory to where the leaf scan images are located
os.chdir("/Users/ian/Library/CloudStorage/OneDrive-RMBL/Documents - Research - Spatial Ecology/General/SpatialEcologyShared/Projects/Meadow_LAI/data/leaf_scans/")

# List all image files in the current directory with the .jpg extension
img_files = [f for f in os.listdir() if f.endswith('.jpg')]

# Define a function to calculate leaf area from an image file
def calc_leaf_area(img_path, out_path = "../leaf_masks/", plot_path = "../leaf_plots/", dpi = 300, h_thresh = (0.01, 0.4), s_thresh = (0, 1), v_thresh = (0, 0.7), write_output = True, overwrite = True):
    # Define the output file name
    out_name = out_path + img_path.replace(".jpg", "_thresh.tif")

    # If the output file already exists and overwrite is set to False, skip processing this image
    if os.path.exists(out_name) and not overwrite:
        print("Output file already exists. Skipping...")
        return None

    # Read the image file into a numpy array
    img = io.imread(img_path)

    # Convert the image from RGB to HSV color space
    img_hsv = color.rgb2hsv(img)

    # Threshold the image based on the HSV values and the specified thresholds
    img_thresh = np.zeros_like(img_hsv[:,:,0])
    img_thresh[(img_hsv[:,:,0] >= h_thresh[0]) & (img_hsv[:,:,0] <= h_thresh[1]) &
               (img_hsv[:,:,1] >= s_thresh[0]) & (img_hsv[:,:,1] <= s_thresh[1]) &
               (img_hsv[:,:,2] >= v_thresh[0]) & (img_hsv[:,:,2] <= v_thresh[1])] = 1

    # Calculate the total leaf area in square centimeters
    img_area = np.sum(img_thresh) * ((1 / dpi) ** 2) * 6.4516
    print(f"Estimated leaf area: {img_area} cm^2")

    # Create a data frame with the file name and leaf area
    area_df = pd.DataFrame({'file': [img_path], 'leaf_area_sqcm': [img_area]})

    # Define the output file name for the plot
    plot_name = plot_path + img_path.replace(".jpg", "_overlay.pdf")

    # Create a plot of the original image with the thresholded image overlaid
    plt.figure(figsize=(8.5, 11))
    plt.imshow(img)
    plt.imshow(img_thresh, alpha=0.5, cmap='gray')
    plt.savefig(plot_name)
    plt.close()

    # Write the thresholded image to a file
    if write_output:
        io.imsave(out_name, img_thresh)

    # Return the thresholded image and the data frame with the leaf area
    return {'img_thresh': img_thresh, 'area_df': area_df}

# Tests function with the 450th image in the list
img_test = calc_leaf_area(img_files[450],
                          out_path="../leaf_masks/",
                          plot_path="../leaf_plots/",
                          dpi=300,
                          h_thresh=(0.0,1),
                          s_thresh=(0.05,1),
                          v_thresh=(0,0.75),
                          write_output=False,
                          overwrite=True)

# Loop through all images, calculate leaf area for each image, and store the results in a list
out_list = []
area_list = []
for i in range(len(img_files)):
    out = calc_leaf_area(img_files[i])
    if out is not None:
        out_list.append(out)
        area_list.append(out['area_df'])

# Combine the results for all images into a single data frame, and add the threshold values as columns
area_df = pd.concat(area_list)
area_df['h_thresh_1'] = h_thresh[0]
area_df['h_thresh_2'] = h_thresh[1]
area_df['s_thresh_1'] = s_thresh[0]
area_df['s_thresh_2'] = s_thresh[1]
area_df['v_thresh_1'] = v_thresh[0]
area_df['v_thresh_2'] = v_thresh[1]

# Write the data frame to a CSV file
area_df.to_csv("../scan_leaf_areas_8_8_2024.csv", index=False)

# Summarize the leaf area for all images in a sample by grouping the data frame by plot and sample, and calculating the sum and count of leaf area
area_df['plot'] = area_df['file'].str.split('-').str[0].str.upper()
area_df['sample'] = area_df['file'].str.split('-').str[1].str.lower()
area_df['scan'] = area_df['file'].str.split('-').str[2].str.replace('.jpg', '')

area_df_sum = area_df.groupby(['plot', 'sample']).agg({'leaf_area_sqcm': ['sum', 'count']})
area_df_sum.columns = ['total_leaf_area_sqcm', 'n_scans']

# Write the summarized data frame to a CSV file
area_df_sum.to_csv("../plot_leaf_areas_8_8_2024.csv")