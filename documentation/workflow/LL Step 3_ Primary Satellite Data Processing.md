# LL Step 3: Primary Satellite Data Processing

Note: You can see how this step fits into the bigger picture of the analysis here: [Leafy Landscapes LAI Workflow](file:///workspace/e094b233-4dce-4599-bc85-a12ab329bda3/7JKZ-alrJKeD8D00oybKY?mode=edgeless\&blockIds=FRYfcvTHOD0CAGHr69rcJ)

## Overview

The objective of this processing step is to find, filter and prepare satellite image data for downstream analysis. The primary source of satellite data is the Harmonized Landsat-Sentinel 2 Surface Reflectance Product Version 2.0 ("HLS"). The HLS data is produced through a collaboration between the US Geological Survey, NASA, and the European Space Agency and takes data from four Earth-observing satellites: Landsat 8 and 9 as well as Sentinel 2a and 2b. HLS data has 30 m resolution observations every 2-3 days from 2015 - present, but many scenes are obscured by clouds. We used the NASA web tool AppEEARS to get a list of scenes intersecting the study area that have <50% cloud cover. We then download the filtered scenes, crop them to the study area, and use the included quality control information to mask out pixels that are obscured by clouds or cloud shadows. We then compute spectral indices for each scene that we can later use to predict LAI.

## Dataset Descriptions

### Primary Datasets

##### 2024 HLS Satellite Scenes

The original dataset is surface reflectance at a set of spectral bands that are in common for the OLI instrument on Landsat 8/9 and the MSI instrument on Sentinel 2a / 2b. These data have been processed to be on the same 30 x 30 m spatial grid and have similar radiometric calibration. For the first stage of the analysis (relating drone data to satellite data) we just need images from 2024 for a small area around the field plots, but later on we will also need data from the full time-period of the dataset (2015 - 2024) for a larger area that we want to predict LAI. These data are stored as single-layer Cloud Optimized Geotiff (COG) raster files publicly available on NASA's Land Processes DAAC site. There is one file for each spectral band, as well as an additional COG dataset containing quality control information.

### Analysis Steps

##### Query, Filter and Download HLS Scenes

In the first step, we use the NASA web-based tool AppEEARS to find all HLS scenes that overlap with the study region and have < 50% cloud cover. We then use the tool to extract and download the data from the NASA Land Processes DAAC repository. We then deposit the cropped HLS scenes in the project OneDrive folder. This step is performed manually by navigating the AppEEARS website.

##### Cloud-mask and Align HLS Scenes

The HLS satellite data is on a consistent 30 m x 30 m grid with each spectral band for each scene stored as a separate COG raster file along with quality-control information. We have a script (Step3a\_Prep\_HLS\_Data.R) which automates the process of masking out clouds and aligning the raster files together.

##### Extract and interpolate spectral information for drone flight areas.

To relate to the drone-based estimates, we need to extract spectral index information for each pixel that is covered by the drone data. Also, because the satellite overpasses are often not on the same day as the drone flights or field data collection, we also need to interpolate the indices to get values that correspond to the exact day of field sampling. A section in the script Step3a\_Prep\_HLS\_Data.R accomplishes this step.

### Derived Datasets

##### Interpolated Spectral Index Values for Drone Flight Areas

This dataset is a single table (2.3 MB) that stores the interpolated reflectance and spectral index values for each 30m pixel covered by the drone dataset on the days of field sampling.

### Notes on Data Storage and Computing

In this processing step, we just need to consider satellite data from 2024 over a small region that encompasses all of the field sites. The NASA AppEEARS system provides the ability to extract just the data we need for this. It's about 15 GB in total and we store it locally in the project OneDrive folder. The processing script takes a few minutes to complete on a laptop computer.
