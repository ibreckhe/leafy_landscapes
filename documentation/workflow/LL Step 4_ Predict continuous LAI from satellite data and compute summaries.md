# LL Step 4: Predict continuous LAI from satellite data and compute summaries

Note: You can see how this step fits into the bigger picture of the analysis here: [Leafy Landscapes LAI Workflow](file:///workspace/e094b233-4dce-4599-bc85-a12ab329bda3/7JKZ-alrJKeD8D00oybKY?mode=edgeless\&blockIds=_regIoMN8N4PDItdjhHrc)

## Overview

In previous steps, we prepared upscaled drone-based estimates of LAI at the same resolution as the HLS sattelite data. We also prepared a set of satellite-based observations that we can use to predict LAI. In this step, we fit a machine-learning model to predict the upscaled LAI using satellite measurements and then apply this model to a large archive of HLS satellite imagery from 2015 - 2024. The resulting LAI estimates for each cloud-free pixel are then temporally interpolated to create continuous gap-free LAI estimates for the entire growing season in each year. Finally, we compute a set of final data products that represent statistical summaries of the interpolated LAI.

## Dataset Descriptions

### Primary Datasets

##### 2015 - 2024 HLS Satellite Scenes

Similar to the primary data from Step 3, the primary dataset for this analysis step are scenes from the Harmonized Landsat-Sentinel (HLS) dataset. In contrast to Step 3, however, we want to process data from a much larger portion of the Gunnison Basin and do so for all the years of available data (2015 - 2024). These data are stored as single-layer Cloud Optimized Geotiff (COG) raster files publicly available on NASA's Land Processes DAAC site, but the full set of scenes are prohibitively large to download locally (\~500 GB). Instead, we use the NASA web tool AppEEARS to get the URLs of the scenes of interest and use the tools in the R package 'terra' to download, extract and process only the parts of the scenes we need.

### Analysis Steps

##### Fit a machine-learning model to predict LAI from HLS satellite data.

We have about 5000 30 m pixels where we have derived upscaled LAI estimates from our previous work, matched up in space and time with spectral indices from the HLS data. We developed a script in R (Step4a\_Fit\_HLS\_LAI\_model.R) to relate the HLS spectral data to the upscaled field-based LAI. Currently, this model is a Generalized Additive Mixed Model (GAMM) fit with the R package 'mgcv'.

##### Prep HLS data and apply the fit ML model

With the fit LAI model, we can then apply it to all the HLS scenes of interest. Because the data is large, we developed optimized code (Step4b\_Apply\_HLS\_LAI\_model.R) that prepares the raw observations and predicts LAI values in a single step. This script completes in several hours on a desktop workstation.

##### Temporal interpolation and smoothing

The resulting LAI predictions for each satellite scene are irregular in time and contain gaps caused by clouds, shadows, and other issues. Temporal interpolation can reconstruct continuous values at a regular interval. We have developed a script (Step4c\_Interpolate\_HLS\_LAI.R) that performs this temporal smoothing to generate continuous daily LAI estimates for each pixel. This script completes in about 18 hours on a desktop workstation.

##### Compute Statistical Summaries

The final quantities of interest from this analysis are seasonally aggregated measures of LAI, including the seasonal maximum (peak LAI), time-integrated LAI, and measures of phenological timing (the day-of-year of 50% green-up and brown-down. The final steps of the script Step4c\_Interpolate\_HLS\_LAI.R perform these operations to compute summaries.

### Derived Datasets

##### Continuous LAI Estimates

The final processing steps that generate aggregated LAI measures for each year require LAI data to be continuous in space (no gaps) and equally spaced in time (e.g. daily). Our processing script generates these interpolated LAI datasets as COG files with many layers, one for each interpolated daily map. 

### Notes on Data Storage and Computing

The primary dataset for this processing step (the HLS scenes from 2015 - 2021) are extremely large (\~500 GB), making it impractical to download the data before processing. Instead, we use a strategy enabled by the data format (Cloud-optimized Geotiff, COG), to extract and work with just the data that we need. We store intermediate datasets locally in the project OneDrive folder. This part of the analysis is currently run locally and takes several days on a desktop workstation. It could be sped up significantly by using a tiled workflow on cloud computing which processes many tiles in parallel. Unfortunately, this dataset is not already on commerical cloud storage, so there would still be a significant bottleneck in data transfer.

