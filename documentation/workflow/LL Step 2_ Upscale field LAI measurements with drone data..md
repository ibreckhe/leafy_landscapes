# LL Step 2: Upscale field LAI measurements with drone data.

Note: You can see how this step fits into the bigger picture of the analysis here: [Leafy Landscapes LAI Workflow](file:///workspace/e094b233-4dce-4599-bc85-a12ab329bda3/7JKZ-alrJKeD8D00oybKY?mode=edgeless\&blockIds=OfOAMVb84z6fWJS4s59kW)

## Overview

The objective of this processing step is to estimate landscape patterns of LAI using the drone imagery that we collect at each of our field sites. This will allow us to scale up from the field data (collected in 1 m^2 plots) to the scale of the satellite imagery (30 x 30 m pixels). This requires us to (1) extract drone data for the exact locations where we have field data, (2) process the drone data to develop measurements that can use to reliably predict LAI ("predictors"), (3) develop a statistical model that related these predictors to the field data, and (4) apply that model to the full drone mosaic to create continuous maps of LAI.

## Dataset Descriptions

### Primary Datasets

##### Geolocations of Plot Boundaries

Most of the datasets for this step in processing are derived data (i.e. they are created through previous processing steps). The only primary dataset that enters the analysis here is one that documents the geolocations of the corners of each of the field sampling plots. These data were collected with a Trimble Catalyst DA2 GNSS system and ESRI Field Maps software. After the completion of field work, the data were then QC'd in ArcGIS Pro and exported as polygons to a Geopackage file.

### Analysis Steps

##### Extract drone data for each plot polygon.

In this step, we are using the geolocations of the field plots to extract "chips" of the drone data for areas that exactly correspond to the areas that were subsequently field sampled. We are then computing a set of derived metrics for each chip (spectral and structure indices) that are then related to field data to produce an LAI estimate for each drone image. We have an R script that automates this procedure (Step2a\_Extract\_Drone\_Predictors.R)

##### Statistical model for LAI based on field and drone data.

Using the field-based estimates of LAI for each plot and the extracted and transformed drone data, we fit a statisitcal model that predicts the field-based LAI using drone data. The current version of this model is a Bayesian Generalized Linear Mixed Model fit with the 'rstanarm' R package. We have an R script that automates this procedure (Step2b\_Drone\_LAI\_Models.R). This script fits a variety of models and compares them, selecting the best model for downstream work. It then predicts LAI values at 1m resolution for each drone image using the best model. 

##### Filter and upscale drone-based LAI estimates.

The HLS satellite data is on a consistent 30 m x 30 m grid, and we need comparable LAI estimates derived from the field and drone data for downstream analysis. We aggregate and filter drone-based predictions from the previous analysis step to form a collection of LAI estimates at the native resolution of the satellite data. We have an R script that automates this procedure (Step2c\_Upscale\_Drone\_LAI.R).

##### Derived Datasets

##### Leaf Area and Mean Traits for Each Plot

For each field plot, we have a reliable estimates of the LAI as well as a few other attributes (aboveground biomass, canopy water content, community mean LMA) that are also of interest. This data is stored as a CSV table in the project OneDrive folder. This is the primary output from the previous processing steps (see [LL Step 1: Assemble and QC field and laboratory data.](file:///workspace/e094b233-4dce-4599-bc85-a12ab329bda3/lbgAKw_fYkaoE7YpIsfgv?mode=page\&blockIds=xfKeXF0XbEKUVBNHCw3XL)). 

##### Drone Derived LAI Predictors

We extracted spectral and canopy structure information for each field plot from the corresponding drone orthomosaics and DSMs. The extracted data is then transformed to be appropriate for relating to the field data in predicting LAI. These data are stored in two forms (1) small image "chips" with extracted data for each field plot as well as (2) continuous mosaics of predictors that proved to be important in the statistical upscaling model. The "chips" are at the original resolution of the data, and the continuous predictors are at 1m resolution. Both datasets are stored in the project OneDrive folder in COG format. These data total 220 MB in size.

##### Predicted LAI for Drone Mosaics

The statistical model produces 1m resolution LAI predictions for each drone mosaic. The current workflow only produces a "best estimate" (posterior mean), but future versions might produce a large ensemble of predictions to propagate uncertainty. This data is stored as single-layer raster data in GeoTIFF format for each site visit and drone flight (32 in total).

##### Upscaled LAI for a sample of HLS Pixels

Our aggregation and filtering procedure produces upscaled LAI estimates for a collection of \~5000 30 x 30m HLS pixels that are covered by the drone data. These data are sparse in space and time (several dozen to several hundred measurements for a given site and 4 intervals per site), so they are stored as a single table in CSV format in our shared OneDrive folder. 

### Notes on Data Storage and Computing

Because we aggregate the drone datasets to 1m resolution (from their native 3cm resolution), data storage needs for this analysis step are modest, and we store all intermediate datasets in the project OneDrive folder (\~500 MB). All scripts are run locally on a scientific workstation, and the most computationally intensive step is preparing the 1m -resolution drone-based LAI predictors from the original 3cm resolution drone imagery. This step takes several hours to complete on a desktop workstation, but only needs to be done once. 
