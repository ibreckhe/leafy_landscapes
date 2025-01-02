# LL Step 1: Assemble and QC field and laboratory data.

Note: You can see how this step fits into the bigger picture of the analysis here: [Leafy Landscapes LAI Workflow](file:///workspace/e094b233-4dce-4599-bc85-a12ab329bda3/7JKZ-alrJKeD8D00oybKY?mode=edgeless\&blockIds=RD9MTKMnSzQttpUQQzZAM)

## Overview

Field-based measurements of Leaf Area Index form "ground truth" of the project. As part of field data collection and laboratory processing, the field team collected extensive plant measurements and harvested all aboveground biomass at several hundred field plots located at 8 sites in the RMBL vicinity. We used these field and lab measurements to reliably measure total LAI, aboveground biomass, leaf water content, and the composition of plant functional types (PFTs).

## Dataset Descriptions

Here are short descriptions of the primary datasets (e.g. "raw" data) that are used by this analysis step.

### Primary Datasets

##### Field Measurements of Plant Cover and Composition

For each field plot, we collected gross plant cover and height data by functional type (shrubs, forbs, grasses). These were recorded on paper datasheets and then entered into to the master lab and field dataset stored in Google Sheets. These data are relatively small (<50 KB) and are stored in a few columns in the master field and lab data table. We use Google Sheets for this data becuase it allows real-time syncing between devices and tracks data versions automatically.

##### Laboratory Measurements of Plant Traits and Biomass

For each plot, we harvested all aboveground biomass for each functional type and brought it back to the lab for processing. After an elaborate lab processing protocol, we entered plant trait information into the master lab and field dataset on Google Sheets. 

##### Leaf Scans for Each Plot

For each sample from each plot, we select 10-50 representative leaves, arrange them on a flatbed scanner, and take high-resolution RGB scans of each. This results in a series of scan images for each plot. For our 2024 work, we collected 832 images stored in JPEG format totalling about 540 MB.

### Analysis Steps

##### Classify and QC Leaf Scans

We have an automated procedure for classifying leaf scans into "leaf" and "non-leaf" components, stored as an R script (Step1a\_Calc\_Leaf\_Areas.R). This is based on a set of color and brightness thresholds for images after they are converted from RGB to HSI color space. The resulting output images are visually examined, and the thresholds are adjusted for images that show poor classification accuracy with the default thresholds. Using information about the resolution of the scanner, we can then estimate the area of leaves in each scan.

##### Mass-based LAI Upscaling for Each Field Plot

With leaf area information for the subset of leaves for each PFT for each plot, we can then use the dry masses from our laboratory procedure to estimate the leaf area of whole plots and for each plant functional type. This procedure is automated as an R script (Step1b\_Upscale\_Field\_LAI.R).

##### QC Field LAI Estimates

We compare the field-based information on plant cover and height for each PFT to the estimated LAI for each field plot. This allows us to flag potential problems with particular lab measurements. After QC, our work provides field-based LAI estimates for each of the 209 field plots sampled in 2024. This analysis step in included in the R script described above.

### Derived Datasets

##### Scan Leaf Areas

Our leaf scan image classification procedure allows us to measure the total area of plant leaves in each scan image. These data are then combined for each scan to develop measurements of the area of representative leaves for each PFT.

### Notes on Data Storage and Computing

The overall data volume is modest (\~600 MB including leaf scans), but these are primary field data and cannot be recreated. We store the primary field measurements in Google Sheets for data entry and QC, after which we create a "master" finalized data table for downstream analysis which is stored in our RMBL OneDrive folder for the project. All computing that integrates the field datasets is done locally on a single machine for these analysis steps. 
