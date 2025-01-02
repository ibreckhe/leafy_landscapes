# LL Step 0: Primary Drone Data Processing

Note: You can see how this step fits into the bigger picture of the analysis here: [Leafy Landscapes LAI Workflow](file:///workspace/e094b233-4dce-4599-bc85-a12ab329bda3/7JKZ-alrJKeD8D00oybKY?mode=edgeless\&blockIds=-afNE3aUeaYckGgi83Txl)

### Overview

We collect drone imagery of each site just prior to field data collection and biomass harvest. For each drone flight, the camera on the drone takes several hundred images with the RGB camera. We then process these images into continuous maps.

### Primary Datasets

##### Raw Drone Images

These are stored in JPEG format with embedded georeferencing information. Each flight produces 3-6 GB of raw RGB imagery. Note that we also collect near-IR camera on the drone, but we don't use that information in this project so-far. The original raw drone images are stored in RMBL's Network Attached Storage (NAS) system.

##### Ground Control Points

At each drone survey site, we have a set of visually distinctive features (usually bucket lids stapled to the soil surface), where we have precise geolocations. This allows us to align the mosaiced imagery "true to the world" and ensure that mosaics from multiple flights are aligned with each other. These data are stored in a single CSV formatted table on the project OneDrive folder.

### Analysis Steps

##### Photogrammetry Processing

Our workflow in Agisoft Metashape is semi-automated, with .xml batch scripts recording the major processing steps. A few steps (QC, georeferencing, data export) are stubbornly resistant to automation, so the workflow is mainly manual point-and-click within Metashape.

### Derived Datasets

##### Drone Image Mosaics and DSMs

We use photogrammetry software (Agisoft Metashape) to transform individual images into georeferenced orthomosaics (continuous imagery from a top-down perspective). This procedure also produces 3D structure information in the form of a Digital Surface Model (DSM). Each drone orthomosaic is a single 3-band raster dataset stored in cloud-optimized GeoTIFF (COG) format. The DSM is stored as a separate COG file. Each orthomosaic is 1-2 GB in size, and we have 32 orthomosaics (4 for each of 8 sites).

### Notes on Data Storage and Computing

The total volume of Drone Image Mosaics and DSMs used for this project is about 80 GB. We are storing this data in the project OneDrive folder to facilitate syncing and backup. The drone data processing requires a license to commercial software and a powerful desktop workstation with discrete GPU.
