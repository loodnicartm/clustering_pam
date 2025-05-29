# clustering_pam
This R script applies Partitioning Around Medoids (PAM) clustering to BSCI time series data in order to identify homogeneous vegetation zones based on their seasonal and interannual dynamics. The approach is tailored for ecological and land surface studies aiming to segment areas by vegetation behavior rather than fixed land cover classes.

-Key Features:
Load preprocessed BSCI time series data for a target region.

Apply temporal smoothing and standardization to prepare data for clustering.

Use PAM clustering to group pixels based on their BSCI temporal profiles.

Visualize spatial distribution of clusters and inspect cluster medoid curves.

Export clustering results for further analysis or integration in GIS workflows.

-Technologies & Libraries Used:
raster, terra, sf – for spatial data management

cluster, factoextra – for PAM clustering and evaluation

ggplot2, reshape2 – for visualization of temporal cluster profiles

-Applications:
Useful for:

Ecosystem classification based on functional behavior

Drought sensitivity analysis

Change detection preprocessing

Land surface phenology studies
