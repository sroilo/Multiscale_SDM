Description
==============

**Roilo, S., Engler, J. O., Václavík, T., & Cord, A. F. (2022). Landscape‐level heterogeneity of agri‐environment measures improves habitat suitability for farmland birds. Ecological Applications, e2720.
  Please acknowledge this paper when using the code or data outputs.


Multiscale_SDM.R
--------------
This R script is provided to facilitate the implementation of the multiscale SDM modelling framework described in Roilo et al. (in review).
Following inputs are required to run the SDMs:
  - bird dataset (shapefile holding georeferenced point data of different bird species);
  - rasters of environmental variables, e.g. elevation, slope, distance from highways, distance from forest edges, land cover maps for grassland, small woody features and urban areas;
  - IACS data (shapefile holding spatially-explicit information on agricultural fields, e.g. grown crops, applied Agri-Environment Measures, etc.).
The code calls other functions saved in the R scripts "Multiscale_SDM_functions.R", and "Multiscale_SDM_runBIOMOD2_function.R"


Multiscale_SDM_functions.R
--------------
This R script contains useful functions for data filtering, variable selection, and calculation of environmental variables within circular windows 
of varying sizes.


Multiscale_SDM_runBIOMOD2_function.R
--------------
This R script is used for running the multiscale ensemble SDM models.


Corresponding author: 
--------------
Stephanie Roilo **stephanie.roilo@tu-dresden.de**

Please get in touch if you have questions about the code or data.

For more information about the BESTMAP project, please see: www.bestmap.eu 


