###################################################
#
# Title: Multiscale_SDM.R
# Purpose: fit multiscale ensemble models to predict breeding habitat suitability for farmland birds in the Mulde River Basin, Germany. 
# Reference: Roilo et al. "Estimating the efficacy of agri-environment measures for farmland birds across spatial scales".
# Author:    Stephanie Roilo, Technische Universität Dresden
# Date:      last modified on October 8th, 2021
#
###################################################

# load packages
library(rgdal)
library(raster)
library(sf)
library(dplyr)
library(units)
library(terra)
library(MuMIn)
library(stringr)
library(mapview)
library(corrplot)

# load functions 
source("R scripts/Multiscale_SDM_functions.R")

### BIRD  DATASET PREPARATION  ----------------------------------
# read the shapefile of the Case Study region and create an inner 1 km buffer
mulde = st_read("Mulde/Mulde_EPSG3035.shp") %>% st_buffer(dist= -1000, singleSide=T)
# read bird dataset (observation points of bird species), already cropped to the study region extent 
# (minus 1 km buffer to ensure correct calculation of all env. variables), and
# filtered to retain only observations of possible, probable or certain reproduction
birds = st_read("Birds_2016-2019_filtered.shp") 

### VARIABLES CALCULATION -------------
pts = birds[,c( "X","Y","species","date","year","source","behaviour","unit","reproduction","geometry")]
# extract values from the topographic layers, and distance to highways and to forest edges
topo = rast(c("Model tests/Projection rasters/Elevation.tif",
              "Model tests/Projection rasters/Slope.tif"))
rfdist = rast(c("Model tests/input rasters/Road_dist_5m.tif",
                "Model tests/input rasters/Forest_dist_5m.tif"))
pts = cbind(pts, terra::extract(y = st_drop_geometry(pts[,c("X", "Y")]), x = topo, method = "simple"))
pts = cbind(pts, terra::extract(y = st_drop_geometry(pts[,c("X", "Y")]), x = rfdist, method = "simple"))
# rename columns
names(pts)[10:14] = c("Elevation", "Slope", "D_road","D_forest", "geometry")
# add a unique identifier for each observation in the dataset
pts = cbind(n = 1:nrow(pts), pts)
## compute land use/cover variables 
# load the binary rasters (values 0 or 1) of presence-absence of different land cover types
lcstack = rast(c("Model tests/input rasters/SWF_2015.tif",    # Small Woody Features
                 "Model tests/input rasters/Gras_2015masked.tif",   # Grassland
                 "Model tests/input rasters/UrbanCover_2016masked.tif"))  # urban cover
names(lcstack) = c("SWF", "Grass", "Urban")

## data needs to be subdivided in the different years to compute the LUV (land use variables) from the land-use data of the corresponding year
# 2019 - subselect the data
pts9  = pts[pts$year==2019,]
# load the IACS data for 2019
in9 = st_read("INVEKOS_data/edited/InVeKoS_2019_Sc.shp")
efa9 = st_read("INVEKOS_data/edited/EFA_2019.shp")
# compute the land cover variables and add them to the point attributes table
# set 3 different radius of 200 m, 500 m, and 1 km 
pts9 = computeLUV(pts9, 200, lcstack, in9, efa9)
pts9 = computeLUV(pts9, 500, lcstack, in9, efa9)
pts9 = computeLUV(pts9, 1000, lcstack, in9, efa9)
# do the same for all other years
# 2018 
pts8  = pts[pts$year==2018,]
in8 = st_read("INVEKOS_data/edited/InVeKoS_2018_Sc_NEU.shp")
efa8 = st_read("INVEKOS_data/edited/EFA_2018_NEU.shp")
pts8 = computeLUV(pts8, 200, lcstack, in8, efa8)
pts8 = computeLUV(pts8, 500, lcstack, in8, efa8)
pts8 = computeLUV(pts8, 1000, lcstack, in8, efa8)
## 2017
pts7 = pts[pts$year==2017,]
in7 = st_read("INVEKOS_data/edited/InVeKoS_2017_Sc.shp")
efa7 = st_read("INVEKOS_data/edited/EFA_2017.shp")
pts7 = computeLUV67(pts7, 200, lcstack, in7, efa7)
pts7 = computeLUV67(pts7, 500, lcstack, in7, efa7)
pts7 = computeLUV67(pts7, 1000, lcstack, in7, efa7)
# 2016
pts6 = pts[pts$year==2016,]
in6 = st_read("INVEKOS_data/edited/InVeKoS_2016_Sc.shp")
efa6 = st_read("INVEKOS_data/edited/EFA_2016.shp")
pts6 = computeLUV67(pts6, 200, lcstack, in6, efa6)
pts6 = computeLUV67(pts6, 500, lcstack, in6, efa6)
pts6 = computeLUV67(pts6, 1000, lcstack, in6, efa6)
# put all data from the different years together again:
pts = rbind(pts9, pts8, pts7, pts6)
# save the variables to file (for reproducibility)
vars = st_drop_geometry(pts)
write.table(vars, paste0("Model tests/Full_dataset_2016-19.csv"), sep=";", dec=",", row.names = F)
rm(birds, pts6, pts7, pts8, pts9, rfdist, topo, lcstack, 
   mulde, in9, in8, in7, in6, efa9, efa8, efa7, efa6)

### SPECIES-SPECIFIC MODELS ---------------------------------
# load the table with presence/absence of birds and values of environmental variables calculated earlier
birds = read.table("Model tests/Full_dataset_2016-19.csv", sep=";", dec=",", header=T)
# set the name of folder where results will be stored
ID <- "25Apr_Quail" 
# set the (scientific) name of species to be modelled
mspecies = "Coturnix coturnix"
#create a folder to contain all results
dir.create(paste0("Model tests/", ID), showWarnings=FALSE, recursive=TRUE)
# convert the table to spatial points dataframe
birds = cbind(birds, birds$X, birds$Y)
birds = st_as_sf(birds, coords = c("birds$X","birds$Y"), crs="epsg:3035")
# remove the NAs
birds = na.omit(birds) 
# subselect only observations of the selected species as presence points, and join the other data to the background dataset
pres = birds %>% filter(species == mspecies)  
# filter the presence records by 20 m distance, to remove observations falling in the same environmental raster cell, regardless of the year:
pres = filterByProximity(pres, dist = 20)  
#visual check of the data, and distribution by year
plot(pres[,"year"], pch=20, nbreaks=19)

### select other birds´ observations as background (absence) points
bkg = birds[which(birds$species!= mspecies),]
# remove points that are closer than 500 m to the presence points of the modelled species
bkg = bkg[which(st_intersects(bkg, st_union(st_buffer(pres, dist=500)), sparse=F) == F),]
# filter by proximity to ensure that background points are > 20 m from each other
bkg = filterByProximity(bkg, dist = 20)
# subselect 10 times as many background points as presence points
bkg = bkg[sample(nrow(bkg), size = nrow(pres)*10, set.seed(13), replace=F),]
# put together presence and background points
pts = rbind(pres, bkg)
pts$presence = ifelse(pts$species== mspecies, 1, 0)
# visual check of the data on interactive map
mapview(pts, zcol = "presence")  # looks OK
vars = st_drop_geometry(pts)
vars = cbind(vars[,1:4], presence = vars$presence, vars[,5:38])
# save single-species dataset to file, for reproducibility
write.table(vars, paste0("Model tests/", ID, "/", ID, "_2016-19.csv"), sep=";", dec=",", row.names = F)
rm(pres, bkg)

### VARIABLE SELECTION ---------------------------------------
# fit univariate GLM models and calculate the AICc
AICc.imp <- apply(vars[,6:ncol(vars)], 2, var.AICc, response = vars$presence)
# Sort the predictors; the lower the AICc the better, and save to file
AICc.imp = sort(AICc.imp)
write.table(AICc.imp, paste0("Model tests/", ID, "/AICc.imp.csv"), sep=";", dec=",", row.names=T)
# check for correlation between variables:
cmat <- cor(vars[,6:ncol(vars)], method="spearman") 
corrplot.mixed(cmat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)
# select the best set of uncorrelated variables, e.g. those with lowest AICc in the univariate linear models
sel.vars = selectVars(pred_names = names(vars[,6:ncol(vars)]), response_name = "presence",
                      data = vars, cor_mat = NULL, threshold = 0.7)
rm (cmat, AICc.imp)

### MODELLING ----------------------------------------------------------------------
library(biomod2)
library(rJava)
# load BIOMOD function
source("R scripts/Multiscale_SDM_runBIOMOD2_function.R")
# set work directory where the results will be stored
setwd("Model tests/")

# create additional directories for output files
varimp.dir <- paste(ID, "/varimp/", sep="")
eval.dir <- paste(ID, "/evaluation/", sep="")
map.dir <- paste(ID, "/maps/", sep="")
dir.create(varimp.dir, showWarnings=FALSE, recursive=TRUE)
dir.create(eval.dir, showWarnings=FALSE, recursive=TRUE)
dir.create(map.dir, showWarnings=FALSE, recursive=TRUE)

# set species name and selected variables (if not loaded beforehand)
species = mspecies
# set number of model runs
nruns = 10
# specify the path to the maxent.jar file
maxent.path = "C:\\Users\\sroilo\\Documents\\R\\R-4.0.2\\library\\dismo\\java"

# run the model!
# the function already projects the model onto the current land-use scenario (CURR), as well as onto a scenario in which
# agri-environment measures are set to zero across the entire landscape. The projections are saved to file in the mapdir.
runBIOMOD2(vars, species, sel.vars, nruns, maxent.path)

### VARIABLE RESPONSE PLOTS --------------------------------------------------------
# load the single-algorithm models in order to extract the data used for modelling
specdot = gsub(" ", ".", species)
Smodels <- get(load(paste0(specdot, "/", specdot, ".", ID, ".models.out")))
data = get_formal_data(Smodels, "expl.var")
# load the ensemble models
Allmod = BIOMOD_LoadModels(get(load(paste0(specdot, "/", specdot, ".", ID,"ensemble.models.out"))))
# make response plots for the ensemble models only (not the single algorithm models), and save image in the varimp directory
resp = response.plot2(Allmod, Data=data, show.variables = sel.vars, 
                      plot = T, save.file = "jpeg", name = paste0(varimp.dir, species, "_respPlots"),
                      ImageSize = 800, col = c(1:length(sel.vars)))
# calculate mean and sd of the response curve across the 10 model runs, and save to file
meanVarResp = data.frame(n = 1:100)
for ( i in 1:length(sel.vars)) {
  df = resp[resp$expl.name == sel.vars[i],]
  meanVarResp[,paste(sel.vars[i], "_expl.val", sep="")] = aggregate(pred.val ~ expl.val, data= df, FUN=mean)[1]
  meanVarResp[,paste(sel.vars[i], "_mean", sep="")] = aggregate(pred.val ~ expl.val, data= df, FUN=mean)[2]  
  meanVarResp[,paste(sel.vars[i], "_SD", sep="")] = aggregate(pred.val ~ expl.val, data= df, FUN=sd)[2] 
}
write.table(meanVarResp, paste0(varimp.dir, species,"_meanVarResponse.csv"), dec=",", sep=";", row.names = F)

### PROJECTIONS INTO NEW SCENARIOS -----------------------------------------------
# project the ensemble model into a new land-use scenario with increased uptake levels of agri-environment measures 
# load the environmental rasters of the selected variables in a stack
varnames = c("Elevation","Slope","D_road","D_forest","SWF_200","Grass_200","Urban_200","Arable_200","CovCrop_200","CropSDI_200",
             "SWF_500","Grass_500","Urban_500","Arable_500","CovCrop_500" ,"SWF_1000" , "Grass_1000", "Urban_1000" ,  
             "Arable_1000" , "CovCrop_1000", "CropSDI_1000")
NEWenvstack <- stack(c( paste("projection rasters/", varnames, ".tif", sep=""),  ## these are the changed rasters with increased AEM proportions
                        paste("projection rasters/Scenario3/", 
                              c("OKO_200", "Buffer_200", "Fallow_200", "GrM_200", 
                                "OKO_500", "Buffer_500", "Fallow_500", "GrM_500",
                                "OKO_1000", "Buffer_1000", "Fallow_1000", "GrM_1000"), ".tif", sep="") ) )

# subselect only the full models from each of the 5 algorithms and project to new scenario
fullmod = Smodels@models.computed[grep("AllData_Full", Smodels@models.computed)]
Allmod = BIOMOD_LoadModels(get(load(paste0(specdot, "/", specdot, ".", ID,".models.out"))))
modvars = read.table(paste0(varimp.dir, species, "_EMwmean.csv"), sep=";", dec=",", header=T)[,1]
projection <- BIOMOD_Projection(modeling.output = Smodels,
                                new.env = NEWenvstack[[modvars]],
                                proj.name = paste0(ID, "_moreAES"),
                                selected.models = fullmod, 
                                binary.meth = NULL,
                                filtered.meth = NULL,
                                compress = "xz",
                                build.clamping.mask = FALSE)
# load projection maps
proj.dict <- paste(gsub(" ", ".", species), "/proj_", ID, "_moreAES", "/", sep="")
proj.files <- list.files(proj.dict, ".grd")
proj.raster <- stack(paste(proj.dict, proj.files, sep=""))
# load mean AUC scores for each algorithm and set those with AUC < 0.7 to zero, use the rest as weights
eval_auc = read.table(paste0(eval.dir, species, "_AUC.csv"), sep=";", dec=",", header=T)
weights = eval_auc[, "mean"]
weights[which(weights < 0.7)] <- 0
## produce an ensemble model projection by weighted average (weights = mean AUC of each algorithm) and save to file
EMproj = raster::weighted.mean(proj.raster, w = weights, na.rm = F, filename = paste0(map.dir, species, " more AES.tif") )

rm(list=ls())
