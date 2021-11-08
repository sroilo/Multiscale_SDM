###################################################
#
# Multiscale_SDM_functions.R
# Purpose:    collection of functions used for running the multiscale SDM models described in the following manuscript:
# Reference:  Roilo et al. "Estimating the efficacy of agri-environment measures for farmland birds across spatial scales".
# Author:     Stephanie Roilo, Technische Universität Dresden
# Date:       last modified on October 8th, 2021
#
###################################################

### var.AICc --------------------------------------------------------
# this function computes the AICc for a univariate GLM model (response ~ predictor)
# package MuMIn is required to run this function
var.AICc <- function (predictor, response) {
  MuMIn::AICc(glm(response ~ predictor, family = "binomial"))
}

### selectVars -----------------------------------------------------
# function to select best set of uncorrelated variables, given a threshold for correlation 
# If variables are calculated at different spatial scales (e.g. Arable_200, Arable_500, ..), only the best-fitting spatial scale is retained 
# requires packages MuMIn and stringr
selectVars <- function(pred_names, response_name, data, cor_mat=NULL, threshold) {
  
  # Function for calculating AICc - we use univariate GLMs with linear and quadratic terms
  var.AICc <- function (predictor, response) {
    MuMIn::AICc(glm(response ~ predictor, family = "binomial"))
  }
  # Calculate AICc for all predictor variables
  aic_imp <- apply(data[pred_names], 2, var.AICc, response= data[,response_name])
  # sort from lowest AICc to highest
  aic_imp = sort(aic_imp)
  # remove the same variables which are calculated at different spatial scales 
  # WARNING: this only works if the variable names are different in their first 4 letters!
  AICc.init = stringr::str_sub(names(aic_imp), 1,4)
  dups = duplicated(AICc.init)
  aic_imp = aic_imp[which(dups==F)]
  # Names of sorted variables
  sort_imp <- names(aic_imp)
  # Calculate correlation matrix if not provided in function call
  if (is.null(cor_mat)) {
    cor_mat <- cor(data[sort_imp], method='spearman')
  }
  # Identify correlated variable pairs
  diag(cor_mat)=NA
  pairs <- which(abs(cor_mat)>= threshold, arr.ind=T) 
  # Identify which variables should be excluded
  exclude <- NULL
  for (i in 1:length(sort_imp))
  {
    if ((sort_imp[i] %in% row.names(pairs))&
        ((sort_imp[i] %in% exclude)==F)) {
      cv <- cor_mat[setdiff(row.names(cor_mat),exclude),sort_imp[i]]
      cv <- cv[setdiff(names(cv),sort_imp[1:i])]
      exclude <- c(exclude,names(which((abs(cv)>=threshold)))) 
    }
  }
  # Select set of uncorrelated predictors and return it
  pred_sel <- sort_imp[!(sort_imp %in% exclude)]
  return(pred_sel)
}

### filterByProximity --------------------------------------------
# function to filter points by proximity, discarding the ones that are closer than a threshold distance 
# modified from https://gis.stackexchange.com/questions/191412/proximity-filter 
filterByProximity <- function(xy, dist) {
  #xy has to be an object of class sf, sfc or sfg
  #dist is the minimum distance wanted between points
  # WARNING: this function is order-sensitive, i.e. it discards all points which fall within the minimum set distance dist from the first one in the order
  # dist has to be given in mapUnits (if using EPSG:3035, in meters)
  d <- st_distance(xy)
  diag(d) <- NA
  d = drop_units(d)
  close <- (d <= dist)
  diag(close) <- NA
  closePts <- which(close,arr.ind=T)
  discard <- matrix(nrow=2,ncol=2)
  if (nrow(closePts) > 0) {
    while (nrow(closePts) > 0) {
      if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
        discard <- rbind(discard, closePts[1,])
        closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
      }
    }
    discard <- discard[complete.cases(discard),]
    if (inherits(discard, "matrix")) {
      return(xy[-discard[, 1], ])
    } else {
      return(xy[-discard[1], ])
    }
  }
  if (nrow(closePts) == 0) {
    return(xy)
  }
}

###computeLUV ----------------------------------------------
# function that computes the proportion of different land cover/land use within circular windows of a given radius based on 
# Integrated Administration and Control System (IACS) data and raster layers
#  Input data are: 
# pts = point shapefile of class sf (locations for which you want to compute the land-use variables), 
# radius = radius in meters of the circular window that you want to use around each point
# lcstack = stack of binary (0/1 values) rasters of different land cover types (e.g. grassland, forest, urban, water, etc.)
# invekosly = shapefile layer with field-level information on application of agri-environment schemes and organic farming
# efaly = shapefile layer with field-level information on Ecological Focus Areas application on or within the field
# output is the sf point object with additional attribute columns for each land cover/use variable calculated
computeLUV = function(pts, radius, lcstack, invekosly, efaly) {
  pb <- txtProgressBar(min = 0, max = nrow(pts), style = 3) # produces a progress bar 
  print("Please be patient while I extract the values from your rasters, then the progress bar will start")
  ncols_input = ncol(pts)
  nly_input = terra::nlyr(lcstack)
  # create a buffer of radius=radius around each point
  buf = st_buffer(pts, dist = radius)
  # extract the sum of number of pixels of each land cover class  
  ext = as.data.frame(terra::extract(lcstack, vect(buf), fun = "sum", na.rm=T, method="simple"))
  # transform from cell numbers to area in m^2 (xres * yres to get the area of one cell) and divide by total area of buffer
  ext = ext[,2:ncol(ext)]*xres(lcstack)*yres(lcstack)/(radius*radius*pi)
  names(ext) = paste(names(ext), radius, sep="_")
  # attach to the pts attribute table
  pts = dplyr::bind_cols(pts, round(ext, digits=3))
  
  # compute intersection between buffers and IACS shapefile data
  int = st_intersection(x = buf[,1],  y = st_make_valid(invekosly))
  int = int[,c("n", "crop_code", "Field_type", "EFA_TYP","Organic", "AES_INFO1", "AES_INFO2", "geometry")]
  int$area = set_units(st_area(int), NULL)
  #compute intersection between points and EFA layer
  efaint = st_intersection(x = buf[,1],  y = st_make_valid(efaly)) 
  efaint$area = set_units(st_area(efaint), NULL)
  
  # loop through each observation point and compute the land-use variables
  for (i in seq_along(pts$n)) {
    setTxtProgressBar(pb, i) 
    
    subset = st_drop_geometry(int[int$n == pts$n[i],])
    subefa = st_drop_geometry(efaint[efaint$n == pts$n[i],])
    if (nrow(subset) == 0) {
      pts[i,c("OKO", "Arable", "CovCrop","GrM","Buffer", "Fallow", "CropSDI")] <- 0 
    } else {
      # compute area of organic farming within the circular window
      pts$OKO[i] = sum(subset[which(subset$Organic == "T"), "area"])
      # compute area of arable land within the circular window
      pts$Arable[i] = sum(subset[subset$Field_type == "AL", "area"])  # AL = arable land
      # compute area of cover crops
      pts$CovCrop[i] = sum(subset[which(subset$EFA_TYP == "052" | subset$AES_INFO1=="AL4" | subset$AES_INFO2=="AL4"), "area"])
      # compute area of permanent grassland maintenance
      pts$GrM[i] = sum(subset[(grepl(pattern="GL1|GL2|GL4|GL5", x = subset$AES_INFO1)| grepl(pattern="GL1|GL2|GL4|GL5", x = subset$AES_INFO2)), "area"])
      # compute area of buffer strips
      pts$Buffer[i] = sum(subset[(grepl(pattern="AL1|AL5c|AL5d", subset$AES_INFO1)| grepl(pattern="AL1|AL5c|AL5d", subset$AES_INFO2)| grepl(pattern="065|066", subset$EFA_TYP)), "area"]) +
        sum(subefa[grepl("054|056|057|058|078", subefa$EFA_TYP), "area"])
      # compute area of fallow land
      pts$Fallow[i] = sum(subset[(grepl(pattern="AL5a|AL5b|GL3", subset$AES_INFO1)| grepl(pattern="AL5a|AL5b|GL3", subset$AES_INFO2)| grepl("062", subset$EFA_TYP)), "area"])
      # compute Shannon diversity index for agricultural land use diversity 
      cropdf = aggregate(area ~ crop_code, data= subset, FUN = sum)
      cropdf$prop = cropdf$area/((radius^2)*pi)
      pts$CropSDI[i] = -sum(sapply(X = cropdf$prop, FUN = function(x) {x * log(x)}, simplify = "vector"))
    }
  }
  # divide the areas by the total buffer area, to get proportions (from 0 to 1) and round to 3 digits
  pts[,c("OKO","Arable", "CovCrop","GrM", "Buffer", "Fallow")] = round(
    st_drop_geometry(pts[,c("OKO","Arable", "CovCrop","GrM", "Buffer", "Fallow")])/(radius*radius*pi), digits=3)
  pts[,c("CropSDI")] = round(st_drop_geometry(pts[,c("CropSDI")]), digits=3)
  
  #rename columns to add the information on the radius (in m)
  names(pts) = c( names(pts)[1:sum(ncols_input, nly_input)], paste(c("OKO", "Arable", "CovCrop","GrM", "Buffer", 
                                                                     "Fallow", "CropSDI"), radius, sep="_") )
  
  close(pb)
  return(pts)
}


