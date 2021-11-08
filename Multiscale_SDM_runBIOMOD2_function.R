###################################################
#
# Multiscale_SDM_runBIOMOD2_function.R
# Purpose:   function that fits 5 model algorithms (GLM, GAM, RF, GBM and MAXENT) 
#            and the ensemble model (weighted mean by AUC with threshold of AUC=0.7 to discard
#            uninformative models). The function returns variable importance values and evaluation 
#            scores (AUC, TSS and KAPPA, sensitivity and specificity based on AUC thresholds).
#            The function also projects the model to the current land-use conditions (CURR scenario), and to 
#            a new scenario in which agri-environment measures are set to zero (NOAE scenario).
# Reference: Roilo et al. "Estimating the efficacy of agri-environment measures for farmland birds across spatial scales".
# Author:    Stephanie Roilo, Technische Universität Dresden
# Date:      last modified on October 8th, 2021
#
###################################################

# input data needed:
#  vars = dataframe with one "presence" column (values either 1 or 0), 
#         and other columns with values of the explanatory variables used for modelling.
# species = Latin name of the modelled species (e.g. "Alauda arvensis").
# sel.vars = vector with names of the selected explanatory variables to use in the model.
# nruns = number of model runs 
# maxent.path = path to the maxent.jar file

#### runBIOMOD2 function -------------------------------------------
runBIOMOD2 <- function(vars, species, sel.vars, nruns, maxent.path) {
  # format input data
  input.data = BIOMOD_FormatingData(resp.var = vars[,"presence"],
                                    expl.var = vars[, sel.vars],
                                    resp.xy = vars[,c("X", "Y")],
                                    resp.name = species,
                                    na.rm = T)
  # modelling options
  options = BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = maxent.path, 
                                                          beta_lqp = -1.0,  ##automatic setting
                                                          linear= TRUE, quadratic = TRUE, 
                                                          product = FALSE,
                                                          threshold=TRUE, hinge=FALSE,  
                                                          maximumiterations= 200),
                                   RF =list(ntree = 500,  # no. of trees
                                            nodesize = 5, # Minimum size of terminal nodes (default is 5 for regression, 1 for classification)
                                            mtry =  1), # mtry = Number of variables randomly sampled as candidates at each split.  
                                   GBM = list(n.trees = 100,
                                              shrinkage = 0.1, #'learning rate
                                              bag.fraction = 0.5,
                                              interaction.depth= 1), # tree complexity
                                   GLM = list( # type = 'quadratic', 
                                     interaction.level = 0,
                                     myFormula = NULL,
                                     test = 'AIC',
                                     family = 'binomial',
                                     #mustart = 0.5,
                                     control = glm.control(epsilon = 1e-08, maxit = 100, trace = FALSE)),
                                   GAM = list(algo='GAM_mgcv',
                                              type='s_smoother',
                                              interaction.level= 0,
                                              myFormula = NULL,
                                              k= 5,
                                              test= 'AIC',
                                              family = 'binomial',
                                              control = glm.control(epsilon = 1e-08, maxit = 100, trace = FALSE)),
  )
  # run models
  output <- BIOMOD_Modeling(input.data,
                            models = c("GLM", "GAM", "RF", "GBM", "MAXENT.Phillips"), 
                            models.options = options,
                            NbRunEval = nruns,  # nr. of runs for each model type
                            DataSplit = 70, # how to split data for calibrating and testing.
                            Prevalence = 0.5, # equal weight for presence and absence locations
                            VarImport = 3, # nr. of permutations used to calculate variable importance
                            models.eval.meth = c("ROC", "KAPPA", "TSS"),  # methods used to evaluate models
                            SaveObj = FALSE,
                            rescal.all.models = FALSE, # scale model predictions?
                            do.full.models = FALSE,
                            modeling.id = ID)
  
  # get the evaluation scores
  eval = get_evaluations(output)
  eval_auc = as.data.frame(eval["ROC","Testing.data",,,])
  eval_auc$mean = round(rowMeans(eval_auc[,1:nruns], na.rm=T),digits=3)
  eval_auc$sd = round(apply(eval_auc[1:nrow(eval_auc),1:nruns],1,sd, na.rm=T), digits=3)
  eval_kappa = as.data.frame(eval["KAPPA","Testing.data",,,])
  eval_kappa$mean = round(rowMeans(eval_kappa[,1:nruns], na.rm=T),digits=3)
  eval_kappa$sd = round(apply(eval_kappa[1:nrow(eval_kappa),1:nruns],1,sd, na.rm=T), digits=3)
  eval_tss = as.data.frame(eval["TSS","Testing.data",,,])
  eval_tss$mean = round(rowMeans(eval_tss[,1:nruns], na.rm=T),digits=3)
  eval_tss$sd = round(apply(eval_tss[1:nrow(eval_tss),1:nruns],1,sd, na.rm=T), digits=3)
  eval_sensitivity = as.data.frame(eval["ROC","Sensitivity",,,])
  eval_sensitivity$mean = round(rowMeans(eval_sensitivity[,1:nruns], na.rm=T),digits=2)
  eval_sensitivity$sd = round(apply(eval_sensitivity[1:nrow(eval_sensitivity),1:nruns],1,sd, na.rm=T), digits=2)
  eval_specificity = as.data.frame(eval["ROC","Specificity",,,])
  eval_specificity$mean = round(rowMeans(eval_specificity[,1:nruns], na.rm=T),digits=2)
  eval_specificity$sd = round(apply(eval_specificity[1:nrow(eval_specificity),1:nruns],1,sd, na.rm=T), digits=2)
  eval_all = cbind(AUC_mean = eval_auc$mean, AUC_sd = eval_auc$sd,
                   KAPPA_mean = eval_kappa$mean, KAPPA_sd = eval_kappa$sd,
                   TSS_mean = eval_tss$mean, TSS_sd = eval_tss$sd,
                   Sensitivity_mean = eval_sensitivity$mean, Sensitivity_sd = eval_sensitivity$sd,
                   Specificity_mean = eval_specificity$mean, Specificity_sd = eval_specificity$sd)
  row.names(eval_all) = row.names(eval_auc)
  # write evaluation metrics to file
  write.table(eval_auc, file=paste0(eval.dir, species, "_AUC.csv"), sep=";", dec=",", row.names=T)
  write.table(eval_kappa, file=paste0(eval.dir, species, "_KAPPA.csv"), sep=";", dec=",", row.names=T)
  write.table(eval_tss, file=paste0(eval.dir, species, "_TSS.csv"), sep=";", dec=",", row.names=T)
  write.table(eval_sensitivity, file=paste0(eval.dir, species, "_Sensitivity.csv"), sep=";", dec=",", row.names=T)
  write.table(eval_specificity, file=paste0(eval.dir, species, "_Specificity.csv"), sep=";", dec=",", row.names=T)
  write.table(eval_all, file=paste0(eval.dir, species, "_all.csv"), sep=";", dec=",", row.names=T)
  
  # extract variable importance scores
  mvi = get_variables_importance(output)
  models = row.names(eval_all)
  vimp_all = data.frame(row.names = row.names(mvi[,1,,]))
  # calculate mean relative variable importance of explanatory variables and save them
  for(jj in 1:length(models)){
    varimp.model <- mvi[,models[jj],,]
    # calculate relative importance scores
    rel.varimp <- lapply(1:ncol(varimp.model), function(x) varimp.model[,x]/sum(varimp.model[,x]))
    rel.varimp <- do.call(cbind, rel.varimp)
    # calculate mean and SD 
    mean.rel.varimp <- sapply(1:nrow(rel.varimp), function(x) mean(rel.varimp[x,], na.rm=TRUE))
    sd.rel.varimp <- sapply(1:nrow(rel.varimp), function(x) sd(rel.varimp[x,], na.rm=TRUE))
    mean.rel.varimp <- data.frame(variable=rownames(rel.varimp), mean_importance=round(mean.rel.varimp,3), sd=round(sd.rel.varimp,3))
    write.table(mean.rel.varimp, file=paste0(varimp.dir, species, "_", models[jj], ".txt"))
    # put the values from all algorithms in the same dataframe
    vimp_all = cbind(vimp_all, mean.rel.varimp[,c("mean_importance", "sd")])
  }
  vimp_all = rbind(rep(models, each=2), vimp_all)
  write.table(vimp_all, file=paste0(varimp.dir, species, "_all.csv"), sep=";", dec=",", row.names=T)
  
  # calculating ensemble models
  Ensemble = BIOMOD_EnsembleModeling(modeling.output = output,
                                     chosen.models = "all",
                                     em.by = 'PA_dataset+repet',  # ensemble models are evaluated on the same part of the data than formal-models they are made with. In this case the way ensemble models are evaluated are fair compared to formal models evaluation.
                                     eval.metric = c('ROC'),
                                     eval.metric.quality.threshold = c(0.7), # AUC threshold value to select the single-run model for the ensemble
                                     prob.mean = F,
                                     committee.averaging = F,
                                     prob.mean.weight = T,
                                     prob.mean.weight.decay = 'proportional',
                                     VarImport = 3 )
  
  # get ensemble model evaluations
  EMeval <- get_evaluations(Ensemble)
  EMmodels = names(EMeval)
  AUC = numeric()
  KAPPA = numeric()
  TSS = numeric()
  sensitivity=numeric()
  specificity=numeric()
  for(jj in 1:length(EMeval)){
    AUC[jj] <- EMeval[[jj]]["ROC","Testing.data"]
    KAPPA[jj] <- EMeval[[jj]]["KAPPA","Testing.data"]
    TSS[jj] <- EMeval[[jj]]["TSS","Testing.data"]
    sensitivity[jj] <- EMeval[[jj]]["ROC","Sensitivity"]
    specificity[jj] <- EMeval[[jj]]["ROC","Specificity"]   
  }
  dfEM = data.frame(EMmodels, AUC, KAPPA, TSS, sensitivity, specificity)
  dfEM = as.data.frame(t(dfEM[,-1]))
  dfEM[,"mean"] = round(rowMeans(dfEM[,1:length(EMeval)], na.rm=T), 3)
  dfEM[,"sd"] = round(apply(dfEM[,1:length(EMeval)], 1, sd, na.rm=T), 3)
  write.table(dfEM, file=paste0(eval.dir, species, "_EMwmean.csv"), sep=";", dec=",", row.names=T)
  
  # get ensemble model variable importance
  EMvarimp <- get_variables_importance(Ensemble)
  dfEMvarimp = data.frame(variable = rownames(EMvarimp))
  for(jj in 1:length(EMmodels)){
    EMvarimp.model <- EMvarimp[,,EMmodels[jj]]
    EMrel.varimp <- lapply(1:ncol(EMvarimp.model), function(x) EMvarimp.model[,x]/sum(EMvarimp.model[,x]))
    EMrel.varimp <- do.call(cbind, EMrel.varimp)
    EMmean.rel.varimp <- sapply(1:nrow(EMrel.varimp), function(x) mean(EMrel.varimp[x,], na.rm=TRUE))
    dfEMvarimp[, jj +1] <- round(EMmean.rel.varimp,3)
  }
  dfEMvarimp$mean_imp = round(rowMeans(dfEMvarimp[,2:(length(EMmodels)+1)]), 3)
  dfEMvarimp$sd = round(apply(dfEMvarimp[,2:(length(EMmodels)+1)], 1, sd), 3)
  write.table(dfEMvarimp, file=paste(varimp.dir, species, "_EMwmean.csv"), sep=";", dec=",", row.names=F)
  
  ### MAKE PROJECTIONS 
  # Project the model to the current land-use conditions (CURR scenario):
  # load the environmental rasters of the selected variables in a stack
  envdir <- "projection rasters/"   # path to directory containing environmental rasters
  envstack <- stack(paste(envdir, sel.vars, ".tif", sep=""))
  # subselect only the full models from each of the 5 algorithms and project to current scenario
  fullmod = output@models.computed[grep("AllData_Full", output@models.computed)]
  projection <- BIOMOD_Projection(modeling.output = output,
                                  new.env = envstack,
                                  proj.name = ID,
                                  selected.models = fullmod, 
                                  binary.meth = NULL,
                                  filtered.meth = NULL,
                                  compress = "xz",
                                  build.clamping.mask = FALSE)
  # load projection maps
  proj.dict <- paste(gsub(" ", ".", species), "/proj_", ID, "/", sep="")
  proj.files <- list.files(proj.dict, ".grd")
  proj.raster <- stack(paste(proj.dict, proj.files, sep=""))
  # load mean AUC scores for each algorithm and set those with AUC < 0.7 to zero, use the rest as weights
  weights = eval_auc[, "mean"]
  weights[which(weights < 0.7)] <- 0
  ## produce an ensemble model projection by weighted average (weights = mean AUC of each algorithm) and save to file
  EMproj = raster::weighted.mean(proj.raster, w = weights, na.rm = F, filename = paste0(map.dir, species, " current.tif") )
  
  # Project the model to a scenario with no AES/EFA (NOAE scenario):
  # set all AES and EFA raster values to zero
  AES = c("Buffer_200","Buffer_500", "Buffer_1000", "CovCrop_200", "CovCrop_500", "CovCrop_1000",
          "GrM_200", "GrM_500", "GrM_1000", "Fallow_200", "Fallow_500", "Fallow_1000",
          "OKO_200","OKO_500","OKO_1000")
  for (i in c(1:nlayers(envstack))) {
    if (names(envstack[[i]]) %in% AES) {values(envstack[[i]]) <- 0}
  }
  # project the full models from each of the 5 algorithms to the NOAE scenario:
  projection <- BIOMOD_Projection(modeling.output = output,
                                  new.env = envstack,
                                  proj.name = paste0(ID, "_noAES"),
                                  selected.models = fullmod, 
                                  binary.meth = NULL,
                                  filtered.meth = NULL,
                                  compress = "xz",
                                  build.clamping.mask = FALSE)
  # load projection maps
  proj.dict <- paste(gsub(" ", ".", species), "/proj_", ID, "_noAES", "/", sep="")
  proj.files <- list.files(proj.dict, ".grd")
  proj.raster <- stack(paste(proj.dict, proj.files, sep=""))
  # weights remain the same as before
  ## produce an ensemble model projection by weighted average (weights = mean AUC of each algorithm) and save to file
  EMproj = raster::weighted.mean(proj.raster, w = weights, na.rm = F, filename = paste0(map.dir, species, " no AES.tif") )
}




