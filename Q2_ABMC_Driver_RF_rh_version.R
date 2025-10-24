library(dplyr)
library(data.table)
library(caret)
library(randomForest)
library(car) 
library(corrplot)
library(pdp)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(stringr)
#LOAD DATA####
rm(list=ls())

#Directories
current_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
data_directory <- (dirname(current_directory))
setwd(data_directory)
rm(current_directory)


data=read.csv(paste0(data_directory,"/Data/CleanedData/FiSL_CompleteData_091524.csv"))
data<-data%>%
  dplyr::select(fire_scar:dec.ind,BUI:CMI_sm.d)
data$abmc<-factor(data$abmc, levels=c("conifer","mixed", "aspen", "birch"))
data<-data%>%
  mutate(prefire.trees=prefire.trees+1)%>% #For three plots with 0 tree C
  mutate(moist=ifelse(moisture_class=="xeric",1,
                      ifelse(moisture_class=="subxeric",2,
                             ifelse(moisture_class=="mesic-subxeric",3,
                                    ifelse(moisture_class=="mesic",4,
                                           ifelse(moisture_class=="mesic-subhygric",5,
                                                  ifelse(moisture_class=="subhygric",6,NA)))))))%>%
  dplyr::select(-moisture_class)
data<-data%>%
  mutate(biomass=ifelse(is.na(biomass),0,biomass))%>%
  mutate(density=ifelse(is.na(density),0,density))
d<-data%>%
  dplyr::select(fire_scar,site,dec.ind)%>%
  na.omit()
d<-d%>%
  group_by(fire_scar,site)%>%
  summarize(mean.dec.ind=mean(dec.ind))

data<-full_join(data,d)
data<-data%>%
  mutate(dec.ind=ifelse(is.na(dec.ind),mean.dec.ind,dec.ind))%>%
  dplyr::select(-mean.dec.ind)
rm(d)

conifer<-read.csv(paste0(data_directory,"/Data/InterimData/JFSPBobydata.csv"))
conifer$plot<-conifer$site
names(conifer)
conifer=conifer%>%
  dplyr::rename(abmc=ABMC.cat,dcm=DCM.cat,dec.ind=DECID_INDEX,prefire.depth=prefire.soil.depth,comb.depth=comb.soil.depth,postfire.depth=postfire.soil.depth,t=temperature,rh=relative_humidity,wspeed=wind_speed,DC=drought_code,DMC=drought_moisture_code,DSR=daily_severity_rank,FWI=fire_weather_index,BUI=buildup_index,ISI=initial_spread_index,FFMC=fine_fuel_moisture_code)%>%
  dplyr::mutate(moist=ifelse(moisture_class=="xeric",1,
                             ifelse(moisture_class=="subxeric",2,
                                    ifelse(moisture_class=="subxeric to mesic",3,
                                           ifelse(moisture_class=="mesic",4,
                                                  ifelse(moisture_class=="submesic",5,
                                                         ifelse(moisture_class=="subhygric",6,NA)))))))%>% #For three plots with 0 tree C
  dplyr::select(fire_scar,site,plot,dcm,abmc,prefire.trees,comb.trees,postfire.trees,comb.trees.prop,prefire.below,comb.below,postfire.below,comb.below.prop,prefire.bgT,comb.bgT,postfire.bgT,comb.bgT.prop,prefire.depth,comb.depth,postfire.depth,elevation,latitude,longitude,ecoregion,slope,aspect,moist,stdage,density,biomass,dec.ind,BUI,DC,DMC,DSR,FFMC,FWI,ISI,rh,t,wspeed)%>%
  mutate(longitude=longitude*-1)
data$site<-as.character(data$site)
conifer_normals30<-read.csv(paste0(data_directory,"/Data/InterimData/JFSPconifer_extractedClimate_normals.csv"))
conifer_normals2004<-read.csv(paste0(data_directory,"/Data/InterimData/JFSPconifer_extractedClimate_2004.csv"))
colnames(conifer_normals2004)<-paste0(colnames(conifer_normals2004),"_2004")
conifer_normals2004$X_2004<-NULL
conifer_normals2004$latitude_2004<-NULL
conifer_normals2004$longitude_2004<-NULL
conifer_normals2004<- conifer_normals2004%>%
  dplyr::rename(site=site_2004)
conifer_normals<-full_join(conifer_normals2004,conifer_normals30)
conifer_normals<-conifer_normals%>%
  mutate(MWMT.d=MWMT_2004-MWMT)%>%
  mutate(MAP.d=MAP_2004-MAP)%>%
  mutate(MSP.d=MSP_2004-MSP)%>%
  mutate(SHM.d=SHM_2004-SHM)%>%
  mutate(CMD.d=CMD_2004-CMD)%>%
  mutate(CMI.d=CMI_2004-CMI)%>%
  mutate(RH.d=RH_2004-RH)%>%
  mutate(CMD_sp.d=CMD_sp_2004-CMD_sp)%>%
  mutate(CMD_sm.d=CMD_sm_2004-CMD_sm)%>%
  mutate(CMI_sp.d=CMI_sp_2004-CMI_sp)%>%
  mutate(CMI_sm.d=CMI_sm_2004-CMI_sm)%>%
  mutate(RH_sp.d=RH_sp_2004-RH_sp)%>%
  mutate(RH_sm.d=RH_sm_2004-RH_sm)%>%
  mutate(Tave_sp.d=Tave_sp_2004-Tave_sp)%>%
  mutate(Tave_sm.d=Tave_sm_2004-Tave_sm)%>%
  mutate(PPT_sp.d=PPT_sp_2004-PPT_sp)%>%
  mutate(PPT_sm.d=PPT_sm_2004-PPT_sm)
conifer_normals<-conifer_normals%>%
  dplyr::select(site,MWMT.d:PPT_sm.d)
conifer<-full_join(conifer,conifer_normals)
rm(conifer_normals2004,conifer_normals30,conifer_normals)

names(data)
names(conifer)

# #Writeout conifer data to get climate numbers
# climatecon<-conifer%>%
#   dplyr::select(site,latitude,longitude)
# write.csv(climatecon,"C:/Users/Test/Documents/MS THESIS/SoilModelRedo/Data/InterimData/JFSPCONIFER_extractionData.csv")

data<-bind_rows(data,conifer)
vegSI<-read.csv(paste0(data_directory,"/Data/CleanedData/ndvi_evi2_slopes_101524.csv"))

data<-full_join(data,vegSI)
data<-filter(data,!is.na(dcm))
rm(conifer,vegSI)
data<-data%>%
  mutate(moist2=ifelse(moist==1|moist==2,1,
                       ifelse(moist==3|moist==4,2,
                              ifelse(moist==5|moist==6,3,NA))))
data$stdage<-as.numeric(data$stdage)
# 
# #Extract mean climate normals
# c=read.csv("C:/Users/Test/Documents/MS THESIS/SoilModelRedo/Data/InterimData/JFSPconifer_extractedClimate_NormalsForRegionMeans.csv")
# c1=read.csv("C:/Users/Test/Documents/MS THESIS/SoilModelRedo/Data/InterimData/JFSPconifer_extractedClimate_NormalMCMTForRegionMeans.csv")
# c=full_join(c,c1)
# 
# names(annual_yearof) #from near end of climate na script: C:\Users\Test\Documents\MS THESIS\Final_RawData\ClimateRawData\ClimateData.R
# d=annual_yearof%>%
#   select(MAT.n:PAS.n)
# names(d) <- sub("\\.n$", "", names(d))
# c=c%>%
#   select(MAT:MAP,PAS:MCMT)
# means=bind_rows(d,c)
# means%>%
#   summarize(MAT=mean(MAT),MWMT=mean(MWMT),MCMT=mean(MCMT),MAP=mean(MAP),PAS=mean(PAS))
# library(plotrix)
# means%>%
#   summarize(MAT=std.error(MAT),MWMT=std.error(MWMT),MCMT=std.error(MCMT),MAP=std.error(MAP),PAS=std.error(PAS))

#START ANALYSIS####

##COLINEARITY####


d<-data%>%
  dplyr::select(comb.bgT,prefire.trees,prefire.below,elevation,latitude,longitude,slope,stdage:evi.yr.prior,moist2)
names(d)

subset<-d%>%
  dplyr::select(comb.bgT,prefire.trees,prefire.below,density,stdage,moist2,evi.yr.prior,CMI_sp.d,CMI_sm.d)

cor <- cor(subset)
# cor=abs(cor)
par(mfrow = c(1, 1)) 
corrplot(cor, method = "number", type="upper")
corrplot(cor, method = "pie", type="upper")
corrplot(cor, method = "shade", type="upper")
corrplot(cor, method = "color", type="upper")
corrplot(cor, method = "circle", type="upper")
corrplot(cor, method = "square", type="upper")

set.seed(123)  # for reproducibility
rf_model <- randomForest(comb.bgT ~ ., data = subset, importance = TRUE)

# View variable importance
importance(rf_model)
varImpPlot(rf_model)

#plot
ggplot(data,aes(x=prefire.below,y=comb.bgT,color=as.factor(moist)))+geom_point(aes(size=1))+theme_minimal()+geom_smooth(method="lm",se=F)+scale_color_viridis_d()
+geom_smooth(method="loess")+
  )

guides(color = guide_legend(override.aes = list(size = 5)))

#Top two ranking vars for climate divergence: CMI_sp.d then CMI_sm.d/MSP.d
#going with evi.yr.prior for spectral
#weather going with wspeed and DC and ISI
#bio/topo: prefire above/below,density,stand age, moisture

##CORRELATION ####
model_data<-data%>%
  dplyr::select(abmc,comb.bgT,prefire.trees,prefire.below,density,stdage,moist,evi.yr.prior,CMI_sp.d,CMI_sm.d,wspeed,DC,ISI)

subset<-filter(model_data,abmc=="conifer")
cor <- cor((subset%>%dplyr::select(-c(abmc))))
# cor=abs(cor)
corrplot(cor, method = "number", type="upper")

#FIT MODEL
subset<-data%>%dplyr::select(c(abmc,comb.bgT,prefire.trees,prefire.below,stdage,moist,CMI_sp.d,CMI_sm.d,wspeed))
subset<-filter(subset,abmc=="conifer")


# 

rf_model <- randomForest(comb.bgT ~ ., data = subset, ntree = 500, mtry = 3, importance = TRUE)
print(rf_model)


rf.best.fit <- rf_model$finalModel
var.imp.dt <- cbind(data.table(var = rownames(rf.best.fit$importance)), data.table(rf.best.fit$importance))[order(IncNodePurity, decreasing = T)]
# 
# 
# pd.list <- list()
# cnt = 1
# 
# vars <- var.imp.dt$var
# var.imp.dt
# 
# # Calculate Mean Squared Error and R-squared of the final model
# mse_final <- mean(rf_model$mse)
# rsq_final <- tail(rf_model$rsq, 1)  # Last value of R-squared after 500 trees
# 
# cat("Mean Squared Error:", mse_final, "\n")
# cat("R-squared:", rsq_final * 100, "% \n")
# 
# par(mfrow = c(3, 3))  # Set up a 3x3 plotting layout for 9 variables
# partialPlot(rf_model, pred.data = subset, x.var = prefire.trees, main = paste("Partial Dependence on", "prefire.trees"))
# partialPlot(rf_model, pred.data = subset, x.var = prefire.below, main = paste("Partial Dependence on", "prefire.below"))
# partialPlot(rf_model, pred.data = subset, x.var = CMI_sp.d, main = paste("Partial Dependence on", "CMI_sp.d"))
# partialPlot(rf_model, pred.data = subset, x.var = CMI_sm.d, main = paste("Partial Dependence on", "CMI_sm.d"))
# partialPlot(rf_model, pred.data = subset, x.var = wspeed, main = paste("Partial Dependence on", "wspeed"))
# partialPlot(rf_model, pred.data = subset, x.var = moist, main = paste("Partial Dependence on", "moist"))
# partialPlot(rf_model, pred.data = subset, x.var = stdage, main = paste("Partial Dependence on", "stdage"))
# partialPlot(rf_model, pred.data = subset, x.var = density, main = paste("Partial Dependence on", "density"))
# partialPlot(rf_model, pred.data = subset, x.var = evi.yr.prior, main = paste("Partial Dependence on", "evi.yr.prior"))
# partialPlot(rf_model, pred.data = subset, x.var = DC, main = paste("Partial Dependence on", "DC"))
# partialPlot(rf_model, pred.data = subset, x.var = ISI, main = paste("Partial Dependence on", "ISI"))
# partialPlot(rf_model, pred.data = subset, x.var = slope, main = paste("Partial Dependence on", "slope"))
# 
# # Set up the plotting layout
# par(mfrow = c(3, 3))  # 3x3 layout for 9 variables
# 
# # Find y-axis limits across all partial plots
# y_limits <- range(
#   partialPlot(rf_model, pred.data = subset, x.var = "prefire.trees", plot = FALSE)$y,
#   partialPlot(rf_model, pred.data = subset, x.var = "prefire.below", plot = FALSE)$y,
#   partialPlot(rf_model, pred.data = subset, x.var = "CMI_sp.d", plot = FALSE)$y,
#   partialPlot(rf_model, pred.data = subset, x.var = "CMI_sm.d", plot = FALSE)$y,
#   partialPlot(rf_model, pred.data = subset, x.var = "wspeed", plot = FALSE)$y,
#   partialPlot(rf_model, pred.data = subset, x.var = "moist", plot = FALSE)$y,
#   partialPlot(rf_model, pred.data = subset, x.var = "stdage", plot = FALSE)$y,
#   partialPlot(rf_model, pred.data = subset, x.var = "density", plot = FALSE)$y,
#   partialPlot(rf_model, pred.data = subset, x.var = "evi.yr.prior", plot = FALSE)$y
# )
# 
# # Generate the partial dependency plots with consistent y-axis limits
# partialPlot(rf_model, pred.data = subset, x.var = "prefire.trees", main = "prefire.trees", ylim = y_limits)
# partialPlot(rf_model, pred.data = subset, x.var = "prefire.below", main = "prefire.below", ylim = y_limits)
# partialPlot(rf_model, pred.data = subset, x.var = "CMI_sp.d", main = "CMI_sp.d", ylim = y_limits)
# partialPlot(rf_model, pred.data = subset, x.var = "CMI_sm.d", main = "CMI_sm.d", ylim = y_limits)
# partialPlot(rf_model, pred.data = subset, x.var = "wspeed", main = "wspeed", ylim = y_limits)
# partialPlot(rf_model, pred.data = subset, x.var = "moist", main = "moist", ylim = y_limits)
# partialPlot(rf_model, pred.data = subset, x.var = "stdage", main = "stdage", ylim = y_limits)
# partialPlot(rf_model, pred.data = subset, x.var = "density", main = "density", ylim = y_limits)
# partialPlot(rf_model, pred.data = subset, x.var = "evi.yr.prior", main = "evi.yr.prior", ylim = y_limits)
# 
# 
# 
# 
# # Define a list to store models and partial dependence data
# models <- list()
# pdp_data <- list()
# 
# # Define the groups and colors
# abmc_groups <- c("conifer", "mixed", "aspen", "birch")
# colors <- c("cornflowerblue", "indianred3", "goldenrod1", "chocolate2")  # Adjust as desired
# 
# # Fit models and extract partial dependence data for each group
# for (i in seq_along(abmc_groups)) {
#   # Filter data for the current group
#   subset <- filter(model_data, abmc == abmc_groups[i]) %>%
#     select(-c(abmc, DC, ISI))  # Drop unnecessary columns
#   
#   # Fit random forest model
#   set.seed(123)  # Ensure reproducibility
#   models[[abmc_groups[i]]] <- randomForest(comb.bgT ~ ., data = subset, ntree = 500, mtry = 3, importance = TRUE)
#   
#   # Extract partial dependence data for a specific variable (e.g., prefire.trees)
#   pd <- partial(models[[abmc_groups[i]]], pred.var = "prefire.trees", grid.resolution = 50, train = subset)
#   pd$abmc_group <- abmc_groups[i]  # Add the group label
#   pdp_data[[i]] <- pd  # Store the partial dependence data
# 
# }
# 
# # Combine all partial dependence data
# pdp_combined <- bind_rows(pdp_data)
# 
# # Plot the partial dependence plot using ggplot2
# ggplot(pdp_combined, aes(x = prefire.trees, y = yhat, color = abmc_group)) +
#   geom_line(linewidth = 1) +  # Updated from size to linewidth
#   labs(x = "Prefire Trees", y = "Partial Dependence on comb.bgT", color = "Stand Type") +
#   scale_color_manual(values = colors) +
#   theme_minimal() +
#   ggtitle("Partial Dependence of C Loss on Prefire Trees by Stand Type")


#FIT FINAL MODELS####

#Specify variables
subset<-data%>%
  # filter(prefire.below<3000)%>%
  dplyr::select(c(abmc,comb.bgT,prefire.trees,prefire.below,stdage,moist,CMI_sp.d,rh)) #Add in dec.ind for 'all' stand type

#Specify stand type
which<-"conifer"
subset<-filter(subset,abmc==which) #Or don't for fu
#specify which
color<-"cornflowerblue" #cornflowerblue chocolate2 indianred4 goldenrod1
subset<-subset%>%dplyr::select(-c(abmc))

# ggplot(subset,aes(x=CMI_sm.d,y=prefire.below,color=fire_scar,shape=abmc))+geom_point()+theme_bw()

# #Quick var importance
# rf_model <- randomForest(comb.bgT ~ ., data = subset, ntree = 500, mtry = 3, importance = TRUE)
# rf_model$importance

pred_vars<-subset%>%dplyr::select(-comb.bgT)
pred_names<-names(pred_vars)
n<-as.numeric(nrow(subset))

seed<-123
set.seed(seed)

# Updated trainControl with more appropriate cross-validation
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, search = "random")

control_method<-control$method
control_number<-control$number
control_repeats<-control$repeats
control_search<-control$search

ntrees<-500

# Define the tuning grid for mtry if you want to limit the range
tune_grid <- expand.grid(mtry = seq(1, length(pred_vars)))

# Train the model
rf.trained <- train(
  comb.bgT ~ ., 
  data = subset, 
  method = "rf", 
  metric = "RMSE", 
  trControl = control, 
  ntree = ntrees,  # Start with 500 trees; increase if needed
  tuneGrid = tune_grid,  # Limit mtry search if desired
  importance = TRUE
)

#Save model if final model
# saveRDS(rf.trained, paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/aspen_compare/final_rf_trained_",which,"_3000.rds"))


# #Spatial autocorrelation?
# subset$residuals <- subset$comb.bgT - predict(rf.trained, subset)
# library(spdep)
# coords <- as.matrix(subset[, c("longitude", "latitude")])
# nb <- dnearneigh(coords, 0, 5000)  # Define neighbors within 5km
# lw <- nb2listw(nb)
# moran.test(subset$residuals, lw)


rf.best.fit <- rf.trained$finalModel

#These performance metrics evaluate performance across cross validation folds and therefore see how weell the model does at predicting on unseen data (reflecting how well the model performs on unseen data (validation folds during training))
getTrainPerf(rf.trained)

#Evaluate split data accuracy####
set.seed(seed)

# 80% for training, 20% for testing
train_index <- createDataPartition(subset$comb.bgT, p = 0.6, list = FALSE)
train_data <- subset[train_index, ]
test_data <- subset[-train_index, ]

# Train on train_data
rf.trained.split <- train(
  comb.bgT ~ ., 
  data = train_data, 
  method = "rf", 
  metric = "RMSE", 
  trControl = control, 
  ntree = ntrees,
  tuneGrid = tune_grid,
  importance = TRUE
)

# Predict on test set
predictions <- predict(rf.trained.split, newdata = test_data)

# Evaluate performance
test_rmse <- RMSE(predictions, test_data$comb.bgT)
test_r2 <- R2(predictions, test_data$comb.bgT)
test_adj_r2 <- summary(lm(test_data$comb.bgT ~ predictions))$adj.r.squared

# Print
cat("Test RMSE:", test_rmse, "\n")
cat("Test R²:", test_r2, "\n")
cat("Test R² adjusted:", test_adj_r2, "\n")

# #Save test/train model
# saveRDS(rf.trained.split, paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/",which,"_test_train_model.rds"))

#Cross vaalidation performance
# Calculate training performance
#These performance metrics evaluste the best fit models' predictions against actual values. Since we are evaluating this using the training data it might overestimate the strength of the model if the model was overfit to the data (evaluates performance on the training dataset, not unseen data, so it reflects training fit rather than generalization.)
train_performance <- postResample(pred = rf.trained$finalModel$predicted, obs = subset$comb.bgT)
train_RMSE<-train_performance[[1]]
train_Rsquared<-train_performance[[2]]
train_MAE<-train_performance[[3]]

# Get cross-validation performance from the caret model
#This evaluates SD of model stats across all folds at the best mtry
best.mtry<-rf.trained$bestTune$mtry
cv_performance <- rf.trained$results[rf.trained$results$mtry == best.mtry, ]
cv_RMSE<-cv_performance[[2]]
cv_RMSESD<-cv_performance[[5]]
cv_Rsquared<-cv_performance[[3]]
cv_RsquaredSD<-cv_performance[[6]]
cv_MAE<-cv_performance[[4]]
cv_MAESD<-cv_performance[[7]]

list(training_performance = train_performance, cv_performance = cv_performance)

#Save model stats####
model_stats<-data.frame(
  stand_type=which,
  n=n,
  control_method=control_method,
  control_number=control_number,
  control_repeats=control_repeats,
  control_search=control_search,
  ntrees=ntrees,
  best.mtry=best.mtry,
  train_RMSE=train_RMSE,
  train_Rsquared=train_Rsquared,
  train_MAE=train_MAE,
  cv_RMSE=cv_RMSE,
  cv_RMSESD=cv_RMSESD,
  cv_Rsquared=cv_Rsquared,
  cv_RsquaredSD=cv_RsquaredSD,
  cv_MAE=cv_MAE,
  cv_MAESD=cv_MAESD
)

# write.csv(model_stats,paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/final_model_stats_rhV.csv"),row.names=F)

#MONTE CARLO VAR IMP PLOT####

# Parameters
n_mc_sims <- 100 #Number of monte carlo simulations to perform

# Monte Carlo Simulations for Variable Importance
simulated_importance <- replicate(n_mc_sims, {
  # Resample the training data (Monte Carlo sampling)
  boot_data <- subset[sample(nrow(subset), replace = TRUE), ]
  
  # Train random forest model
  rf_sim <- randomForest(
    comb.bgT ~ ., 
    data = boot_data, 
    ntree = 500, 
    mtry = best.mtry, 
    importance = TRUE
  )
  
  # Extract variable importance
  importance(rf_sim, type = 1)  # Type 1 for %IncMSE, type = 2 for IncNodePurity
}, simplify = FALSE)

# Combine all results into a data frame
importance_df <- do.call(rbind, lapply(seq_along(simulated_importance), function(i) {
  data.frame(
    Variable = rownames(simulated_importance[[i]]),
    Importance = simulated_importance[[i]][, 1],  # Extract the single column of importance values
    Iteration = i  # Add an iteration identifier
  )
}))

# Summarize results
importance_summary <- importance_df %>%
  group_by(Variable) %>%
  summarise(
    mean_importance = mean(Importance),
    lower_ci = mean(Importance) - 1.96 * sd(Importance) / sqrt(n_distinct(Iteration)),  # 95% CI lower bound
    upper_ci = mean(Importance) + 1.96 * sd(Importance) / sqrt(n_distinct(Iteration))   # 95% CI upper bound
  )

# write.csv(importance_summary,paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/MC_importance_summary_rhV.csv"),row.names=F)

# Plotting Variable Importance
variable_colors <- c(
  "prefire.trees" = "darkgreen",
  "prefire.below" = "brown4",
  "stdage" = "orange",
  "moist" = "purple",
  "CMI_sp.d" = "red",
  # "CMI_sm.d" = "cyan",
  "rh" = "cyan" #,
  # "dec.ind" = "brown"
)

varimp.labs <- c('prefire.trees' = expression('Above C'),
                 'prefire.below' = expression('Below C'),
                 # 'density' = expression('Density (n m'^-2*')'),
                 'stdage' = expression('Stand Age'),
                 'moist' = expression('Moisture'),
                 # 'evi.yr.prior' = expression('Mean Previous Summer EVI-2'),
                 'CMI_sp.d' = expression('dCMI'),
                 # 'CMI_sm.d' = expression('Summer CMI'),
                 # 'wspeed' = expression('Wind Speed at DOB (m / s)'),
                 'rh' = expression('RH') #,
                 # 'dec.ind' = expression('Deciduous Index')
                 )

var.imp.fig.mc <- ggplot(importance_summary, aes(x = reorder(Variable, mean_importance), y = mean_importance)) +
  # geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),linewidth=.5, width = 0.2, color = "black") +
  # geom_point(size = 2,color="black",fill=color) +  # Mean values
  geom_point(position = position_dodge(width = 0.5), shape = 21, size = 2, fill = color, color = 'black', na.rm=T) + 
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12)
  )+ 
  scale_x_discrete(labels = varimp.labs)
# scale_color_manual(values = variable_colors)+
# guides(color="none")
print(var.imp.fig.mc)

# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/importance_MC_incMSE_rhV.png"), plot = var.imp.fig.mc, width = 2, height = 1.5, dpi = 300)




#BEST FIT VAR IMP PLOTS####
# var.imp.dt <- cbind(data.table(var = rownames(rf.best.fit$importance)), data.table(rf.best.fit$importance))[order(MeanDecreaseAccuracy, decreasing = T)]
var.imp.dt <- cbind(data.table(var = rownames(rf.best.fit$importance)), data.table(rf.best.fit$importance))[order(IncNodePurity, decreasing = T)]

vars <- var.imp.dt$var
var.imp.dt

#Save var imp best fit data
# write.csv(var.imp.dt,paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/bestfit_varimportance_rhV.csv"),row.names=F)
var.imp.dt<-read.csv(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/bestfit_varimportance_rhV.csv"))
names(var.imp.dt)<-c("var","%IncMSE","IncNodePurity")
var.imp.dt=as.data.table(var.imp.dt)

var.imp.smry.top6.dt <- var.imp.dt[1:6,]
var.imp.smry.top6.dt[, var := factor(var, levels = rev(var))]                    

# var.imp.fig 
var.imp.fig.inc <- ggplot(var.imp.smry.top6.dt, aes(x = var, y = IncNodePurity))  +
  geom_point(position = position_dodge(width = 0.5), shape = 21, size = 2, fill = color, color = 'black', na.rm=T) + 
  coord_flip() + labs(x = NULL, y = NULL) + 
  scale_x_discrete(labels = varimp.labs) + 
  theme_classic() + theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 8),
                          axis.title = element_text(size = 12))
var.imp.fig.inc

# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/importance_bestfit_incNP_rhV2.png"), plot = var.imp.fig.inc, width = 2, height = 1.5, dpi = 300)

var.imp.smry.top6.dt.mse<-var.imp.smry.top6.dt%>%arrange((`%IncMSE`))%>%
  mutate(var = factor(var, levels = var)) 
var.imp.fig.mse <- ggplot(var.imp.smry.top6.dt.mse, aes(x = var, y = `%IncMSE`))  +
  geom_point(position = position_dodge(width = 0.5), shape = 21, size = 2, fill = color, color = 'black', na.rm=T) + 
  coord_flip() + labs(x =NULL, y = NULL) + 
  scale_x_discrete(labels = varimp.labs) + 
  theme_classic() + theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 8),
                          axis.title = element_text(size = 12))
var.imp.fig.mse

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/importance_bestfit_incMSE_rhV2.png"), plot = var.imp.fig.mse, width = 2, height = 1.5, dpi = 300)


var.imp.fig.mse.wide <- var.imp.fig.mse +
  theme(plot.margin = margin(5, 10, 5, 5))
ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/importance_bestfit_incMSE_rhV2wide.png"), plot = var.imp.fig.mse.wide, width = 2.05, height = 1.5, dpi = 300)

#Partial dependency plots####

#First extract pd estimates
pd.list <- list()
cnt = 1
for (i in vars){
  pd <- data.table(pdp::partial(object = rf.trained, pred.var = i, train = subset, smooth = T, type = 'regression', trim.outliers = T))
  names(pd) <- c('x.value','y.value')
  pd[, var := i]
  pd.list[[cnt]] <- pd
  cnt = cnt + 1
  print(paste('finished',i,sep=' '))
}

pd.dt <- rbindlist(pd.list)   
summary(pd.dt)

# PARTIAL DEPENDENCIES -----------------------------------------


# Define important variables and their labels
top_vars <- as.character(var.imp.smry.top6.dt[1:6]$var)  # Get top variables based on importance
top_vars
pd.labs <- c(   
  'Below C',
  'Above C',
  'Stand Age',
  'dCMI',
  'RH',
  'Moisture') #COnifer
#'Deciduous Index','Spring CMI','Summer CMI','Prefire SOL C','Prefire Tree C','Stand Age','Moisture'
names(pd.labs) <- top_vars  # Associate labels with variable names

# Calculate partial dependence for each important variable, handling factors
pd_list <- lapply(top_vars, function(var) {
  if (is.numeric(subset[[var]])) {
    # Numeric variables: use partial function directly
    pd <- pdp::partial(rf.trained$finalModel, pred.var = var, train = subset, grid.resolution = 50)
  } else {
    # Factor variables: set custom grid for levels with the correct column name
    levels_grid <- data.frame(setNames(list(levels(subset[[var]])), var))
    pd <- pdp::partial(rf.trained$finalModel, pred.var = var, train = subset, pred.grid = levels_grid)
  }
  pd$var <- var  # Add variable name for identification
  return(pd)
})

# Combine partial dependence data
pd_combined <- rbindlist(pd_list, use.names = TRUE, fill = TRUE)
pd_combined$var <- factor(pd_combined$var, levels = top_vars, labels = pd.labs)

# Define `trend.cols` with colors for each variable or class
trend.cols <- c("mediumseagreen", "hotpink", "dodgerblue", "coral", "purple", 
                "darkorange", "deepskyblue",  "firebrick","hotpink")

# Reshape pd_combined to long format
pd_long <- melt(pd_combined, 
                id.vars = c("yhat", "var"), 
                measure.vars = top_vars, 
                variable.name = "x.name", 
                value.name = "x.value")

# Filter out rows with NA in x.value
pd_long <- pd_long[!is.na(pd_long$x.value), ]

# Determine the overall y-axis range for `yhat`
y_range <- range(pd_long$yhat, na.rm = TRUE)

# Plot using ggplot2 with labeller set to as_labeller
pd_fig <- ggplot(pd_long, aes(x = x.value, y = yhat, color = var)) + 
  geom_line(size = 1.5, show.legend = FALSE) +  # Remove legend
  facet_wrap(~ var, scales = "free_x") +  # Custom labels in facet boxes
  xlab("Value of predictor variable") + 
  ylab("Partial Dependence on C Loss") +
  ylim(y_range) +  # Consistent y-axis range
  scale_color_manual(values = trend.cols) + 
  theme_bw() + 
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 10), 
    axis.title = element_text(size = 12)
  )

pd_fig

# Make confidence intervals for PDPs####

#PRINT PDPs####
#Read in model (and to rerun make sure 'subset' data doesn't contain abmc variable)
which<-"aspen"
subset<-data%>%
  filter(prefire.below<3000)%>%
  dplyr::select(abmc,comb.bgT,prefire.trees,prefire.below,stdage,moist,CMI_sp.d,rh)%>%
  filter(abmc==which)%>%
  dplyr::select(comb.bgT,prefire.trees,prefire.below,stdage,moist,CMI_sp.d,rh)

# subset<-filter(subset,prefire.below<2500)
# rf.trained<-readRDS( paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/final_rf_trained_",which,"_.rds"))
# best.mtry<-rf.trained$bestTune$mtry
# 
#conifer
min<-2268
max<-5571
color<-"cornflowerblue"
#mixed
min<-1572
max<-3656
color<-"indianred4"
#aspen
# min<-522
min<-699
# max<-1996
max<-2241
color<-"goldenrod1"
#birch
min<-505
max<-2617
color<-"chocolate2"


##DECIDUOUS INDEX####

# # Calculate the main PDP line (smoothed)
# main_pdp <- partial(
#   rf.trained$finalModel,
#   pred.var = "dec.ind",
#   train = subset,
#   grid.resolution = 50,
#   smooth = TRUE
# )
# 
# Bootstrap partial dependence calculations for confidence intervals
# set.seed(123)
# n_boot <- 100
# partial_list <- replicate(n_boot, {
#   boot_data <- subset[sample(nrow(subset), replace = TRUE), ]
#   rf_boot <- randomForest(comb.bgT ~ ., data = boot_data, ntree = 500,mtry=best.mtry,importance=T)
#   partial(rf_boot, pred.var = "dec.ind", train = boot_data, grid.resolution = 50)
# }, simplify = FALSE)
# 
# # Combine bootstrapped results into one data frame
# pd_combined <- do.call(rbind, partial_list)
# pd_combined <- as.data.frame(pd_combined)
# 
# # Calculate confidence intervals for each predictor value
# pd_summary <- pd_combined %>%
#   group_by(dec.ind) %>%
#   summarise(
#     lower_ci = quantile(yhat, 0.025),
#     upper_ci = quantile(yhat, 0.975)
#   )
# 
# # Add smoothed main PDP line to the summary
# pd_summary.dec.ind <- pd_summary %>%
#   left_join(main_pdp, by = "dec.ind") %>%
#   rename(main_line = yhat)
# pd_summary.dec.ind<-filter(pd_summary.dec.ind,!is.na(main_line))
# 
# write.csv(pd_summary.dec.ind,paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/pd_CIs_dec.ind.csv"),row.names=F)
# pd_summary.dec.ind<-read.csv(paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/pd_CIs_dec.ind.csv"))
# 
# # Plot the PDP with confidence intervals
# pdp.dec.ind<-
#   ggplot(pd_summary.dec.ind, aes(x = dec.ind)) +
#   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = color) +  # CI shading
#   geom_line(aes(y = main_line), size = 1.5, color = color) +  # Main PDP line
#   labs(x = "Deciduous Index", y = NULL) +
#   theme_classic() +
#   theme(
#     strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
#     axis.text = element_text(size = 6),
#     axis.title = element_text(size = 8)
#   )+
#   theme(rect = element_rect(fill = "transparent"),
#         # plot.background = element_rect(fill = "transparent"),
#         panel.background = element_rect(fill = "transparent",
#                                         colour = NA_character_), # necessary to avoid drawing panel outline
#         panel.grid.major = element_blank(), # get rid of major grid
#         panel.grid.minor = element_blank(), # get rid of minor grid
#         plot.background = element_rect(fill = "transparent",
#                                        colour = NA_character_)) # necessary to avoid drawing plot outline)

##SPRING CMI####

# Calculate the main PDP line (smoothed)
main_pdp <- pdp::partial(
  rf.trained$finalModel,
  pred.var = "CMI_sp.d",
  train = subset,
  grid.resolution = 50,
  smooth = TRUE
)

# Bootstrap partial dependence calculations for confidence intervals
set.seed(123)
n_boot <- 100
partial_list <- replicate(n_boot, {
  boot_data <- subset[sample(nrow(subset), replace = TRUE), ]
  rf_boot <- randomForest(comb.bgT ~ ., data = boot_data, ntree = 500,mtry=best.mtry,importance=T)
  pdp::partial(rf_boot, pred.var = "CMI_sp.d", train = boot_data, grid.resolution = 50)
}, simplify = FALSE)

# Combine bootstrapped results into one data frame
pd_combined <- do.call(rbind, partial_list)
pd_combined <- as.data.frame(pd_combined)

# Calculate confidence intervals for each predictor value
pd_summary <- pd_combined %>%
  group_by(CMI_sp.d) %>%
  summarise(
    lower_ci = quantile(yhat, 0.025),
    upper_ci = quantile(yhat, 0.975)
  )

# Add smoothed main PDP line to the summary
pd_summary.cmi.sp <- pd_summary %>%
  left_join(main_pdp, by = "CMI_sp.d") %>%
  rename(main_line = yhat)
pd_summary.cmi.sp<-filter(pd_summary.cmi.sp,!is.na(main_line))
# 
# write.csv(pd_summary.cmi.sp,paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/aspen_compare/pd_CIs_cmi.sp_3000.csv"),row.names=F)
# pd_summary.cmi.sp<-read.csv(paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/pd_CIs_cmi.sp_rhV.csv"))

pd_summary.cmi.sp<-pd_summary.cmi.sp%>%
  mutate(lower_ci=lower_ci/1000)%>%
  mutate(upper_ci=upper_ci/1000)%>%
  mutate(main_line=main_line/1000)
min<-min/1000
max<-max/1000

# x_range <- quantile(pd_summary.cmi.sp$CMI_sp.d, probs = c(0.025, 0.975), na.rm = TRUE)

ggplot(pd_summary.cmi.sp, aes(x = CMI_sp.d, y = main_line)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = color) +
  geom_line(size = 1.5, color = color) +
  # coord_cartesian(xlim = x_range) +  # Clipping extreme values
  theme_classic()+
  ylim(min,max)

# Plot the PDP with confidence intervals
pdp.cmi.sp<-ggplot(pd_summary.cmi.sp, aes(x = CMI_sp.d)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = color) +  # CI shading
  geom_line(aes(y = main_line), size = 1.5, color = color) +  # Main PDP line
  # labs(x = NULL, y = NULL) + #if fully inside
  # labs(x = "Spring dCMI", y = NULL) + #if bottom
  # labs(x = NULL, y = "<2500")+ #paste0(str_to_title(which))) + #if left-most
  labs(x = "dCMI", y=expression(atop("<3 kg m"^-2, "Below C")))+ #paste0(str_to_title(which))) + #if left-most and bottom
  theme_classic() + 
  theme(
    axis.text.y = element_text(size = 6), #If left-most
    axis.text.x = element_text(size = 6),
    # axis.text.y = element_blank(), #if inside (right)
    axis.title = element_text(size = 10) 
  )+
  # scale_y_continuous(breaks = seq(min, max, length.out = 5))+
  ylim(min,max)+ 
  theme(rect = element_rect(fill = "transparent"),
        # plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_)) # necessary to avoid drawing plot outline)
  # +xlim(x_range[[1]],x_range[[2]])
pdp.cmi.sp

##RH####

# Calculate the main PDP line (smoothed)
main_pdp <- pdp::partial(
  rf.trained$finalModel,
  pred.var = "rh",
  train = subset,
  grid.resolution = 50,
  smooth = TRUE
)

# Bootstrap partial dependence calculations for confidence intervals
set.seed(123)
n_boot <- 100
partial_list <- replicate(n_boot, {
  boot_data <- subset[sample(nrow(subset), replace = TRUE), ]
  rf_boot <- randomForest(comb.bgT ~ ., data = boot_data, ntree = 500,mtry=best.mtry,importance=T)
  pdp::partial(rf_boot, pred.var = "rh", train = boot_data, grid.resolution = 50)
}, simplify = FALSE)

# Combine bootstrapped results into one data frame
pd_combined <- do.call(rbind, partial_list)
pd_combined <- as.data.frame(pd_combined)

# Calculate confidence intervals for each predictor value
pd_summary <- pd_combined %>%
  group_by(rh) %>%
  summarise(
    lower_ci = quantile(yhat, 0.025),
    upper_ci = quantile(yhat, 0.975)
  )

# Add smoothed main PDP line to the summary
pd_summary.rh <- pd_summary %>%
  left_join(main_pdp, by = "rh") %>%
  rename(main_line = yhat)
pd_summary.rh<-filter(pd_summary.rh,!is.na(main_line))
# 
# write.csv(pd_summary.rh,paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/aspen_compare/pd_CIs_rh_3000.csv"),row.names=F)
# pd_summary.rh<-read.csv(paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/pd_CIs_rh_rhV.csv"))

#If you want on kg rather than g scale
pd_summary.rh<-pd_summary.rh%>%
  mutate(lower_ci=lower_ci/1000)%>%
  mutate(upper_ci=upper_ci/1000)%>%
  mutate(main_line=main_line/1000)

# Plot the PDP with confidence intervals
pdp.rh<-ggplot(pd_summary.rh, aes(x = rh)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = color) +  # CI shading
  geom_line(aes(y = main_line), size = 1.5, color =color) +  # Main PDP line
  # labs(x = NULL, y = NULL) + #if fully inside
  labs(x = "RH", y = NULL) + #if bottom
  # labs(x = NULL, y = paste0(str_to_title(which))) + #if left-most
  # labs(x = "RH", y = paste0(str_to_title(which))) + #if left-most and bottom
  theme_classic() + 
  theme(
    # axis.text.y = element_text(size = 6), #If left-most
    axis.text.x = element_text(size = 6),
    axis.text.y = element_blank(), #if inside (right)
    axis.title = element_text(size = 10)
  )+ylim(min,max)+ 
  theme(rect = element_rect(fill = "transparent"),
        # plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_)) # necessary to avoid drawing plot outline)
pdp.rh

##PREFIRE SOL C####

# Calculate the main PDP line (smoothed)
main_pdp <- pdp::partial(
  rf.trained$finalModel,
  pred.var = "prefire.below",
  train = subset,
  grid.resolution = 50,
  smooth = TRUE
)

# Bootstrap partial dependence calculations for confidence intervals
set.seed(123)
n_boot <- 100
partial_list <- replicate(n_boot, {
  boot_data <- subset[sample(nrow(subset), replace = TRUE), ]
  rf_boot <- randomForest(comb.bgT ~ ., data = boot_data, ntree = 500,mtry=best.mtry,importance=T)
  pdp::partial(rf_boot, pred.var = "prefire.below", train = boot_data, grid.resolution = 50)
}, simplify = FALSE)

# Combine bootstrapped results into one data frame
pd_combined <- do.call(rbind, partial_list)
pd_combined <- as.data.frame(pd_combined)

# Calculate confidence intervals for each predictor value
pd_summary <- pd_combined %>%
  group_by(prefire.below) %>%
  summarise(
    lower_ci = quantile(yhat, 0.025),
    upper_ci = quantile(yhat, 0.975)
  )

# Add smoothed main PDP line to the summary
pd_summary.sol.c <- pd_summary %>%
  left_join(main_pdp, by = "prefire.below") %>%
  rename(main_line = yhat)
pd_summary.sol.c<-filter(pd_summary.sol.c,!is.na(main_line))
# 
# write.csv(pd_summary.sol.c,paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/aspen_compare/pd_CIs_solC_3000.csv"),row.names=F)
# pd_summary.sol.c<-read.csv(paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/pd_CIs_solC_rhV.csv"))

#If you want on kg rather than g scale
pd_summary.sol.c<-pd_summary.sol.c%>%
  mutate(lower_ci=lower_ci/1000)%>%
  mutate(upper_ci=upper_ci/1000)%>%
  mutate(main_line=main_line/1000)%>%
  mutate(prefire.below=prefire.below/1000)


# x_range <- quantile(pd_summary.sol.c$prefire.below, probs = c(0.025, 0.975), na.rm = TRUE)

ggplot(pd_summary.sol.c, aes(x = prefire.below, y = main_line)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = color) +
  geom_line(size = 1.5, color = color) +
  # xlim(x_range[[1]],x_range[[2]]) +  # Clipping extreme values
  theme_classic()+
  ylim(min,max)

# Plot the PDP with confidence intervals
pdp.sol.c<-ggplot(pd_summary.sol.c, aes(x = prefire.below)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = color) +  # CI shading
  geom_line(aes(y = main_line), size = 1.5, color = color) +  # Main PDP line
  # labs(x = NULL, y = NULL) + #if fully inside
  labs(x = "Below C", y = NULL) + #if bottom
  # labs(x = NULL, y = paste0(str_to_title(which)," stands")) + #if left-most
  # labs(x = "BG SOL C", y = paste0(str_to_title(which)," stands")) + #if left-most and bottom
  theme_classic() + 
  theme(
    # axis.text.y = element_text(size = 6), #If left-most
    axis.text.x = element_text(size = 6),
    axis.text.y = element_blank(), #if inside (right)
    axis.title = element_text(size = 10)
  )+ylim(min,max)+ 
  theme(rect = element_rect(fill = "transparent"),
        # plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))#+ # necessary to avoid drawing plot outline)
  # xlim(x_range[[1]],x_range[[2]])
pdp.sol.c

##PREFIRE TREE C####

# Calculate the main PDP line (smoothed)
main_pdp <- pdp::partial(
  rf.trained$finalModel,
  pred.var = "prefire.trees",
  train = subset,
  grid.resolution = 50,
  smooth = TRUE
)

# Bootstrap partial dependence calculations for confidence intervals
set.seed(123)
n_boot <- 100
partial_list <- replicate(n_boot, {
  boot_data <- subset[sample(nrow(subset), replace = TRUE), ]
  rf_boot <- randomForest(comb.bgT ~ ., data = boot_data, ntree = 500,mtry=best.mtry,importance=T)
  pdp::partial(rf_boot, pred.var = "prefire.trees", train = boot_data, grid.resolution = 50)
}, simplify = FALSE)

# Combine bootstrapped results into one data frame
pd_combined <- do.call(rbind, partial_list)
pd_combined <- as.data.frame(pd_combined)

# Calculate confidence intervals for each predictor value
pd_summary <- pd_combined %>%
  group_by(prefire.trees) %>%
  summarise(
    lower_ci = quantile(yhat, 0.025),
    upper_ci = quantile(yhat, 0.975)
  )

# Add smoothed main PDP line to the summary
pd_summary.tree.c <- pd_summary %>%
  left_join(main_pdp, by = "prefire.trees") %>%
  rename(main_line = yhat)
pd_summary.tree.c<-filter(pd_summary.tree.c,!is.na(main_line))
# 
# write.csv(pd_summary.tree.c,paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/aspen_compare/pd_CIs_treeC_3000.csv"),row.names=F)
# pd_summary.tree.c<-read.csv(paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/pd_CIs_treeC_rhV.csv"))

#If you want on kg rather than g scale
pd_summary.tree.c<-pd_summary.tree.c%>%
  mutate(lower_ci=lower_ci/1000)%>%
  mutate(upper_ci=upper_ci/1000)%>%
  mutate(main_line=main_line/1000)%>%
  mutate(prefire.trees=prefire.trees/1000)

# Plot the PDP with confidence intervals
pdp.tree.c<-ggplot(pd_summary.tree.c, aes(x = prefire.trees)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = color) +  # CI shading
  geom_line(aes(y = main_line), size = 1.5, color = color) +  # Main PDP line
  # labs(x = NULL, y = NULL) + #if fully inside
  labs(x = "Above C", y = NULL) + #if bottom
  # labs(x = NULL, y = paste0(str_to_title(which)," stands")) + #if left-most
  # labs(x = "AG SOL C", y = paste0(str_to_title(which)," stands")) + #if left-most and bottom
  theme_classic() + 
  theme(
    # axis.text.y = element_text(size = 6), #If left-most
    axis.text.x = element_text(size = 6),
    axis.text.y = element_blank(), #if inside (right)
    axis.title = element_text(size = 10)
  )+ylim(min,max)+ 
  theme(rect = element_rect(fill = "transparent"),
        # plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_)) # necessary to avoid drawing plot outline)
pdp.tree.c

##STAND AGE####

# Calculate the main PDP line (smoothed)
main_pdp <- pdp::partial(
  rf.trained$finalModel,
  pred.var = "stdage",
  train = subset,
  grid.resolution = 50,
  smooth = TRUE
)
#
# # Bootstrap partial dependence calculations for confidence intervals
set.seed(123)
n_boot <- 100
partial_list <- replicate(n_boot, {
  boot_data <- subset[sample(nrow(subset), replace = TRUE), ]
  rf_boot <- randomForest(comb.bgT ~ ., data = boot_data, ntree = 500,mtry=best.mtry,importance=T)
  pdp::partial(rf_boot, pred.var = "stdage", train = boot_data, grid.resolution = 50)
}, simplify = FALSE)

# Combine bootstrapped results into one data frame
pd_combined <- do.call(rbind, partial_list)
pd_combined <- as.data.frame(pd_combined)

# Calculate confidence intervals for each predictor value
pd_summary <- pd_combined %>%
  group_by(stdage) %>%
  summarise(
    lower_ci = quantile(yhat, 0.025),
    upper_ci = quantile(yhat, 0.975)
  )

# Add smoothed main PDP line to the summary
pd_summary.stdage <- pd_summary %>%
  left_join(main_pdp, by = "stdage") %>%
  rename(main_line = yhat)
pd_summary.stdage<-filter(pd_summary.stdage,!is.na(main_line))
# 
# write.csv(pd_summary.stdage,paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/aspen_compare/pd_CIs_stdage_3000.csv"),row.names=F)
# pd_summary.stdage<-read.csv(paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/pd_CIs_stdage_rhV.csv"))

#If you want on kg rather than g scale
pd_summary.stdage<-pd_summary.stdage%>%
  mutate(lower_ci=lower_ci/1000)%>%
  mutate(upper_ci=upper_ci/1000)%>%
  mutate(main_line=main_line/1000)

# Plot the PDP with confidence intervals
pdp.stdage<-ggplot(pd_summary.stdage, aes(x = stdage)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = color) +  # CI shading
  geom_line(aes(y = main_line), size = 1.5, color = color) +  # Main PDP line
  # labs(x = NULL, y = NULL) + #if fully inside
  labs(x = "Stand Age", y = NULL) + #if bottom
  # labs(x = NULL, y = paste0(str_to_title(which)," stands")) + #if left-most
  # labs(x = "STand Age", y = paste0(str_to_title(which)," stands")) + #if left-most and bottom
  theme_classic() + 
  theme(
    # axis.text.y = element_text(size = 6), #If left-most
    axis.text.x = element_text(size = 6),
    axis.text.y = element_blank(), #if inside (right)
    axis.title = element_text(size = 10)
  )+ylim(min,max)+ 
  theme(rect = element_rect(fill = "transparent"),
        # plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_)) # necessary to avoid drawing plot outline)
pdp.stdage

##MOISTURE####

# Calculate the main PDP line (smoothed)
main_pdp <- pdp::partial(
  rf.trained$finalModel,
  pred.var = "moist",
  train = subset,
  grid.resolution = 50,
  smooth = TRUE
)

# Bootstrap partial dependence calculations for confidence intervals
set.seed(123)
n_boot <- 100
partial_list <- replicate(n_boot, {
  boot_data <- subset[sample(nrow(subset), replace = TRUE), ]
  rf_boot <- randomForest(comb.bgT ~ ., data = boot_data, ntree = 500,mtry=best.mtry,importance=T)
  pdp::partial(rf_boot, pred.var = "moist", train = boot_data, grid.resolution = 50)
}, simplify = FALSE)

# Combine bootstrapped results into one data frame
pd_combined <- do.call(rbind, partial_list)
pd_combined <- as.data.frame(pd_combined)

# Calculate confidence intervals for each predictor value
pd_summary <- pd_combined %>%
  group_by(moist) %>%
  summarise(
    lower_ci = quantile(yhat, 0.025),
    upper_ci = quantile(yhat, 0.975)
  )

# Add smoothed main PDP line to the summary
pd_summary.moist <- pd_summary %>%
  left_join(main_pdp, by = "moist") %>%
  rename(main_line = yhat)
pd_summary.moist<-filter(pd_summary.moist,!is.na(main_line))
# 
# write.csv(pd_summary.moist,paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/aspen_compare/pd_CIs_moist_3000.csv"),row.names=F)
# pd_summary.moist<-read.csv(paste0(data_directory,"/model_outputs/Q2_RF_models/",which,"/pd_CIs_moist_rhV.csv"))

#If you want on kg rather than g scale
pd_summary.moist<-pd_summary.moist%>%
  mutate(lower_ci=lower_ci/1000)%>%
  mutate(upper_ci=upper_ci/1000)%>%
  mutate(main_line=main_line/1000)

# Plot the PDP with confidence intervals
pdp.moist<-ggplot(pd_summary.moist, aes(x = moist)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = color) +  # CI shading
  geom_line(aes(y = main_line), size = 1.5, color = color) +  # Main PDP line
  # labs(x = NULL, y = NULL) + #if fully inside
  labs(x = "Moisture", y = NULL) + #if bottom
  # labs(x = NULL, y = paste0(str_to_title(which)," stands")) + #if left-most
  # labs(x = "MMoisture", y = paste0(str_to_title(which)," stands")) + #if left-most and bottom
  theme_classic() + 
  theme(
    # axis.text.y = element_text(size = 6), #If left-most
    axis.text.x = element_text(size = 6),
    axis.text.y = element_blank(), #if inside (right)
    axis.title = element_text(size = 10)
  )+ylim(min,max)+ 
  theme(rect = element_rect(fill = "transparent"),
        # plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_)) # necessary to avoid drawing plot outline)
pdp.moist
# 
# #For arranging grid extra plot
# 
# # pdp.dec.ind.con<-pdp.dec.ind
# pdp.cmi.sp.bir<-pdp.cmi.sp
# pdp.cmi.sm.bir<-pdp.cmi.sm
# pdp.sol.c.bir<-pdp.sol.c
# pdp.tree.c.bir<-pdp.tree.c
# pdp.stdage.bir<-pdp.stdage
# pdp.moist.bir<-pdp.moist
# 
# #Add numbers
# library(cowplot)
# pdp.cmi.sp.test<-pdp.cmi.sp+
#   draw_label(
#     "1", 
#     # x = 0.02,   # x position relative to plot (0 = far left, 1 = far right)
#     # y = 0.98,   # y position relative to plot (0 = bottom, 1 = top)
#     hjust = 0,  # Left-align the label
#     vjust = 1,  # Top-align the label
#     size = 4   # Text size
#   )
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/test.png"), plot = pdp.cmi.sp.test, width = 1.5, height = 1.5, dpi = 300)
# 
# pdp.cmi.sp.bir<-pdp.cmi.sp.bir+annotate("text",x = Inf, y = Inf,label = "1",size = 4,hjust = 1.1, vjust = 1.1)
# pdp.cmi.sm.bir<-pdp.cmi.sm.bir+annotate("text",x = Inf, y = Inf,label = "4",size = 4,hjust = 1.1, vjust = 1.1)
# pdp.sol.c.bir<-pdp.sol.c.bir+annotate("text",x = Inf, y = Inf,label = "2",size = 4,hjust = 1.1, vjust = 1.1)
# pdp.tree.c.bir<-pdp.tree.c.bir+annotate("text",x = Inf, y = Inf,label = "3",size = 4,hjust = 1.1, vjust = 1.1)
# pdp.stdage.bir<-pdp.stdage.bir+annotate("text",x = Inf, y = Inf,label = "5",size = 4,hjust = 1.1, vjust = 1.1)
# pdp.moist.bir<-pdp.moist.bir+annotate("text",x = Inf, y = Inf,label = "6",size = 4,hjust = 1.1, vjust = 1.1)
# 
# library(gridExtra)
# grid<-grid.arrange(
#   pdp.tree.c.con,pdp.tree.c.mix,pdp.tree.c.asp,pdp.tree.c.bir,
#   pdp.sol.c.con,pdp.sol.c.mix,pdp.sol.c.asp,pdp.sol.c.bir,
#   pdp.moist.con,pdp.moist.mix,pdp.moist.asp,pdp.moist.bir,
#   pdp.stdage.con,pdp.stdage.mix,pdp.stdage.asp,pdp.stdage.bir,
#   pdp.cmi.sp.con,pdp.cmi.sp.mix,pdp.cmi.sp.asp,pdp.cmi.sp.bir,
#   pdp.cmi.sm.con,pdp.cmi.sm.mix,pdp.cmi.sm.asp,pdp.cmi.sm.bir,
#   ncol = 4, nrow = 6
# )
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/plotgrid.png"), plot = grid, width = 5.1, height = 7, dpi = 300)
# 
# #Adjust gridExtra#Adjust axes and save PDPs####
# # 
# # #All
# # min<-1300
# # max<-3100
# 
# #Across all plots
# min<-400
# max<-5450
# 
# pdp.dec.ind<-pdp.dec.ind+ylim(min,max)
# pdp.cmi.sp<-pdp.cmi.sp+ylim(min,max)
# pdp.cmi.sm<-pdp.cmi.sm+ylim(min,max)
# pdp.sol.c<-pdp.sol.c+ylim(min,max)
# pdp.tree.c<-pdp.tree.c+ylim(min,max)
# pdp.stdage<-pdp.stdage+ylim(min,max)
# pdp.moist<-pdp.moist+ylim(min,max)
# 
# ##V1####
# #Save versions with updated ylims
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp2_dec.ind.png"), plot = pdp.dec.ind, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp2_CMI.sp.png"), plot = pdp.cmi.sp, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp2_CMI.sm.png"), plot = pdp.cmi.sm, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp2_sol.c.png"), plot = pdp.sol.c, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp2_tree.c.png"), plot = pdp.tree.c, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp2_stdage.png"), plot = pdp.stdage, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp2_moist.png"), plot = pdp.moist, width = 1.5, height = 1.5, dpi = 300)
# 
# #################################################
# #V2####
# #Saver original versions
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp_dec.ind.png"), plot = pdp.dec.ind, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp_CMI.sp.png"), plot = pdp.cmi.sp, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp_CMI.sm.png"), plot = pdp.cmi.sm, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp_sol.c.png"), plot = pdp.sol.c, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp_tree.c.png"), plot = pdp.tree.c, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp_stdage.png"), plot = pdp.stdage, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp_moist.png"), plot = pdp.moist, width = 1.5, height = 1.5, dpi = 300)
# 
# 
# #V3####
# #Saver only captions left and bottom versions
# # ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_dec.ind.png"), plot = pdp.dec.ind, width = 1.5, height = 1.5, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sp.png"), plot = pdp.cmi.sp, width = 1.237, height = 1.11, dpi = 300) #if left axis does NNOT have a decimal pt
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sp.png"), plot = pdp.cmi.sp, width = 1.3, height = 1.11, dpi = 300) #if left axis has a decimal pt
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sm.png"), plot = pdp.cmi.sm, width = 1, height = 1.11, dpi = 300, bg = 'transparent')
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_sol.c.png"), plot = pdp.sol.c, width = 1, height = 1.11, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_tree.c.png"), plot = pdp.tree.c, width = 1, height = 1.11, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_stdage.png"), plot = pdp.stdage, width = 1, height = 1.11, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_moist.png"), plot = pdp.moist, width = 1, height = 1.11, dpi = 300)
# 
# #For bottom row:
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sp.png"), plot = pdp.cmi.sp, width = 1.22, height = 1.27, dpi = 300) #if left axis does NNOT have a decimal pt
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sp.png"), plot = pdp.cmi.sp, width = 1.3, height = 1.27, dpi = 300) #if left axis has a decimal pt
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sm.png"), plot = pdp.cmi.sm, width = 1, height = 1.27, dpi = 300, bg = 'transparent')
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_sol.c.png"), plot = pdp.sol.c, width = 1, height = 1.27, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_tree.c.png"), plot = pdp.tree.c, width = 1, height = 1.27, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_stdage.png"), plot = pdp.stdage, width = 1, height = 1.27, dpi = 300)
# 
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_moist.png"), plot = pdp.moist, width = 1, height = 1.27, dpi = 300)



#V4####
#Saver only captions left and bottom versions
# ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_dec.ind.png"), plot = pdp.dec.ind, width = 1.5, height = 1.5, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sp_rhV.png"), plot = pdp.cmi.sp, width = 1.237, height = 1.11, dpi = 300) #if left axis does NNOT have a decimal pt
ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sp_rhV.png"), plot = pdp.cmi.sp, width = 1.3, height = 1.11, dpi = 300) #if left axis has a decimal pt

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_rh_rhV.png"), plot = pdp.rh, width = 1, height = 1.11, dpi = 300, bg = 'transparent')

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_sol.c_rhV.png"), plot = pdp.sol.c, width = 1, height = 1.11, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_tree.c_rhV.png"), plot = pdp.tree.c, width = 1, height = 1.11, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_stdage_rhV.png"), plot = pdp.stdage, width = 1, height = 1.11, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_moist_rhV.png"), plot = pdp.moist, width = 1, height = 1.11, dpi = 300)

#For bottom row:
ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sp_rhV.png"), plot = pdp.cmi.sp, width = 1.22, height = 1.27, dpi = 300) #if left axis does NNOT have a decimal pt
ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_CMI.sp_rhV.png"), plot = pdp.cmi.sp, width = 1.3, height = 1.27, dpi = 300) #if left axis has a decimal pt

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_rh_rhV.png"), plot = pdp.rh, width = 1, height = 1.27, dpi = 300, bg = 'transparent')

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_sol.c_rhV.png"), plot = pdp.sol.c, width = 1, height = 1.27, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_tree.c_rhV.png"), plot = pdp.tree.c, width = 1, height = 1.27, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_stdage_rhV.png"), plot = pdp.stdage, width = 1, height = 1.27, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/pdp3_moist_rhV.png"), plot = pdp.moist, width = 1, height = 1.27, dpi = 300)

#EXPORT COMPARE####


ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_CMI.sp_2500.png"), plot = pdp.cmi.sp, width = 1.237, height = 1.11, dpi = 300) #if left axis does NNOT have a decimal pt
ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_CMI.sp_2500.png"), plot = pdp.cmi.sp, width = 1.3, height = 1.11, dpi = 300) #if left axis has a decimal pt

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_rh_2500.png"), plot = pdp.rh, width = 1, height = 1.11, dpi = 300, bg = 'transparent')

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_sol.c_2500.png"), plot = pdp.sol.c, width = 1, height = 1.11, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_tree.c_2500.png"), plot = pdp.tree.c, width = 1, height = 1.11, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_stdage_2500.png"), plot = pdp.stdage, width = 1, height = 1.11, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_moist_2500.png"), plot = pdp.moist, width = 1, height = 1.11, dpi = 300)

#For bottom row:
ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_CMI.sp_3000.png"), plot = pdp.cmi.sp, width = 1.22, height = 1.27, dpi = 300) #if left axis does NNOT have a decimal pt
ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_CMI.sp_3000v2.png"), plot = pdp.cmi.sp, width = 1.35, height = 1.27, dpi = 300) #if left axis has a decimal pt

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_rh_3000.png"), plot = pdp.rh, width = 1, height = 1.27, dpi = 300, bg = 'transparent')

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_sol.c_3000.png"), plot = pdp.sol.c, width = 1, height = 1.27, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_tree.c_3000.png"), plot = pdp.tree.c, width = 1, height = 1.27, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_stdage_3000.png"), plot = pdp.stdage, width = 1, height = 1.27, dpi = 300)

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/aspen_compare/pdp3_moist_3000.png"), plot = pdp.moist, width = 1, height = 1.27, dpi = 300)



#Raw data plots
which<-"all"

data$abmc<-factor(data$abmc, levels=c("conifer","mixed", "aspen", "birch"))
g=ggplot(data,aes(x=dec.ind,y=comb.bgT))+
  geom_point(aes(color=data$abmc))+
  geom_smooth(method="loess",color="black")+
  scale_color_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  labs(x = "Deciduous Index", y = "Total C Loss") +
  theme(
    strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 6), 
    axis.title = element_text(size = 10)
  )+
  guides(color="none")
g

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/raw_dec.ind.png"), plot = g, width = 3, height = 2, dpi = 300)

g=ggplot(data,aes(x=CMI_sp.d,y=comb.bgT,color=abmc,fill=abmc))+
  geom_point(aes(color=abmc))+
  geom_smooth(method="lm")+
  # geom_smooth(data,aes(color=abmc),method="loess")+
  scale_color_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  labs(x = "Spring CMI", y = "Total C Loss") +
  theme(
    strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 6), 
    axis.title = element_text(size = 10)
  )+
  guides(color="none",fill="none")
g

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/raw_cmi.sp.png"), plot = g, width = 3, height = 2, dpi = 300)

g=ggplot(data,aes(x=CMI_sm.d,y=comb.bgT,color=abmc,fill=abmc))+
  geom_point(aes(color=abmc))+
  geom_smooth(method="lm")+
  # geom_smooth(data,aes(color=abmc),method="loess")+
  scale_color_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  labs(x = "Summer CMI", y = "Total C Loss") +
  theme(
    strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 6), 
    axis.title = element_text(size = 10)
  )+
  guides(color="none",fill="none")
g

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/raw_cmi.sm.png"), plot = g, width = 3, height = 2, dpi = 300)

g=ggplot(data,aes(x=stdage,y=comb.bgT,color=abmc,fill=abmc))+
  geom_point(aes(color=abmc))+
  geom_smooth(method="lm")+
  # geom_smooth(data,aes(color=abmc),method="loess")+
  scale_color_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  labs(x = "Stand Age", y = "Total C Loss") +
  theme(
    strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 6), 
    axis.title = element_text(size = 10)
  )+
  guides(color="none",fill="none")
g

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/raw_stdage.png"), plot = g, width = 3, height = 2, dpi = 300)

g=ggplot(data,aes(x=prefire.trees,y=comb.bgT,color=abmc,fill=abmc))+
  geom_point(aes(color=abmc))+
  geom_smooth(method="lm")+
  # geom_smooth(data,aes(color=abmc),method="loess")+
  scale_color_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  labs(x = "Prefire Tree C", y = "Total C Loss") +
  theme(
    strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 6), 
    axis.title = element_text(size = 10)
  )+
  guides(color="none",fill="none")
g

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/raw_trees.png"), plot = g, width = 3, height = 2, dpi = 300)



g=ggplot(data,aes(x=prefire.below,y=comb.bgT,color=abmc,fill=abmc))+
  geom_point(aes(color=abmc))+
  geom_smooth(method="lm")+
  # geom_smooth(data,aes(color=abmc),method="loess")+
  scale_color_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  labs(x = "Prefire SOL C", y = "Total C Loss") +
  theme(
    strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 6), 
    axis.title = element_text(size = 10)
  )+
  guides(color="none",fill="none")
g

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/raw_sol.png"), plot = g, width = 3, height = 2, dpi = 300)

g=ggplot(data,aes(x=moist,y=comb.bgT,color=abmc,fill=abmc))+
  geom_point(aes(color=abmc))+
  geom_smooth(method="lm")+
  # geom_smooth(data,aes(color=abmc),method="loess")+
  scale_color_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  labs(x = "Moisture", y = "Total C Loss") +
  theme(
    strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 6), 
    axis.title = element_text(size = 10)
  )+
  guides(color="none",fill="none")
g

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/",which,"/raw_moist.png"), plot = g, width = 3, height = 2, dpi = 300)

#EXTRACT SUMMARY CSV####
# Define rounding rules for each variable
rounding_rules <- list(
  comb.bgT = 0,
  prefire.trees = 0,
  prefire.below = 0,
  stdage = 0,
  dec.ind = 0,
  moist = 1,
  CMI_sp.d = 1,
  CMI_sm.d = 1
)

# Define a function to calculate mean ± SE and min-max with variable-specific rounding
calculate_stats <- function(data, rounding_rules) {
  data %>%
    summarise(across(
      where(is.numeric),
      list(
        mean_se = ~ {
          var_name <- cur_column()
          digits <- rounding_rules[[var_name]]
          paste0(
            round(mean(.), digits), 
            " ± ", 
            round(sd(.) / sqrt(n()), digits)
          )
        },
        min_max = ~ {
          var_name <- cur_column()
          digits <- rounding_rules[[var_name]]
          paste0(
            "(", 
            round(min(.), digits), 
            "-", 
            round(max(.), digits), 
            ")"
          )
        }
      ),
      .names = "{col}_{fn}"
    ))
}

# Group by 'abmc' and calculate statistics
summary_stats <- subset %>%
  group_by(abmc) %>%
  calculate_stats(rounding_rules) %>%
  ungroup()

# Save the results to a CSV file
write.csv(summary_stats,paste0(data_directory,"/model_outputs/Q2_RF_models/driver_summary.csv"), row.names = FALSE)


#Visual plots

new_labels <- c(
  "Aggie" = "Aggie\nCreek",
  "Beaver" = "2019BC001",
  "Chena" = "Munson\nCreek",
  "Ethel" = "2019MA014",
  "Isom" = "Isom\nCreek",
  "Livengood" = "Hess\nCreek",
  "Manley" = "Baker",
  "Shovel" = "Shovel\nCreek"
)
g=ggplot(filter(data,!abmc=="conifer"),aes(x=fire_scar,y=comb.bgT,color=abmc))+geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA)+scale_fill_manual(values=c("indianred4","goldenrod1","chocolate2"))+scale_color_manual(values=c("indianred4","goldenrod1","chocolate2"))+geom_jitter( position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1.5, alpha = 0.6)+theme_bw() +
  labs(color = NULL, fill = NULL)+
  theme(
    legend.position = c(0.5, 0.85),       # Move legend to 80% (x) and 80% (y) of the plot area
    legend.background = element_rect(fill = "white", color = "black"),  # Optional: Add background
    legend.box.background = element_blank(),  # Remove extra box
    legend.key = element_rect(fill = NA),
    # strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 8,face="bold"), 
    axis.title = element_text(size = 10)
  )+
  labs(x = "Fire scar", y = "Total C Loss") + 
  scale_x_discrete(labels = new_labels)
g

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/comb_by_fire.png"), plot = g, width = 6.5, height = 4, dpi = 300)



new_labels <- c(
  "Aggie" = "Aggie\nCreek\n06/22/2015",
  "Beaver" = "2019BC005\n07/06/2019",
  "Chena" = "Munson\nCreek\n06/18/2021",
  "Ethel" = "2019MA014\n07/19/2019",
  "Isom" = "Isom\nCreek\n06/05/2020",
  "Livengood" = "Hess\nCreek\n06/21/2019",
  "Manley" = "Baker\n06/21/2015",
  "Shovel" = "Shovel\nCreek\n06/21/2019"
)

data$fire_scar<-factor(data$fire_scar, levels=c("Manley","Aggie", "Livengood", "Shovel","Beaver","Ethel","Isom","Chena"))

sample_sizes <- data %>%
  filter(!abmc == "conifer") %>%
  group_by(fire_scar) %>%
  summarize(n = n())

g=ggplot(filter(data,!abmc=="conifer"),aes(x=fire_scar,y=comb.bgT,color=fire_scar))+geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA)+geom_jitter( position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8), size = 1.5, alpha = 0.6)+theme_bw() +
  labs(color = NULL, fill = NULL)+
  theme(
    # legend.position = c(0.5, 0.85),       # Move legend to 80% (x) and 80% (y) of the plot area
    legend.background = element_rect(fill = "white", color = "black"),  # Optional: Add background
    legend.box.background = element_blank(),  # Remove extra box
    legend.key = element_rect(fill = NA),
    # strip.text = element_text(size = 6, face = "bold"),  # Adjust facet title text in gray box
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 10)
  )+
  labs(x = "Fire Names & Ignition Dates", y = "Total C Loss") + 
  scale_x_discrete(labels = new_labels)+
  guides(color="none",fill="none")+
  geom_text(
    data = sample_sizes, 
    aes(x = fire_scar, y = max(filter(data, !abmc == "conifer")$comb.bgT) + 5, label = paste0("n = ", n)),
    inherit.aes = FALSE,  # Prevent the `aes()` in ggplot() from overriding this layer
    size = 3, color = "black"
  ) 
g

ggsave(paste0(data_directory, "/model_outputs/Q2_RF_models/comb_by_fire.png"), plot = g, width = 6.5, height = 4, dpi = 300)

#Make point plots
ggplot(filter(data,abmc=="aspen"), aes(x = prefire.below)) +
  geom_line(aes(y = comb.below.prop), color = "blue") +
  geom_line(aes(y = comb.below), color = "red") +  # scale z to match y-axis
  scale_y_continuous(
    name = "Y Axis (y)",
    sec.axis = sec_axis(~ . * 10, name = "Z Axis (z)")
  ) +
  theme_minimal()