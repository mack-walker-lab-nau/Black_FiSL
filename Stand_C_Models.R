library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(DHARMa)
library(spaMM)
library(rcdd)
library(glmmTMB)
library(multcomp)
rm(list=ls())

set.seed(123)

#Set theme for figs
theme_set(theme_bw(base_family = "Arial"))
library(extrafont)
loadfonts(device = "pdf")

#Directory
data_directory <- "C:/Users/Test/Documents/MS THESIS/SoilModelRedo"

#Data
data<-read.csv("C:/Users/Test/Documents/MS THESIS/BNZ Data Archiving/FiSL_FullModelingDataset_BB2023.csv")

#Filter Data for Modeling
data<-data %>% 
  filter(!is.na(postfire.below)) %>%  #Missing SOL sample
  filter(!is.na(prefire.below)) #FiSL poplar, mixed-poplar, and conifer data we didn't model prefire SOL for
  
#One lat/long taken per site but plots were 10m apart, simulate 10m intervals for Fisl plots (plots of A, B, or C)
data=data%>%
  mutate(longitude=ifelse(data=="FiSL"&plot=="A",longitude,
                          ifelse(data=="FiSL"&plot=="B",(longitude+((.01/6378.1)*(180/pi))/cos(latitude*pi/108)),
                                 ifelse(data=="FiSL"&plot=="C",(longitude+((.02/6378.1)*(180/pi))/cos(latitude*pi/108)),longitude))))

data=data%>%
  mutate(latitude=ifelse(data=="FiSL"&plot=="A",latitude,
                         ifelse(data=="FiSL"&plot=="B",(latitude+((.01/6378.1)*(180/pi))),
                                ifelse(data=="FiSL"&plot=="C",(latitude+((.02/6378.1)*(180/pi))),latitude
                                ))))

#Set stand type factor levels
data$abmc<-factor(data$abmc, levels=c("conifer","mixed", "aspen", "birch"))

#Create data frame to store outputs
full_output_data<-data.frame()

#ABOVEGROUND####

##PREFIRE TREES####
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,prefire.trees)%>%
  na.omit()

d1<-d%>%
  mutate(prefire.trees=prefire.trees+.1) #Simulate non-zero response metric

###Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(prefire.trees ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(prefire.trees ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(prefire.trees ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(prefire.trees ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)

#Add a random effect of fire scar because spatial effect is not going great
m4 <- fitme(prefire.trees ~ abmc + Matern(1 | latitude + longitude) + (1|fire_scar),
            family = gaussian(),
            data = d,
            resid.model = ~ abmc)
m4_AIC<-extractAIC(m4)

# m5 <- fitme(prefire.trees ~ abmc + Cauchy(1 | latitude + longitude),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m5_AIC<-extractAIC(m5)

# m6 <- fitme(prefire.trees ~ abmc + Cauchy(1 | latitude + longitude) + (1|fire_scar),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m6_AIC<-extractAIC(m6)


m7 <- fitme(prefire.trees ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

# m8 <- fitme(prefire.trees ~ abmc + Cauchy(1 | latitude + longitude),
#             family = Gamma(log),
#             data = d1,
#             resid.model = ~ abmc)
# m8_AIC<-extractAIC(m8)

m9 <- fitme(prefire.trees ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(prefire.trees ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)

# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$prefire.trees<-d_nb$prefire.trees+.001
# d_nb$prefire.trees<-round(d_nb$prefire.trees,0)
# d_nb$prefire.trees<-as.integer(d_nb$prefire.trees)

# m11 <-  fitme(prefire.trees ~ abmc,
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(prefire.trees ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(prefire.trees ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(prefire.trees ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)

###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m7)
testResiduals(m7)

plotResiduals(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"


###Outlier test####
testOutliers(m7)
out_res<-residuals(residuals)

# Check for residuals beyond a certain threshold (depending on the distribution of residuals)
threshold <- quantile(out_res, probs = c(0.025, 0.975))  # for 95% interval
outliers <- which(out_res < threshold[1] | out_res > threshold[2])

# Get the actual values that are outliers
outlier_values <- d1[outliers, ]

d_no_outliers<-d1[-outliers,]

#Re-run model with no outliers
no_outlier_model <-  fitme(prefire.trees ~ abmc + Matern(1 | latitude + longitude),
                           family = Gamma(log),
                           data = d_no_outliers,
                           resid.model = ~ abmc)
no_outlier_model_AIC<-extractAIC(no_outlier_model)



### Extract summary####
# Extract the response variable
winning_model<-m7
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/prefire.trees_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- as.numeric(DoF(winning_model)["p_fixef"])  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_prefire_trees = mean(prefire.trees, na.rm = TRUE),
    se_prefire_trees = sd(prefire.trees, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_prefire_trees, se = se_prefire_trees) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)



##POSTFIRE TREES####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,m11_AIC,m12_AIC,m13_AIC,m14_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,d_nb,d_no_outliers,m12_out,outlier_values,m12_out_AIC,out_res,outliers,threshold,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,postfire.trees)%>%
  na.omit()

d1<-d%>%
  mutate(postfire.trees=postfire.trees+.1)


####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(postfire.trees ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(postfire.trees ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(postfire.trees ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(postfire.trees ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)

#Add a random effect of fire scar because spatial effect is not going great
m4 <- fitme(postfire.trees ~ abmc + Matern(1 | latitude + longitude) + (1|fire_scar),
            family = gaussian(),
            data = d,
            resid.model = ~ abmc)
m4_AIC<-extractAIC(m4) # "did not always converge"
# 
# m5 <- fitme(postfire.trees ~ abmc + Cauchy(1 | latitude + longitude),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m5_AIC<-extractAIC(m5)
# 
# m6 <- fitme(postfire.trees ~ abmc + Cauchy(1 | latitude + longitude) + (1|fire_scar),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m6_AIC<-extractAIC(m6)


m7 <- fitme(postfire.trees ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

# m8 <- fitme(postfire.trees ~ abmc + Cauchy(1 | latitude + longitude),
#             family = Gamma(log),
#             data = d1,
#             resid.model = ~ abmc)
# m8_AIC<-extractAIC(m8)

m9 <- fitme(postfire.trees ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(postfire.trees ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)
# 
# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$postfire.trees<-d_nb$postfire.trees+.001
# d_nb$postfire.trees<-round(d_nb$postfire.trees,0)
# # d_nb$postfire.trees<-as.integer(d_nb$postfire.trees)
# 
# m11 <-  fitme(postfire.trees ~ abmc,
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(postfire.trees ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(postfire.trees ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(postfire.trees ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)


###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m7)
testResiduals(m7)
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))


###Outlier test####
testOutliers(m7)
out_res<-residuals(residuals)

# Check for residuals beyond a certain threshold (depending on the distribution of residuals)
threshold <- quantile(out_res, probs = c(0.025, 0.975))  # for 95% interval
outliers <- which(out_res < threshold[1] | out_res > threshold[2])

# Get the actual values that are outliers
outlier_values <- d1[outliers, ]

d_no_outliers<-d1[-outliers,]

#Re-run model with no outliers
no_outlier_model <-  fitme(postfire.trees ~ abmc + Matern(1 | latitude + longitude),
                           family = Gamma(log),
                           data = d_no_outliers,
                           resid.model = ~ abmc)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d_no_outliers$abmc)),
  latitude = mean(d_no_outliers$latitude),  # Use the mean latitude for prediction
  longitude = mean(d_no_outliers$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(no_outlier_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(no_outlier_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))


### Extract summary####
# Extract the response variable
winning_model<-m7
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/postfire.trees_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(m7)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_postfire_trees = mean(postfire.trees, na.rm = TRUE),
    se_postfire_trees = sd(postfire.trees, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_postfire_trees, se = se_postfire_trees) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)
#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)

##COMB TREES####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,new_data,combined_df,contrast_df,grouped_means_se,predicted_means,log_predictions_with_variances,response_scale_se,log_variances)

d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,comb.trees)%>%
  na.omit()

d1<-d%>%                 
  mutate(comb.trees=comb.trees+.1)


####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(comb.trees ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(comb.trees ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(comb.trees ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(comb.trees ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)

#Add a random effect of fire scar because spatial effect is not going great
m4 <- fitme(comb.trees ~ abmc + Matern(1 | latitude + longitude) + (1|fire_scar),
            family = gaussian(),
            data = d,
            resid.model = ~ abmc)
m4_AIC<-extractAIC(m4) 

# m5 <- fitme(comb.trees ~ abmc + Cauchy(1 | latitude + longitude),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m5_AIC<-extractAIC(m5)

# m6 <- fitme(comb.trees ~ abmc + Cauchy(1 | latitude + longitude) + (1|fire_scar),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m6_AIC<-extractAIC(m6)


m7 <- fitme(comb.trees ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

# m8 <- fitme(comb.trees ~ abmc + Cauchy(1 | latitude + longitude),
#             family = Gamma(log),
#             data = d1,
#             resid.model = ~ abmc)
# m8_AIC<-extractAIC(m8)

m9 <- fitme(comb.trees ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(comb.trees ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)


###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m7)
testResiduals(m7)
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))



### Extract summary####
# Extract the response variable
winning_model<-m7
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.trees_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(m7)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_comb_trees = mean(comb.trees, na.rm = TRUE),
    se_comb_trees = sd(comb.trees, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_comb_trees, se = se_comb_trees) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)
#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)



##PREFIRE BELOW####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,prefire.below)%>%
  na.omit()

d1<-d%>%                 
  mutate(prefire.below=prefire.below+.1)


####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(prefire.below ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(prefire.below ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(prefire.below ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(prefire.below ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)

#Add a random effect of fire scar because spatial effect is not going great
m4 <- fitme(prefire.below ~ abmc + Matern(1 | latitude + longitude) + (1|fire_scar),
            family = gaussian(),
            data = d,
            resid.model = ~ abmc)
m4_AIC<-extractAIC(m4) 

# m5 <- fitme(prefire.below ~ abmc + Cauchy(1 | latitude + longitude),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m5_AIC<-extractAIC(m5)

# m6 <- fitme(prefire.below ~ abmc + Cauchy(1 | latitude + longitude) + (1|fire_scar),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m6_AIC<-extractAIC(m6)


m7 <- fitme(prefire.below ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

# m8 <- fitme(prefire.below ~ abmc + Cauchy(1 | latitude + longitude),
#             family = Gamma(log),
#             data = d1,
#             resid.model = ~ abmc)
# m8_AIC<-extractAIC(m8)

m9 <- fitme(prefire.below ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(prefire.below ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)


###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m7)
testResiduals(m7)
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))



### Extract summary####
# Extract the response variable
winning_model<-m7
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/prefire.below_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(m7)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)


# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_prefire_below = mean(prefire.below, na.rm = TRUE),
    se_prefire_below = sd(prefire.below, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_prefire_below, se = se_prefire_below) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)

##POSTFIRE BELOW####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,postfire.below)%>%
  na.omit()

d1<-d%>%                 
  mutate(postfire.below=postfire.below+.1)
# d2<-d%>%                 
#   mutate(postfire.below=postfire.below+1)
# d2$postfire.below<-round(d2$postfire.below,0)

####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(postfire.below ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(postfire.below ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(postfire.below ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)

#Add a random effect of fire scar because spatial effect is not going great
m4 <- fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude) + (1|fire_scar),
            family = gaussian(),
            data = d,
            resid.model = ~ abmc)
m4_AIC<-extractAIC(m4) 

# m5 <- fitme(postfire.below ~ abmc + Cauchy(1 | latitude + longitude),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m5_AIC<-extractAIC(m5)

# m6 <- fitme(postfire.below ~ abmc + Cauchy(1 | latitude + longitude) + (1|fire_scar),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m6_AIC<-extractAIC(m6)


m7 <- fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

# m8 <- fitme(postfire.below ~ abmc + Cauchy(1 | latitude + longitude),
#             family = Gamma(log),
#             data = d1,
#             resid.model = ~ abmc)
# m8_AIC<-extractAIC(m8)

m9 <- fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(postfire.below ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)
# 
# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$postfire.below<-d_nb$postfire.below+.001
# d_nb$postfire.below<-round(d_nb$postfire.below,0)
# # d_nb$postfire.below<-as.integer(d_nb$postfire.below)
# 
# m11 <-  fitme(postfire.below ~ abmc,
#                           family = negbin1(),
#                           data = d_nb,
#                           resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(postfire.below ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)
# 
# m15<- fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude), 
#                          family = tweedie(var.power = 1.5, link.power = 0),  # var.power and link.power may need adjusting
#                          data = d, 
#                          resid.model = ~ abmc, 
#                          control = list(check.CGHQ = FALSE))
# 
# m16<-fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude), 
#            family = poisson(), 
#            data = d_nb, 
#            resid.model = ~ abmc)
# m16_AIC<-extractAIC(m16)
# 
# m17<-fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude), 
#            family = COMPoisson(), 
#            data = d_nb, 
#            resid.model = ~ abmc)
# m17_AIC<-extractAIC(m17)
# 
# m18<-fitme(postfire.below ~ abmc + Matern(1 | latitude + longitude), 
#            family = Tnegbin(), 
#            data = d2, 
#            resid.model = ~ abmc)
# m18_AIC<-extractAIC(m18)


###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m7)
testResiduals(m12)
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))

### Extract summary####
# Extract the response variable
winning_model<-m7
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/postfire.below_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(m7)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_postfire_below = mean(postfire.below, na.rm = TRUE),
    se_postfire_below = sd(postfire.below, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_postfire_below, se = se_postfire_below) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)



##COMB BELOW####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,m11_AIC,m12_AIC,m13_AIC,m14_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,d_nb,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,comb.below)%>%
  na.omit()

d1<-d%>%                 
  mutate(comb.below=comb.below+.1)


####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(comb.below ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(comb.below ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(comb.below ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(comb.below ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)

#Add a random effect of fire scar because spatial effect is not going great
m4 <- fitme(comb.below ~ abmc + Matern(1 | latitude + longitude) + (1|fire_scar),
            family = gaussian(),
            data = d,
            resid.model = ~ abmc)
m4_AIC<-extractAIC(m4) 

# m5 <- fitme(comb.below ~ abmc + Cauchy(1 | latitude + longitude),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m5_AIC<-extractAIC(m5)

# m6 <- fitme(comb.below ~ abmc + Cauchy(1 | latitude + longitude) + (1|fire_scar),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m6_AIC<-extractAIC(m6)


m7 <- fitme(comb.below ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

# m8 <- fitme(comb.below ~ abmc + Cauchy(1 | latitude + longitude),
#             family = Gamma(log),
#             data = d1,
#             resid.model = ~ abmc)
# m8_AIC<-extractAIC(m8)

m9 <- fitme(comb.below ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(comb.below ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)
# 
# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$comb.below<-d_nb$comb.below+.001
# d_nb$comb.below<-round(d_nb$comb.below,0)
# # d_nb$comb.below<-as.integer(d_nb$comb.below)
# 
# m11 <-  fitme(comb.below ~ abmc,
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(comb.below ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(comb.below ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(comb.below ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)

# m15 <- glmmTMB(
#   comb.below ~ abmc , # Fixed effect of abmc, random effect of fire_scar
#   family = tweedie(link = "log"),      # Tweedie family with log link
#   data = d,
#   dispformula = ~ abmc                # Optional: allow dispersion to vary by abmc
# )
# 
# # Start with a simpler model, e.g., Gaussian family
# m_simple <- glmmTMB(
#   comb.below ~ abmc + (1 | fire_scar),
#   family = gaussian(link = "log"),  # Start with Gaussian
#   data = d
# )
# 
# # Now fit the Tweedie model using starting values from the simpler model
# m_tweedie <- update(m_simple, family = tweedie(link = "log"))


###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m12)
testResiduals(m12)
testResiduals(m7)
testResiduals(m8)
testResiduals(m10)
testResiduals(m14)
testResiduals(m11)
testResiduals(m13)
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))


### Extract summary####
# Extract the response variable
winning_model<-m7
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.below_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(m7)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_comb_below = mean(comb.below, na.rm = TRUE),
    se_comb_below = sd(comb.below, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_comb_below, se = se_comb_below) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)


##PREFIRE BGT####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,m11_AIC,m12_AIC,m13_AIC,m14_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,d_nb,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,prefire.bgT)%>%
  na.omit()

d1<-d%>%                 
  mutate(prefire.bgT=prefire.bgT+.1)


####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(prefire.bgT ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(prefire.bgT ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(prefire.bgT ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(prefire.bgT ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)

#Add a random effect of fire scar because spatial effect is not going great
m4 <- fitme(prefire.bgT ~ abmc + Matern(1 | latitude + longitude) + (1|fire_scar),
            family = gaussian(),
            data = d,
            resid.model = ~ abmc)
m4_AIC<-extractAIC(m4) 

# m5 <- fitme(prefire.bgT ~ abmc + Cauchy(1 | latitude + longitude),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m5_AIC<-extractAIC(m5)
# 
# m6 <- fitme(prefire.bgT ~ abmc + Cauchy(1 | latitude + longitude) + (1|fire_scar),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m6_AIC<-extractAIC(m6)


m7 <- fitme(prefire.bgT ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

# m8 <- fitme(prefire.bgT ~ abmc + Cauchy(1 | latitude + longitude),
#             family = Gamma(log),
#             data = d1,
#             resid.model = ~ abmc)
# m8_AIC<-extractAIC(m8)

m9 <- fitme(prefire.bgT ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(prefire.bgT ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)
# 
# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$prefire.bgT<-d_nb$prefire.bgT+.001
# d_nb$prefire.bgT<-round(d_nb$prefire.bgT,0)
# # d_nb$prefire.bgT<-as.integer(d_nb$prefire.bgT)
# 
# m11 <-  fitme(prefire.bgT ~ abmc,
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(prefire.bgT ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(prefire.bgT ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(prefire.bgT ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)

# m15 <- glmmTMB(
#   prefire.bgT ~ abmc , # Fixed effect of abmc, random effect of fire_scar
#   family = tweedie(link = "log"),      # Tweedie family with log link
#   data = d,
#   dispformula = ~ abmc                # Optional: allow dispersion to vary by abmc
# )
# 
# # Start with a simpler model, e.g., Gaussian family
# m_simple <- glmmTMB(
#   prefire.bgT ~ abmc + (1 | fire_scar),
#   family = gaussian(link = "log"),  # Start with Gaussian
#   data = d
# )
# 
# # Now fit the Tweedie model using starting values from the simpler model
# m_tweedie <- update(m_simple, family = tweedie(link = "log"))


###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m7)
testResiduals(m4)
testResiduals(m11)
testResiduals(m13)
testResiduals(m7)
testResiduals(m8)
testResiduals(m11)
testResiduals(m13)
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))


### Extract summary####
# Extract the response variable
winning_model<-m7
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/prefire.bgT_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(m7)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_prefire_bgT = mean(prefire.bgT, na.rm = TRUE),
    se_prefire_bgT = sd(prefire.bgT, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_prefire_bgT, se = se_prefire_bgT) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)

##POSTFIRE BGT####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,m11_AIC,m12_AIC,m13_AIC,m14_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,d_nb,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,postfire.bgT)%>%
  na.omit()

d1<-d%>%                 
  mutate(postfire.bgT=postfire.bgT+.1)

####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(postfire.bgT ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(postfire.bgT ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(postfire.bgT ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(postfire.bgT ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)

#Add a random effect of fire scar because spatial effect is not going great
m4 <- fitme(postfire.bgT ~ abmc + Matern(1 | latitude + longitude) + (1|fire_scar),
            family = gaussian(),
            data = d,
            resid.model = ~ abmc)
m4_AIC<-extractAIC(m4) 

# m5 <- fitme(postfire.bgT ~ abmc + Cauchy(1 | latitude + longitude),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m5_AIC<-extractAIC(m5)
# 
# m6 <- fitme(postfire.bgT ~ abmc + Cauchy(1 | latitude + longitude) + (1|fire_scar),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m6_AIC<-extractAIC(m6)


m7 <- fitme(postfire.bgT ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

# m8 <- fitme(postfire.bgT ~ abmc + Cauchy(1 | latitude + longitude),
#             family = Gamma(log),
#             data = d1,
#             resid.model = ~ abmc)
# m8_AIC<-extractAIC(m8)

m9 <- fitme(postfire.bgT ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(postfire.bgT ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)
# 
# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$postfire.bgT<-d_nb$postfire.bgT+.001
# d_nb$postfire.bgT<-round(d_nb$postfire.bgT,0)
# # d_nb$postfire.bgT<-as.integer(d_nb$postfire.bgT)
# 
# m11 <-  fitme(postfire.bgT ~ abmc,
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(postfire.bgT ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(postfire.bgT ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(postfire.bgT ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)

# m15 <- glmmTMB(
#   postfire.bgT ~ abmc , # Fixed effect of abmc, random effect of fire_scar
#   family = tweedie(link = "log"),      # Tweedie family with log link
#   data = d,
#   dispformula = ~ abmc                # Optional: allow dispersion to vary by abmc
# )
# 
# # Start with a simpler model, e.g., Gaussian family
# m_simple <- glmmTMB(
#   postfire.bgT ~ abmc + (1 | fire_scar),
#   family = gaussian(link = "log"),  # Start with Gaussian
#   data = d
# )
# 
# # Now fit the Tweedie model using starting values from the simpler model
# m_tweedie <- update(m_simple, family = tweedie(link = "log"))



###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m12)
testResiduals(m12)
testResiduals(m11)
testResiduals(m13)
testResiduals(m7)
testResiduals(m8)
testResiduals(m11)
testResiduals(m13)
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))


###Outlier test####
testOutliers(m12)
residuals <- simulateResiduals(m12)
out_res<-residuals(residuals)

# Check for residuals beyond a certain threshold (depending on the distribution of residuals)
threshold <- quantile(out_res, probs = c(0.025, 0.975))  # for 95% interval
outliers <- which(out_res < threshold[1] | out_res > threshold[2])

# Get the actual values that are outliers
outlier_values <- d_nb[outliers, ]

d_no_outliers<-d_nb[-outliers,]

#Re-run model with no outliers
m12_out <-  fitme(postfire.bgT ~ abmc + Matern(1 | latitude + longitude),
                  family = negbin1(),
                  data = d_no_outliers,
                  resid.model = ~ abmc)
m12_out_AIC<-extractAIC(m12_out)


### Extract summary####
# Extract the response variable
winning_model<-m7
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/postfire.bgT_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(m7)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_postfire_bgT = mean(postfire.bgT, na.rm = TRUE),
    se_postfire_bgT = sd(postfire.bgT, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_postfire_bgT, se = se_postfire_bgT) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)

##COMB BGT####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,m11_AIC,m12_AIC,m13_AIC,m14_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,d_nb,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,comb.bgT)%>%
  na.omit()

d1<-d%>%                 
  mutate(comb.bgT=comb.bgT+.1)

####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(comb.bgT ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(comb.bgT ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(comb.bgT ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(comb.bgT ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)

#Add a random effect of fire scar because spatial effect is not going great
m4 <- fitme(comb.bgT ~ abmc + Matern(1 | latitude + longitude) + (1|fire_scar),
            family = gaussian(),
            data = d,
            resid.model = ~ abmc)
m4_AIC<-extractAIC(m4) 

# m5 <- fitme(comb.bgT ~ abmc + Cauchy(1 | latitude + longitude),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m5_AIC<-extractAIC(m5)
# 
# m6 <- fitme(comb.bgT ~ abmc + Cauchy(1 | latitude + longitude) + (1|fire_scar),
#             family = gaussian(), 
#             data = d, 
#             resid.model = ~ abmc)
# m6_AIC<-extractAIC(m6)


m7 <- fitme(comb.bgT ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

# m8 <- fitme(comb.bgT ~ abmc + Cauchy(1 | latitude + longitude),
#             family = Gamma(log),
#             data = d1,
#             resid.model = ~ abmc)
# m8_AIC<-extractAIC(m8)

m9 <- fitme(comb.bgT ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(comb.bgT ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)
# 
# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$comb.bgT<-d_nb$comb.bgT+.001
# d_nb$comb.bgT<-round(d_nb$comb.bgT,0)
# # d_nb$comb.bgT<-as.integer(d_nb$comb.bgT)
# 
# m11 <-  fitme(comb.bgT ~ abmc,
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(comb.bgT ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(comb.bgT ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(comb.bgT ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)

# m15 <- glmmTMB(
#   comb.bgT ~ abmc , # Fixed effect of abmc, random effect of fire_scar
#   family = tweedie(link = "log"),      # Tweedie family with log link
#   data = d,
#   dispformula = ~ abmc                # Optional: allow dispersion to vary by abmc
# )
# 
# # Start with a simpler model, e.g., Gaussian family
# m_simple <- glmmTMB(
#   comb.bgT ~ abmc + (1 | fire_scar),
#   family = gaussian(link = "log"),  # Start with Gaussian
#   data = d
# )
# 
# # Now fit the Tweedie model using starting values from the simpler model
# m_tweedie <- update(m_simple, family = tweedie(link = "log"))



###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m12)
testResiduals(m12)
testResiduals(m14)
testResiduals(m13)
testResiduals(m7)
testResiduals(m8)
testResiduals(m11)
testResiduals(m13)
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))

### Extract summary####
# Extract the response variable
winning_model<-m7
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.bgT_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(m7)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_comb_bgT = mean(comb.bgT, na.rm = TRUE),
    se_comb_bgT = sd(comb.bgT, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_comb_bgT, se = se_comb_bgT) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)


## % COMB TREES####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,m11_AIC,m12_AIC,m13_AIC,m14_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,d_nb,d_no_outliers,m12_out,outlier_values,m12_out_AIC,out_res,outliers,threshold,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,comb.trees.prop)%>%
  na.omit()

# d1<-d%>%                 
#   mutate(comb.trees.prop=comb.trees.prop+.00001) #No zero values for this one


####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(comb.trees.prop ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(comb.trees.prop ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(comb.trees.prop ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(comb.trees.prop ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)


m7 <- fitme(comb.trees.prop ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

m9 <- fitme(comb.trees.prop ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(comb.trees.prop ~ abmc ,
             family = Gamma(log),
             data = d,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)


# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$comb.trees.prop<-d_nb$comb.trees.prop+.001
# d_nb$comb.trees.prop<-100*(d_nb$comb.trees.prop)
# d_nb$comb.trees.prop<-round(d_nb$comb.trees.prop,0)
# d_nb$comb.trees.prop<-as.integer(d_nb$comb.trees.prop)

# m11 <-  fitme(comb.trees.prop ~ abmc,
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(comb.trees.prop ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(comb.trees.prop ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(comb.trees.prop ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)

#Beta resp and binomila familes for % data
m15 <-  fitme(comb.trees.prop ~ abmc,
              family = beta_resp(link=logit),
              data = d,
              resid.model = ~ abmc)
m15_AIC<-extractAIC(m15)

m16 <-  fitme(comb.trees.prop ~ abmc ,
              family = beta_resp(),
              data = d)
m16_AIC<-extractAIC(m16)


###CHECKS####
# Simulate residuals and check for spatial autocorrelation
# residuals <- simulateResiduals(m12)
testResiduals(m7)
plot(simulateResiduals(m7))
testResiduals(m15)
plot(simulateResiduals(m15))
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))

### Extract summary####
# Extract the response variable
winning_model<-m15
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.trees.prop_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(winning_model)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_comb_trees_prop = mean(comb.trees.prop, na.rm = TRUE),
    se_comb_trees_prop = sd(comb.trees.prop, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_comb_trees_prop, se = se_comb_trees_prop) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)


## % COMB BELOW####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,m11_AIC,m12_AIC,m13_AIC,m14_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,d_nb,d_no_outliers,m12_out,outlier_values,m12_out_AIC,out_res,outliers,threshold,m15,m16,m15_AIC,m16_AIC,d2,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,comb.below.prop)%>%
  na.omit()

d1<-d%>%
  mutate(comb.below.prop=comb.below.prop+.00001)


####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(comb.below.prop ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(comb.below.prop ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(comb.below.prop ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(comb.below.prop ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)


m7 <- fitme(comb.below.prop ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

m9 <- fitme(comb.below.prop ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d1) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(comb.below.prop ~ abmc ,
             family = Gamma(log),
             data = d1,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)

# 
# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$comb.below.prop<-d_nb$comb.below.prop+.001
# d_nb$comb.below.prop<-100*(d_nb$comb.below.prop)
# d_nb$comb.below.prop<-round(d_nb$comb.below.prop,0)
# # d_nb$comb.below.prop<-as.integer(d_nb$comb.below.prop)
# 
# m11 <-  fitme(comb.below.prop ~ abmc,
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(comb.below.prop ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(comb.below.prop ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(comb.below.prop ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)

d2<-d%>%
  mutate(comb.below.prop=ifelse(comb.below.prop==0,.00001,comb.below.prop))%>%
  mutate(comb.below.prop=ifelse(comb.below.prop==1,0.99999,comb.below.prop))
#Beta resp and binomila familes for % data
m15 <-  fitme(comb.below.prop ~ abmc,
              family = beta_resp(),
              data = d2,
              resid.model = ~ abmc)
m15_AIC<-extractAIC(m15)

m16 <-  fitme(comb.below.prop ~ abmc ,
              family = beta_resp(),
              data = d2)
m16_AIC<-extractAIC(m16)



###CHECKS####
# Simulate residuals and check for spatial autocorrelation
residuals <- simulateResiduals(m12)
testResiduals(m15)
plot(simulateResiduals(m15))
testResiduals(m7)
plot(simulateResiduals(m7))
testResiduals(m9)
plot(simulateResiduals(m9)) #passes QQ,  but not uniformity
testResiduals(m15)
plot(simulateResiduals(m15))
testResiduals(m12)
plot(simulateResiduals(m12)) #passes QQ and uniformity, but not leven test
testResiduals(m11)
plot(simulateResiduals(m11))
testResiduals(m13)
plot(simulateResiduals(m13))
testResiduals(m14)
plot(simulateResiduals(m14))
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))

### Extract summary####
# Extract the response variable
winning_model<-m15
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.below.prop_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(winning_model)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_comb_below_prop = mean(comb.below.prop, na.rm = TRUE),
    se_comb_below_prop = sd(comb.below.prop, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_comb_below_prop, se = se_comb_below_prop) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey)

## % COMB BGT####
rm(d,d1,m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,estimates_df,residuals,m0_AIC,m1_AIC,m2_AIC,m3_AIC,m4_AIC,m5_AIC,m6_AIC,m7_AIC,m8_AIC,m9_AIC,m10_AIC,m11_AIC,m12_AIC,m13_AIC,m14_AIC,new_data,combined_df,contrast_df,grouped_means_se,log_predictions_with_variances,predicted_means,response_scale_se,log_variances,d_nb,d_no_outliers,m12_out,outlier_values,m12_out_AIC,out_res,outliers,threshold,m15,m16,m15_AIC,m16_AIC,d2,df,df_contrasts)
d<-data%>%
  dplyr::select(site,latitude,longitude,fire_scar,abmc,comb.bgT.prop)%>%
  na.omit()

d1<-d%>%
  mutate(comb.bgT.prop=comb.bgT.prop+.00001)


####Models####
# Plain (no correlation structure, gaussian response distribution)
m0 <- fitme(comb.bgT.prop ~ abmc, 
            family = gaussian(), 
            data = d) 
m0_AIC<-extractAIC(m0)

#Spatial model with Gaussian random effects for spatial correlation
m1 <- fitme(comb.bgT.prop ~ abmc + Matern(1|latitude + longitude), 
            family = gaussian(), 
            data = d) 
m1_AIC<-extractAIC(m1)

# Fit a model with heteroscedasticity
m2 <- fitme(comb.bgT.prop ~ abmc, 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m2_AIC<-extractAIC(m2)


# Fit a Matern spatial correlation model with heteroscedasticity
m3 <- fitme(comb.bgT.prop ~ abmc + Matern(1 | latitude + longitude), 
            family = gaussian(), 
            data = d, 
            resid.model = ~ abmc)  # This accounts for heteroscedasticity
m3_AIC<-extractAIC(m3)


m7 <- fitme(comb.bgT.prop ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d,
            resid.model = ~ abmc) 
m7_AIC<-extractAIC(m7)

m9 <- fitme(comb.bgT.prop ~ abmc + Matern(1 | latitude + longitude),
            family = Gamma(log),
            data = d) 
m9_AIC<-extractAIC(m9)

m10 <- fitme(comb.bgT.prop ~ abmc ,
             family = Gamma(log),
             data = d,
             resid.model = ~ abmc)
m10_AIC<-extractAIC(m10)

# 
# # Fit a spaMM model using the Tweedie distribution
# d_nb<-d
# # d_nb$comb.bgT.prop<-d_nb$comb.bgT.prop+.001
# d_nb$comb.bgT.prop<-100*(d_nb$comb.bgT.prop)
# d_nb$comb.bgT.prop<-round(d_nb$comb.bgT.prop,0)
# # d_nb$comb.bgT.prop<-as.integer(d_nb$comb.bgT.prop)
# 
# m11 <-  fitme(comb.bgT.prop ~ abmc,
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m11_AIC<-extractAIC(m11)
# 
# m12 <-  fitme(comb.bgT.prop ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin1(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m12_AIC<-extractAIC(m12)
# 
# m13 <-  fitme(comb.bgT.prop ~ abmc,
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m13_AIC<-extractAIC(m13)
# 
# m14 <-  fitme(comb.bgT.prop ~ abmc + Matern(1 | latitude + longitude),
#               family = negbin2(),
#               data = d_nb,
#               resid.model = ~ abmc)
# m14_AIC<-extractAIC(m14)

d2<-d%>%
  mutate(comb.bgT.prop=ifelse(comb.bgT.prop==0,.00001,comb.bgT.prop))%>%
  mutate(comb.bgT.prop=ifelse(comb.bgT.prop==1,0.99999,comb.bgT.prop))
#Beta resp and binomila familes for % data
m15 <-  fitme(comb.bgT.prop ~ abmc,
              family = beta_resp(),
              data = d2,
              resid.model = ~ abmc)
m15_AIC<-extractAIC(m15)

m16 <-  fitme(comb.bgT.prop ~ abmc +Matern(1|longitude+latitude),
              family = beta_resp(),
              data = d2,
              resid.model = ~ abmc)
m16_AIC<-extractAIC(m16)

###CHECKS####
# Simulate residuals and check for spatial autocorrelation
testResiduals(m7)
plot(simulateResiduals(m7))
testResiduals(m15)
plot(simulateResiduals(m15))
plot(residuals)
testSpatialAutocorrelation(residuals, x = d$latitude, y = d$longitude)
1 - var(residuals$scaledResiduals) #"explained variance"
# check_model(m9, check=c("linearity", "normality", "vif", "qq", "homogeneity", "outliers", "ncv"))


### Extract summary####
# Extract the response variable
winning_model<-m15
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]

#Save model
# saveRDS(winning_model,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_model.rds"))
winning_model<-readRDS(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.bgT.prop_model.rds"))


AIC<-extractAIC(winning_model)[[2]]
response_variable <- attr(winning_model$main_terms_info$fixef_terms, "variables")[[2]]
df<- DoF(winning_model)["p_fixef"]  # Get DF for all fixed effects

# Step 1: Create new data for each abmc group (adjust for actual factor levels)
new_data <- data.frame(
  abmc = factor(c("conifer", "mixed", "aspen", "birch"), levels = levels(d$abmc)),
  latitude = mean(d$latitude),  # Use the mean latitude for prediction
  longitude = mean(d$longitude)  # Use the mean longitude for prediction
)

# Step 2: Get predictions on the response scale (biomass scale)
predicted_means <- predict(winning_model, newdata = new_data, re.form = NA, type = "response")

# Step 3: Get predictions on the log scale, along with variances
log_predictions_with_variances <- predict(winning_model, newdata = new_data, re.form = NA, type = "link", variances = list(predVar = TRUE))

# Extract the prediction variances (on the log scale)
log_variances <- attr(log_predictions_with_variances, "predVar")

# Step 4: Use the delta method to convert log scale variances to response scale standard errors
# Standard error on response scale = exp(predicted_mean) * sqrt(log_variances)
response_scale_se <- sqrt(log_variances) * predicted_means

# Step 5: Create a data frame with predictions and standard errors on the response scale
estimates_df <- data.frame(
  abmc = new_data$abmc,
  predicted_mean = predicted_means,
  standard_error = response_scale_se
)

# Print the response scale means and standard errors
estimates_df$pred.mean<-round(estimates_df$predicted_mean,0)
estimates_df$se<-round(estimates_df$standard_error,0)
print(estimates_df)

# Extract degrees of freedom from DoF()
df_contrasts <- DoF(winning_model)["p_fixef"]  # Extract DF for fixed effects

# Function to compute contrasts & adjusted p-values
compute_contrast_tukey <- function(mean1, se1, mean2, se2, df) {
  contrast <- mean1 - mean2
  combined_se <- sqrt(se1^2 + se2^2)
  t_value <- contrast / combined_se
  
  # Compute p-values
  p_value <- 2 * pt(-abs(t_value), df = Inf)  # Inf-based p-value
  p_value_df <- 2 * pt(-abs(t_value), df = df)  # DoF-based p-value
  
  return(list(contrast = contrast, t_value = t_value, p_value = p_value, p_value_df = p_value_df))
}

# Compute pairwise contrasts
pairwise_contrasts <- data.frame(
  contrast = c("mixed vs conifer", "aspen vs conifer", "birch vs conifer", 
               "aspen vs mixed", "birch vs mixed", "birch vs aspen"),
  
  contrast_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$contrast,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$contrast
  ),
  
  t_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$t_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$t_value
  ),
  
  p_value = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value
  ),
  
  p_value_df = c(
    compute_contrast_tukey(predicted_means[2], response_scale_se[2], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[1], response_scale_se[1], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[3], response_scale_se[3], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[2], response_scale_se[2], df_contrasts)$p_value_df,
    compute_contrast_tukey(predicted_means[4], response_scale_se[4], predicted_means[3], response_scale_se[3], df_contrasts)$p_value_df
  )
)

# Apply Tukey-Kramer adjustment
pairwise_contrasts$p_adjusted_inf <- p.adjust(pairwise_contrasts$p_value, method = "holm")  
pairwise_contrasts$p_adjusted_df <- p.adjust(pairwise_contrasts$p_value_df, method = "holm")  

# Round values for better readability
pairwise_contrasts$pval_inf <- round(pairwise_contrasts$p_value, 3)
pairwise_contrasts$pval_df <- round(pairwise_contrasts$p_value_df, 3)
pairwise_contrasts$pval_adj_inf <- round(pairwise_contrasts$p_adjusted_inf, 3)
pairwise_contrasts$pval_adj_df <- round(pairwise_contrasts$p_adjusted_df, 3)

# Print final results
print(pairwise_contrasts)

####Make figure####
# Assuming your dataframe is named 'd'
grouped_means_se <- d %>%
  group_by(abmc) %>%
  summarise(
    mean_comb_bgT_prop = mean(comb.bgT.prop, na.rm = TRUE),
    se_comb_bgT_prop = sd(comb.bgT.prop, na.rm = TRUE) / sqrt(n())
  )

# Create a combined data frame for raw and modeled means and SEs
combined_df <- bind_rows(
  grouped_means_se %>%
    mutate(source = "Raw", mean = mean_comb_bgT_prop, se = se_comb_bgT_prop) %>%
    dplyr::select(abmc, source, mean, se),
  
  estimates_df %>%
    mutate(source = "Modeled", mean = predicted_mean, se = standard_error) %>%
    dplyr::select(abmc, source, mean, se)
)

####Write out results####
write.csv(estimates_df,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_estimated_means.csv"),row.names=F)
write.csv(pairwise_contrasts,paste0(data_directory,"/model_outputs/Q1_ABMC_models/",response_variable,"_contrasts.csv"),row.names=F)

# Create the plot
g<-ggplot(combined_df, aes(x = abmc, y = mean, fill = abmc,alpha=as.factor(source))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(y = response_variable, x = "ABMC Group", fill = "Source") +
  theme_bw() +
  ggtitle(paste0("Raw vs Modeled ",response_variable," by ABMC Group"))+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_data.png"), plot = g, width = 10, height = 6, dpi = 300)

#Extract R2
r2<-pseudoR2(winning_model)

# Extract the main distribution family and link function
distribution_family <- (winning_model$family$family)[[1]]
link_function <- winning_model$family$link

# Check for spatial structure
spatial_structure <- any(attr(winning_model$spatial_terms, "type") == "Matern")

#Extract model residuals
residuals<-simulateResiduals(winning_model)

#Test global uniformity, dispersion, and outliers
checks<-testResiduals(winning_model)
KStest_significant <- ifelse(checks$uniformity$p.value < 0.05, "Yes", "No")
dispersion_significant <- ifelse(checks$dispersion$p.value < 0.05, "Yes", "No")
outlier_significant <- ifelse(checks$outliers$p.value < 0.05, "Yes", "No")
remove_outlier_change<-"no" #Adjust for outlier chanegs

#Test within group uniformity
group_uniformity<-testCategorical(residuals,d$abmc)
levene_significant<-ifelse((group_uniformity$homogeneity$`Pr(>F)`)[[1]]< 0.05, "Yes", "No")
group_uniformity_pvalues<-group_uniformity$uniformity$p.value
conifer_uniformity_significant<-ifelse((group_uniformity_pvalues)[[1]]< 0.05, "Yes", "No")
mixed_uniformity_significant<-ifelse((group_uniformity_pvalues)[[2]]< 0.05, "Yes", "No")
aspen_uniformity_significant<-ifelse((group_uniformity_pvalues)[[3]]< 0.05, "Yes", "No")
birch_uniformity_significant<-ifelse((group_uniformity_pvalues)[[4]]< 0.05, "Yes", "No")

# Extract fixed effects from the summary
model_summary <- summary(winning_model)
fixed_effects <- as.data.frame(model_summary$beta_table)
fixed_effects$df <- DoF(winning_model) 

# Calculate p-values for each fixed effect
fixed_effects$p_value <- 2 * pnorm(-abs(fixed_effects$`t-value`))

# Adjust p-values for multiple comparisons (Holm's method recommended)
fixed_effects$p_value_adjusted <- p.adjust(fixed_effects$p_value, method = "holm")

# Compute p-values using the extracted DF
fixed_effects$p_value_df <- 2 * pt(-abs(fixed_effects$`t-value`), df = fixed_effects$df)

# Adjust p-values using Holm’s method
fixed_effects$p_value_adjusted_df <- p.adjust(fixed_effects$p_value_df, method = "holm")

# Calculate 95% Confidence Intervals
fixed_effects$CI_lower <- fixed_effects$Estimate - 1.96 * fixed_effects$`Cond. SE`
fixed_effects$CI_upper <- fixed_effects$Estimate + 1.96 * fixed_effects$`Cond. SE`

# Rename rows for clarity
rownames(fixed_effects) <- c("conifer", "mixed", "aspen", "birch")

# Create a single-row dataframe with appropriate column names
output_data <- data.frame(
  Response_Variable = as.character(response_variable),
  Distribution_Family = distribution_family,
  Link_Function = link_function,
  Spatial_Structure = spatial_structure,
  R2=r2,
  AIC=AIC,
  KStest_significant = KStest_significant,
  dispersion_significant = dispersion_significant,
  outlier_significant = outlier_significant,
  remove_outlier_change = remove_outlier_change,
  levene_significant = levene_significant,
  conifer_uniformity_significant = conifer_uniformity_significant,
  mixed_uniformity_significant = mixed_uniformity_significant,
  aspen_uniformity_significant = aspen_uniformity_significant,
  birch_uniformity_significant = birch_uniformity_significant,
  # Fixed effects estimates, SEs, t-values, and p-values
  
  # Fixed effects estimates, SEs, t-values, p-values, and confidence intervals
  conifer_Estimate = fixed_effects["conifer", "Estimate"],
  conifer_SE = fixed_effects["conifer", "Cond. SE"],
  conifer_tvalue = fixed_effects["conifer", "t-value"],
  conifer_pvalue = round(fixed_effects["conifer", "p_value"], 3),
  conifer_pvalue_adj = round(fixed_effects["conifer", "p_value_adjusted"], 3),
  conifer_pvalue_df = round(fixed_effects["conifer", "p_value_df"], 3),
  conifer_pvalue_adj_df = round(fixed_effects["conifer", "p_value_adjusted_df"], 3),
  conifer_CI_lower = round(fixed_effects["conifer", "CI_lower"], 4),
  conifer_CI_upper = round(fixed_effects["conifer", "CI_upper"], 4),
  
  mixed_Estimate = fixed_effects["mixed", "Estimate"],
  mixed_SE = fixed_effects["mixed", "Cond. SE"],
  mixed_tvalue = fixed_effects["mixed", "t-value"],
  mixed_pvalue = round(fixed_effects["mixed", "p_value"], 3),
  mixed_pvalue_adj = round(fixed_effects["mixed", "p_value_adjusted"], 3),
  mixed_pvalue_df = round(fixed_effects["mixed", "p_value_df"], 3),
  mixed_pvalue_adj_df = round(fixed_effects["mixed", "p_value_adjusted_df"], 3),
  mixed_CI_lower = round(fixed_effects["mixed", "CI_lower"], 4),
  mixed_CI_upper = round(fixed_effects["mixed", "CI_upper"], 4),
  
  aspen_Estimate = fixed_effects["aspen", "Estimate"],
  aspen_SE = fixed_effects["aspen", "Cond. SE"],
  aspen_tvalue = fixed_effects["aspen", "t-value"],
  aspen_pvalue = round(fixed_effects["aspen", "p_value"],  3),
  aspen_pvalue_adj = round(fixed_effects["aspen", "p_value_adjusted"], 3),
  aspen_pvalue_df = round(fixed_effects["aspen", "p_value_df"], 3),
  aspen_pvalue_adj_df = round(fixed_effects["aspen", "p_value_adjusted_df"], 3),
  aspen_CI_lower = round(fixed_effects["aspen", "CI_lower"], 4),
  aspen_CI_upper = round(fixed_effects["aspen", "CI_upper"], 4),
  
  birch_Estimate = fixed_effects["birch", "Estimate"],
  birch_SE = fixed_effects["birch", "Cond. SE"],
  birch_tvalue = fixed_effects["birch", "t-value"],
  birch_pvalue = round(fixed_effects["birch", "p_value"],   3),
  birch_pvalue_adj = round(fixed_effects["birch", "p_value_adjusted"], 3),
  birch_pvalue_df = round(fixed_effects["birch", "p_value_df"], 3),
  birch_pvalue_adj_df = round(fixed_effects["birch", "p_value_adjusted_df"], 3),
  birch_CI_lower = round(fixed_effects["birch", "CI_lower"], 4),
  birch_CI_upper = round(fixed_effects["birch", "CI_upper"], 4)
)

#Print plots
residuals <- simulateResiduals(winning_model)

png(paste0(data_directory, "/model_outputs/Q1_ABMC_models/", response_variable, "_residuals.png"), 
    width = 10, height = 6, units = "in", res = 300)
plot(residuals)
dev.off()

full_output_data<-bind_rows(full_output_data,output_data)
rm(checks,fixed_effects,g,model_summary,output_data,winning_model,dispersion_significant,distribution_family,link_function,outlier_significant,remove_outlier_change,response_variable,spatial_structure,conifer_uniformity_significant,mixed_uniformity_significant,aspen_uniformity_significant,birch_uniformity_significant,group_uniformity,group_uniformity_pvalues,KStest_significant,levene_significant,r2,AIC,pairwise_contrasts,compute_contrast_tukey,df,df_contrasts)


#WRITE OUT MODEL STATS####
write.csv(full_output_data,paste0(data_directory,"/model_outputs/Q1_ABMC_models/full_modelstats.csv"),row.names=F)



#ABSOLUTE POOL FIGURES####

##CALL MODEL OUTPUTS FOR FIGURES####
x1<-read.csv(paste0(data_directory,"/model_outputs/Q1_ABMC_models/prefire.trees_estimated_means.csv"))
x1<-x1%>%
  dplyr::rename(SE=standard_error,
                mean=predicted_mean)%>%
  dplyr::select(abmc,SE,mean)
x1$prepost<-"Pre-fire C"
x1$abovebelow<-"Aboveground"
x2<-read.csv(paste0(data_directory,"/model_outputs/Q1_ABMC_models/postfire.trees_estimated_means.csv"))
x2<-x2%>%
  dplyr::rename(SE=standard_error,
                mean=predicted_mean)%>%
  dplyr::select(abmc,SE,mean)
x2$prepost<-"Post-fire C"
x2$abovebelow<-"Aboveground"
x3<-read.csv(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.trees_estimated_means.csv"))
x3<-x3%>%
  dplyr::rename(SE=standard_error,
                mean=predicted_mean)%>%
  dplyr::select(abmc,SE,mean)
x3$prepost<-"C Loss"
x3$abovebelow<-"Aboveground"
x4<-read.csv(paste0(data_directory,"/model_outputs/Q1_ABMC_models/prefire.below_estimated_means.csv"))
x4<-x4%>%
  dplyr::rename(SE=standard_error,
                mean=predicted_mean)%>%
  dplyr::select(abmc,SE,mean)
x4$prepost<-"Pre-fire C"
x4$abovebelow<-"Belowground"
x5<-read.csv(paste0(data_directory,"/model_outputs/Q1_ABMC_models/postfire.below_estimated_means.csv"))
x5<-x5%>%
  dplyr::rename(SE=standard_error,
                mean=predicted_mean)%>%
  dplyr::select(abmc,SE,mean)
x5$prepost<-"Post-fire C"
x5$abovebelow<-"Belowground"
x6<-read.csv(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.below_estimated_means.csv"))
x6<-x6%>%
  dplyr::rename(SE=standard_error,
                mean=predicted_mean)%>%
  dplyr::select(abmc,SE,mean)
x6$prepost<-"C Loss"
x6$abovebelow<-"Belowground"

datayc1<-bind_rows(x1,x2)
datayc1<-bind_rows(datayc1,x3)
datayc1<-bind_rows(datayc1,x4)
datayc1<-bind_rows(datayc1,x5)
datayc1<-bind_rows(datayc1,x6)

datayc1<-datayc1%>%
  mutate(abmc = str_to_title(abmc))

dyc<-datayc1%>%
  group_by(abmc,prepost)%>%
  summarize(summean=sum(mean))

datayc1<-full_join(datayc1,dyc)
datayc1<-datayc1%>%
  mutate(summean=ifelse(abovebelow=="Belowground",mean,summean))
rm(dyc)

datayc1$abmc<-factor(datayc1$abmc, levels=c("Conifer","Mixed", "Aspen", "Birch"))
datayc1$prepost<-factor(datayc1$prepost, levels=c("Pre-fire C","Post-fire C", "C Loss"))
datayc1$abovebelow<-factor(datayc1$abovebelow, levels=c("Aboveground","Belowground"))

# # Un-altered p-value
# datayc1$above_letter <- c("a", "b", "b", "c", "a", "b", "b", "c", "a", "a", "b", "b", "","","","","","","","","","","","")
#adjusted pvalue
datayc1$above_letter <- c("a", "b", "b", "b", "a", "b", "b", "b", "a", "a", "b", "b", "","","","","","","","","","","","")
# #adjusted and df pvalue
# datayc1$above_letter <- c("a", "b", "b", "b", "a", "b", "b", "b", "a", "a", "a", "a", "","","","","","","","","","","","")


# # Un-altered p-value                          
# datayc1$below_letter <- c("","","","","","","","","","","","","a", "b", "c", "b", "a", "bc", "c", "ab", "a", "ab", "c", "bc")
#
#adjusted pvalue
datayc1$below_letter <- c("","","","","","","","","","","","",
                          "a", "b", "c", "b", "a", "ab", "b", "a", "a", "ab", "b", "b")

# #adjusted and df pvalue
# datayc1$below_letter <- c("","","","","","","","","","","","",
#                           "a", "a", "b", "ab", "a", "a", "a", "a", "a", "a", "a", "a")
# Assign colors for below_letter
datayc1$below_letter_color <-  c(  "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black","white", "black", "black", "black", "white", "black", "black", "black", "white", "black", "black")
# 
# # Un-altered p-value    
# datayc1$top_letter <- c("A", "B", "A", "C", "A", "B", "B", "C", "A", "A", "B", "B", "","","","","","","","","","","","")
#adjusted pvalue
datayc1$top_letter <- c("A", "B", "A", "B", "A", "B", "B", "C", "A", "A", "B", "B", "","","","","","","","","","","","")
# #adjusted and df pvalue
# datayc1$top_letter <- c("A", "AB", "AB", "B", "A", "AB", "AB", "B", "A", "AB", "B", "B", "","","","","","","","","","","","")

datayc1<-datayc1%>%
  mutate(mean=mean/1000)%>%
  mutate(summean=summean/1000)%>%
  mutate(SE=SE/1000)

##MODELED MEANS PLOT####
b=ggplot()+
  geom_bar(data=datayc1,aes(x=datayc1$abmc, y=mean,fill=abmc,alpha=as.factor(abovebelow)),color="black",stat="summary",size=.3)+
  facet_wrap(~prepost)+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  geom_errorbar(data=datayc1,aes(x=datayc1$abmc,ymin=summean-SE, ymax=summean+SE,y=summean),fill=datayc1$abmc,alpha=as.factor(datayc1$abovebelow),width=0.25, size=0.2)+
  # Add letters above the SE bars
  geom_text(data = datayc1,
            aes(x = abmc, y = summean + SE + .300, label = above_letter),
            size = 4, vjust = 0) +
  # Add letters below the SE bars
  geom_text(data = datayc1,
            aes(x = abmc, y = summean - SE - .100, label = below_letter, color = below_letter_color),
            size = 4, vjust = 1)+
  scale_color_identity()+# Use colors as specified in the data
  # Add letters above the SE bars
  geom_text(data = datayc1,
            aes(x = abmc, y = 12.200, label = top_letter),
            size = 4, vjust = 0)+
  ylim(0,12.550)
b
b1=b+theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text( size=10, vjust=1.5), axis.text.x  = element_text(size=10)) +
  theme(axis.title.y = element_text( size=13, vjust=1.5), axis.text.y  = element_text(size=9))+
  ylab(expression('Carbon (kg C m'^-2*')'))+ xlab("")+
  # theme(legend.position = c(1, .8), legend.justification = c(1, 1))
  theme(legend.position="top", legend.spacing.x = unit(0.1, 'cm'))
b1
b2=b1+theme(legend.background=element_blank()) + guides(fill="none")
b2
#g3=g2+scale_y_continuous(breaks=number_ticks(8), limits=c(0,9100))
b3=b2+theme(strip.text.x = element_text(size = 12)) +theme(strip.background = element_rect(fill="white", size=0.18))+
  theme(legend.text=element_text(size=12))+
  theme(axis.title.x = element_blank())+ theme(panel.border = element_rect(colour = "black", fill=NA, size=0.18))+
  theme(axis.ticks.length=unit(.05, "cm"), axis.ticks = element_line(colour = "black", size = 0.18)) +
  theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.18))
b3
b4=b3+
  theme(
    strip.text.x = element_text(margin = ggplot2::margin(0.05,0,0.05,0, "cm")),
    plot.margin = unit(c(0.15,0.15,0.05,0), "cm"),
    axis.text.y = element_text(hjust = 1), 
    axis.title.y = element_text(vjust = 0),
    legend.box.margin=ggplot2::margin(0,-1,-5,-1),
    legend.margin=ggplot2::margin(0.15,0,-0.15,0),
    legend.key.width = unit(0.3, 'cm'), 
    legend.key.height = unit(0.0, 'cm'),
    legend.position = "top",                  # Keep legend above the graph
    legend.justification = "center",  #c(1, 0),           # Align legend to the right (top-right)
    # legend.box.just = "right",                 # Push legend further to the right
    # axis.text.x = element_blank(),             # Remove x-axis text
    # axis.text.x =element_text(size = 13),         #JK add it back in
    axis.title.x = element_blank())

b4
rm(b,b1,b2,b3)

b_final <- b4 +
  theme(
    text = element_text(family = "Arial", size = 6.5),
    axis.text.x = element_text(size = 11) ,
    axis.text.y = element_text(size = 8),
    # axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    # legend.title = element_text(size = 6),
    strip.text.x = element_text(size = 12)
  )
b_final
fig1=b_final

ggsave(
  filename = "C:/Users/Test/Documents/MS THESIS/PaperDrafts/FinalFigs/Figure2.pdf",
  plot = b_final,
  width = 180, height = 120, units = "mm",
  device = cairo_pdf, # Ensures embedded fonts and true vector output
  dpi = 300
)

# ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/absolute_pools_full_figure_modeled.png"), plot = b4, width = 7.5, height = 5, dpi = 300)
ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/absolute_pools_full_figure_modeled_pval_adj.png"), plot = b4, width = 6.5, height = 4.5, dpi = 300)
# ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/absolute_pools_full_figure_modeled_pval_adj_df.png"), plot = b4, width = 6.5, height = 4.5, dpi = 300)

rm(b4,datayc1,x1,x2,x3,x4,x5,x6)

##RAW MEANS PLOT####
library(dplyr)
library(tidyr)
library(plotrix)
library(ggplot2)
names(data)
dataxc=dplyr::select(data,site,plot,abmc,prefire.below,comb.below,postfire.below,prefire.trees,comb.trees,postfire.trees)

dataxc=reshape2::melt(dataxc,by=c("site","plot","abmc"))

dataxc=dataxc%>%
  mutate(prepost=ifelse(variable=="prefire.trees"|variable=="prefire.below", "Pre-fire C",
                        ifelse(variable=="postfire.trees"|variable=="postfire.below", "Post-fire C","C Loss")))
dataxc=dataxc%>%
  mutate(abovebelow=ifelse(variable=="prefire.below" | variable=="postfire.below"|variable=="comb.below", "Belowground", "Aboveground"))

dataxc=subset(dataxc,!is.na(value))
datayc=dataxc%>%
  group_by(abmc,prepost, abovebelow)%>%
  summarize(SE=std.error(value), mean=mean(value))

dataxc$abmc<-factor(dataxc$abmc,levels=c("conifer","mixed","aspen","birch"))
dataxc$prepost<-factor(dataxc$prepost,levels=c("Pre-fire C","Post-fire C","C Loss"))


dataxc=subset(dataxc,!is.na(abmc))
datayc=subset(datayc,!is.na(abmc))
dataxc=subset(dataxc,!is.na(value))
dyc=datayc%>%
  group_by(abmc,prepost)%>%
  summarise(summean=sum(mean))
datayc=full_join(datayc,dyc)

datayc=datayc%>%
  mutate(mean=ifelse(abovebelow=="Aboveground",summean,mean))
#
dataxc$prepost=factor(dataxc$prepost,levels=c('Pre-fire C','Post-fire C','C Loss'))
datayc$prepost=factor(datayc$prepost,levels=c('Pre-fire C','Post-fire C','C Loss'))


rm(b)
b=ggplot()+
  geom_bar(data=dataxc,aes(x=dataxc$abmc, y=value,fill=abmc,alpha=as.factor(abovebelow)),color="black",stat="summary")+
  facet_wrap(~prepost)+
  scale_alpha_discrete(range=c(.4,1))+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  geom_errorbar(data=datayc,aes(x=datayc$abmc,ymin=mean-SE, ymax=mean+SE,y=mean),fill=datayc$abmc,alpha=as.factor(datayc$abovebelow),width=0.25, size=0.2)
b
b1=b+theme(legend.title=element_blank()) + theme(legend.position="top", legend.spacing.x = unit(0.1, 'cm'))+
  theme(axis.title.x = element_text( size=10, vjust=1.5), axis.text.x  = element_text(size=14)) +
  theme(axis.title.y = element_text( size=14, vjust=1.5), axis.text.y  = element_text(size=14))+
  ylab(expression('Carbon (g C m'^-2*')'))+ xlab("")
b1
b2=b1+theme(legend.background=element_blank())+ guides(fill="none")
b2
#g3=g2+scale_y_continuous(breaks=number_ticks(8), limits=c(0,9100))
b3=b2+theme(strip.text.x = element_text(size = 14)) +theme(strip.background = element_rect(fill="white", size=0.18))+
  theme(legend.text=element_text(size=14))+
  theme(axis.title.x = element_blank())+ theme(panel.border = element_rect(colour = "black", fill=NA, size=0.18))+
  theme(axis.ticks.length=unit(.05, "cm"), axis.ticks = element_line(colour = "black", size = 0.18)) +
  theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.18))
b3
b4=b3+
  theme(strip.text.x = element_text(margin = ggplot2::margin(0.05,0,0.05,0, "cm"))) + theme(plot.margin = unit(c(0.15,0.15,0.05,0), "cm"))+
  theme(axis.text.y = element_text(hjust = 1), axis.title.y = element_text(vjust = 0)) +
  theme(legend.box.margin=ggplot2::margin(0,-1,-5,-1))+
  theme(legend.margin=ggplot2::margin(0.15,0,-0.15,0), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.0, 'cm'))
b4
rm(b,b1,b2,b3)

#PROPORTIONAL POOL FIGURES####

##CALL MODEL OUTPUTS FOR FIGURES####
x1<-read.csv(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.trees.prop_estimated_means.csv"))
x1 <- x1 %>%
  mutate(across(where(is.numeric), ~ . * 100)) 
x1<-x1%>%
  dplyr::rename(SE=standard_error,
                mean=predicted_mean)%>%
  dplyr::select(abmc,SE,mean)
# x1$prepost<-"Pre-fire C"
x1$abovebelow<-"Aboveground"
x2<-read.csv(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.below.prop_estimated_means.csv"))
x2 <- x2 %>%
  mutate(across(where(is.numeric), ~ . * 100)) 
x2<-x2%>%
  dplyr::rename(SE=standard_error,
                mean=predicted_mean)%>%
  dplyr::select(abmc,SE,mean)
# x2$prepost<-"Post-fire C"
x2$abovebelow<-"Belowground"
x3<-read.csv(paste0(data_directory,"/model_outputs/Q1_ABMC_models/comb.bgT.prop_estimated_means.csv"))
x3 <- x3 %>%
  mutate(across(where(is.numeric), ~ . * 100))
x3<-x3%>%
  dplyr::rename(SE=standard_error,
                mean=predicted_mean)%>%
  dplyr::select(abmc,SE,mean)
# x3$prepost<-"C Loss"
x3$abovebelow<-"Total"

datayc1<-bind_rows(x1,x2)
datayc1<-bind_rows(datayc1,x3)

datayc1<-datayc1%>%
  mutate(abmc = str_to_title(abmc))

datayc1$abmc<-factor(datayc1$abmc, levels=c("Conifer","Mixed", "Aspen", "Birch"))
# datayc1$prepost<-factor(datayc1$prepost, levels=c("Pre-fire C","Post-fire C", "C Loss"))
datayc1$abovebelow<-factor(datayc1$abovebelow, levels=c("Total","Aboveground","Belowground"))

# #raw pval
# datayc1$above_letter <- c("a", "b", "c", "d", "ab", "ab", "a", "b", "a", "a", "b", "b")
#pval adjusted
datayc1$above_letter <- c("a", "b", "c", "d", "a", "a", "a", "a", "a", "a", "b", "b")

#Add fraction below

datayc1$frac_below <- c(0,0,0,0, 1,1,1,1, 0.89, 0.89, 0.88, 0.87)

##MODELED MEANS PLOT####
b <- ggplot() +
  # Bar plot of mean proportional carbon loss by abmc group
  geom_bar(data = datayc1, aes(x = abmc, y = mean*frac_below, fill = abmc), 
           color = NA, stat = "summary",size=.3,alpha=1) +
  geom_bar(data = datayc1, aes(x = abmc, y = mean, fill = abmc), 
           color = "black", stat = "summary",size=.3,alpha=.4) +
  # Bar plot of mean proportional carbon loss by abmc group
  # Facet by abovebelow, arranging panels vertically
  facet_wrap(~abovebelow, ncol = 1) +
  # Custom fill colors for abmc groups
  scale_fill_manual(values = c("cornflowerblue", "indianred4", "goldenrod1", "chocolate2")) +
  theme_bw() +
  # Error bars showing mean ± SE
  geom_errorbar(data = datayc1, 
                aes(x = abmc, ymin = mean - SE, ymax = mean + SE), 
                width = 0.25, size = 0.2) +
  # Annotate bars with significance letters above SE bars
  geom_text(data = datayc1, 
            aes(x = abmc, y = mean + SE + 10, label = above_letter), 
            size = 4, vjust = 0) +
  ylim(0, 100)
b
# Improve legend positioning and axis text sizes
b1 <- b + theme(axis.title.x = element_blank(), 
                axis.text.x = element_text(size = 11),
                axis.title.y = element_text(size = 13, vjust = 1.5), 
                axis.text.y = element_text(size = 9)) +
  ylab(expression('Proportional Carbon Loss (%)'))
b1
# Remove legend fill
b2 <- b1 + theme(legend.background = element_blank()) + guides(fill = "none")
b2
# Adjust panel and text formatting
b3 <- b2 + theme(strip.text.x = element_text(size = 12), 
                 strip.background = element_rect(fill = "white", size = 0.18),
                 # legend.text = element_text(size = 14),
                 panel.border = element_rect(colour = "black", fill = NA, size = 0.18),
                 axis.ticks.length = unit(.05, "cm"), 
                 axis.ticks = element_line(colour = "black", size = 0.18),
                 panel.grid.minor = element_line(size = 0.1), 
                 panel.grid.major = element_line(size = 0.18))
b3
# Final layout tweaks for margins and alignment
b4 <- b3 + theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")),
                 plot.margin = unit(c(0.15, 0.15, 0.05, 0), "cm"),
                 axis.text.y = element_text(hjust = 1), 
                 axis.title.y = element_text(vjust = 0),
                 legend.box.margin = margin(0, -1, -5, -1),
                 legend.margin = margin(0.15, 0, -0.15, 0), 
                 legend.key.width = unit(0.3, 'cm'), 
                 legend.key.height = unit(0.0, 'cm'))

b4

b_final <- b4 +
  theme(
    text = element_text(family = "Arial", size = 6.5),
    axis.text.x = element_text(size = 12) ,
    axis.text.y = element_text(size = 10),
    # axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    # legend.title = element_text(size = 6),
    strip.text.x = element_text(size = 12,margin=margin(t=3,r=3,b=3,l=3,"pt"))
  )
b_final
fig2=b_final
ggsave(
  filename = "C:/Users/Test/Documents/MS THESIS/PaperDrafts/FinalFigs/Figure3.pdf",
  plot = b_final,
  width = 88, height = 180, units = "mm",
  device = cairo_pdf, # Ensures embedded fonts and true vector output
  dpi = 300
)

# 
# ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/proportional_pools_full_figure_modeled.png"), plot = b4, width = 2.5, height = 4.5, dpi = 300,units="in")

ggsave(paste0(data_directory, "/model_outputs/Q1_ABMC_models/proportional_pools_full_figure_modeled_pval_adj.png"), plot = b4, width = 2.5, height = 4.5, dpi = 300,units="in")

##RAW MEANS PLOT####
library(dplyr)
library(tidyr)
library(plotrix)
library(ggplot2)
names(data)
dataxc=dplyr::select(data,site,plot,abmc,comb.trees.prop,comb.below.prop,comb.bgT.prop)

dataxc=reshape2::melt(dataxc,by=c("site","plot","abmc"))
# 
# dataxc=dataxc%>%
#   mutate(prepost=ifelse(variable=="prefire.trees"|variable=="prefire.below", "Pre-fire C",
#                         ifelse(variable=="postfire.trees"|variable=="postfire.below", "Post-fire C","C Loss")))
dataxc=dataxc%>%
  mutate(abovebelow=ifelse(variable=="comb.trees.prop", "Aboveground",
                           ifelse(variable=="comb.below.prop", "Belowground", "Total")))

dataxc=subset(dataxc,!is.na(value))
datayc=dataxc%>%
  group_by(abmc,abovebelow)%>%
  summarize(SE=std.error(value), mean=mean(value))

dataxc$abmc<-factor(dataxc$abmc,levels=c("conifer","mixed","aspen","birch"))
# dataxc$prepost<-factor(dataxc$prepost,levels=c("Pre-fire C","Post-fire C","C Loss"))


dataxc=subset(dataxc,!is.na(abmc))
datayc=subset(datayc,!is.na(abmc))
dataxc=subset(dataxc,!is.na(value))

rm(b)
b=ggplot()+
  geom_bar(data=dataxc,aes(x=dataxc$abmc, y=value,fill=abmc),color="black",stat="summary")+
  facet_wrap(~abovebelow)+
  scale_fill_manual(values=c("cornflowerblue" ,"indianred4","goldenrod1","chocolate2"))+
  theme_bw()+
  geom_errorbar(data=datayc,aes(x=datayc$abmc,ymin=mean-SE, ymax=mean+SE,y=mean),fill=datayc$abmc,width=0.25, size=0.2)
b
b1=b+theme(legend.title=element_blank()) + theme(legend.position="top", legend.spacing.x = unit(0.1, 'cm'))+
  theme(axis.title.x = element_text( size=10, vjust=1.5), axis.text.x  = element_text(size=14)) +
  theme(axis.title.y = element_text( size=14, vjust=1.5), axis.text.y  = element_text(size=14))+
  ylab(expression('Proportional Combustion (%)'))+ xlab("")
b1
b2=b1+theme(legend.background=element_blank())+ guides(fill="none")
b2
#g3=g2+scale_y_continuous(breaks=number_ticks(8), limits=c(0,9100))
b3=b2+theme(strip.text.x = element_text(size = 14)) +theme(strip.background = element_rect(fill="white", size=0.18))+
  theme(legend.text=element_text(size=14))+
  theme(axis.title.x = element_blank())+ theme(panel.border = element_rect(colour = "black", fill=NA, size=0.18))+
  theme(axis.ticks.length=unit(.05, "cm"), axis.ticks = element_line(colour = "black", size = 0.18)) +
  theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.18))
b3
b4=b3+
  theme(strip.text.x = element_text(margin = ggplot2::margin(0.05,0,0.05,0, "cm"))) + theme(plot.margin = unit(c(0.15,0.15,0.05,0), "cm"))+
  theme(axis.text.y = element_text(hjust = 1), axis.title.y = element_text(vjust = 0)) +
  theme(legend.box.margin=ggplot2::margin(0,-1,-5,-1))+
  theme(legend.margin=ggplot2::margin(0.15,0,-0.15,0), legend.key.width = unit(0.3, 'cm'), legend.key.height = unit(0.0, 'cm'))
b4
rm(b,b1,b2,b3)

#PATCHWORK FIGURE 1 + 2####
library(patchwork)
library(grid)
b_final <- b4 +
  theme(
    text = element_text(family = "Arial", size = 6.5),
    axis.text.x = element_text(size = 11) ,
    axis.text.y = element_text(size = 8),
    # axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    # legend.title = element_text(size = 6),
    strip.text.x = element_text(size = 12)
  )
b_final <- b4 +
  theme(
    text = element_text(family = "Arial", size = 6.5),
    axis.text.x = element_text(size = 12) ,
    axis.text.y = element_text(size = 10),
    # axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    # legend.title = element_text(size = 6),
    strip.text.x = element_text(size = 12,margin=margin(t=3,r=3,b=3,l=3,"pt"))
  )

fig1_mod=fig1+theme(axis.text.x=element_text(angle=45,hjust=1),strip.text.x = element_text(size = 12,margin=margin(t=3,r=3,b=3,l=3,"pt")))
fig2_mod=fig2+theme(axis.text.x=element_text(angle=45,hjust=1,size = 11),axis.text.y = element_text(size = 8))
combined_fig <- fig1_mod + fig2_mod +
  plot_layout(ncol = 2, widths = c(3, 1)) #+
# plot_annotation(tag_levels = '1') &
# theme(plot.tag = element_text(size = 14, face = "bold")) # adds labels “a”, “b”
combined_fig
# combined_fig <- fig1 + fig2 +
#   plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") & 
#   theme(plot.margin = unit(c(0.1, 0.4, 0.1, 0.1), "cm")) +
#   plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")")
ggsave(
  filename = "C:/Users/Test/Documents/MS THESIS/PaperDrafts/FinalFigs/Figure2_3_combined.pdf",
  plot = combined_fig,
  width = 180, height = 120, units = "mm",
  device = cairo_pdf,
  dpi = 300
)


