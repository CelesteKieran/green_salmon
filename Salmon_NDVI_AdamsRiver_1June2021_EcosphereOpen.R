#-------------------------------------

# Title: Adams River script, Open access
# Authors: Celeste Kieran, Debora Obrist, Nico Muñoz
# Start date: 1 June 2021

#-------------------------------------
getwd()
rm(list=ls())


# We analyzed spatial and temporal associations between salmon abundance and riparian forest productivity (measured as NDVI) from a 36-year time series at the Adams River, British Columbia, Canada.
# This script shows our main models, model diagnostics, and predictions.
# Main datasets can be accessed upon request or at https://github.com/CelesteKieran/Babine_Scale_Isotopes_2021
# Additional information and data is freely available upon request, contact Celeste Kieran at celestekieran@gmail.com




#### ------------------- Data collection and preparation -------------- ####



# Salmon abundance data from 1983 - 2018 was retrieved from the DFOs salmon escapement 
# database (NuSEDS) at: https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6



# Weather data was retrieved using the R package weathercan:
# SILVER CREEK	(Station ID: 1327), 	50.55	-119.35, Active: 1989-2020, Elev: 419.0
# EAGLE BAY	(Station ID: 1258),	50.93	-119.22, Active: 1924-1995, Elev: 411.50	



# Plots were selected using Google Earth Pro, with 'close' plots directly along the river bank 
# and 'far' plots 65 metres behind the close plots. Each close far pair constitutes a 'block'  



# Distance downstream (metres) and elevation (metres) were measured using Google Earth Pro.



# NDVI imagery data from 1984 - 2019 was retrieved using Google Earth Engine from Landsat 5 and 7, and U.S. Geological Survey LEDAPS atmospherically corrected Tier 1 surface reflectance with the pixel QA band was used to remove interference. For growing season analysis, NDVI from July - September was used, for fall analysis NDVI from September - November was used. In order to retain only densely vegetated plots and to avoid disturbances to vegetation, individual NDVI observations below 0.3 were removed and all plots with a mean growing season NDVI of less than 0.5 were removed from analyses. We then took the mean NDVI for each plot for both summer (July - September) and fall (September - November). Our un-processed NDVI dataset is available upon request. 


#### ------------------- Main datasets -------------------------------- ####


# Summer dataframe
grow_df <- read.csv("growing_season_df_adams_ecosphere_2june2021.csv")




# Fall dataframe
fall_df <- read.csv("fall_season_df_adams_ecosphere_2june2021.csv")






#### ------------------- Call packages -------------------------------- ####
library(performance)
library(glmmTMB)
library(AICcmodavg)
library(DHARMa)
#### ------------------- Summer models -------------------------------- ####



## Summer models: testing effects of salmon on riparian greenness (NDVI) of the following year



glmm1 <-  glmmTMB(growing_season_mean_NDVI ~ combined_salmon_abundance_scaled + distance_from_bank + elevation_change_scaled + summer_mean_temp_scaled + summer_total_precip_scaled + metres_downstream_scaled  + (1|block/plot.id),
                  family=beta_family(),
                  data = grow_df)


glmm2 <-  glmmTMB(growing_season_mean_NDVI ~ combined_salmon_abundance_scaled * distance_from_bank + combined_salmon_abundance_scaled * elevation_change_scaled + summer_mean_temp_scaled + summer_total_precip_scaled + metres_downstream_scaled + 
                    (1|block/plot.id),
                  family=beta_family(),
                  data = grow_df)


glmm3 <-  glmmTMB(growing_season_mean_NDVI ~ pulse_year + distance_from_bank + elevation_change_scaled + summer_mean_temp_scaled + summer_total_precip_scaled + metres_downstream_scaled  + (1|block/plot.id),
                  family=beta_family(),
                  data = grow_df)


glmm4 <-  glmmTMB(growing_season_mean_NDVI ~ pulse_year * distance_from_bank + pulse_year * elevation_change_scaled + summer_mean_temp_scaled + summer_total_precip_scaled + metres_downstream_scaled + 
                    (1|block/plot.id),
                  family=beta_family(),
                  data = grow_df)


glmm5 <-  glmmTMB(growing_season_mean_NDVI ~ distance_from_bank + elevation_change_scaled + summer_mean_temp_scaled + summer_total_precip_scaled + metres_downstream_scaled + (1|block/plot.id),
                  family=beta_family(),
                  data = grow_df)






## Save RDS
saveRDS(glmm1, "mod1_add-2june2021.rds")
saveRDS(glmm2, "mod2_int-2june2021.rds")
saveRDS(glmm3, "mod3_pulse-2june2021.rds")
saveRDS(glmm4, "mod4_pulse-2june2021.rds")
saveRDS(glmm5, "mod5_control-2june2021.rds")



#### Diagnostics ####

# Variance Inflation Factors (VIFs) for all models was checked to avoid multicollinearity 
# among independent variables.
# Akaike Information Criterion (AIC) was used to estimate model quality 
# The DHARMa package (by Florian Hartig) was used for residual diagnostics of our top model
# (as selected by AIC)


## VIF ##


check_collinearity(glmm1) 
check_collinearity(glmm2) 
check_collinearity(glmm3) 
check_collinearity(glmm4)
check_collinearity(glmm5)




## AIC ##


glmm1 <- readRDS( "mod1_add-2june2021.rds")
glmm2 <- readRDS("mod2_int-2june2021.rds")
glmm3 <- readRDS("mod3_pulse-2june2021.rds")
glmm4 <- readRDS("mod4_pulse-2june2021.rds")
glmm5 <- readRDS( "mod5_control-2june2021.rds")

## setup a subset of models of Table 1
Cand.models <- list( )
Cand.models[[1]] <- glmm1
Cand.models[[2]] <- glmm2
Cand.models[[3]] <- glmm3
Cand.models[[4]] <- glmm4
Cand.models[[5]] <- glmm5

# create a vector of names to trace back models in set
Modnames <- c("1. salmon abundance + distance_from_bank + elev_change", 
              "2. salmon abundance * distance_from_bank + salmon abundance * elev_change",
              "3. sockeye pulse + distance_from_bank + elev_change",
              "4. sockeye pulse * distance_from_bank + pulse_year * elev_change", 
              "5. distance_from_bank + elev_change")

## generate AIC table
aictab(cand.set = Cand.models, modnames = Modnames, second.ord = FALSE)
## round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames),
      digits = 4, LL = TRUE, second.ord = FALSE)


# AIC suggests model 3, which tests additive effects of sockeye pulse, elevation change, 
# and distance from the riverbank
summary(glmm3)






## DHARMa ##
# "a simulation-based approach to create readily interpretable scaled (quantile) residuals for fitted (generalized) linear mixed models" - https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html



# Diagnostics for model 3, our top AIC model.
mod3 <- glmm3


# Calculate residuals
simulationOutput <- simulateResiduals(fittedModel = mod3, plot = T)



# Plot scaled residuals
plot(simulationOutput)
# Signifiant deviation from assumed distribution (ie. KS test) was decided to be acceptable given the large number of data points (5472 NDVI observations)


# Plot residuals against predictors 
# (add quantreg = T for higher diagnostics)
plotResiduals(simulationOutput, grow_df$pulse_year, quantreg = T)
plotResiduals(simulationOutput, grow_df$combined_salmon_abundance_scaled, quantreg = T)
plotResiduals(simulationOutput, grow_df$summer_mean_temp_scaled, quantreg = T)
plotResiduals(simulationOutput, grow_df$summer_total_precip_scaled, quantreg = T)
plotResiduals(simulationOutput, grow_df$distance_from_bank, quantreg = T)
plotResiduals(simulationOutput, grow_df$metres_downstream, quantreg = T)
plotResiduals(simulationOutput, grow_df$elev_change, quantreg = T)
plotResiduals(simulationOutput, grow_df$sat, quantreg = T)
plotResiduals(simulationOutput, grow_df$plot.id, quantreg = T)
# All continuous variable show quantile deviations

# Histogram of residuals
hist(simulationOutput)







#### ------------------- Fall models ---------------------------------- ####



## Fall models: testing effects of salmon on riparian greenness (NDVI) of the same year
fall_glmm1 <-  glmmTMB(fall_mean_NDVI ~ combined_salmon_abundance_scaled + distance_from_bank + elevation_change_scaled + mean_fall_temp_scaled + total_fall_precip_scaled + metres_downstream_scaled  + (1|block/plot.id),
                  family=beta_family(),
                  data = fall_df)


fall_glmm2 <-  glmmTMB(fall_mean_NDVI ~ combined_salmon_abundance_scaled * distance_from_bank + combined_salmon_abundance_scaled * elevation_change_scaled + mean_fall_temp_scaled + total_fall_precip_scaled + metres_downstream_scaled + 
                    (1|block/plot.id),
                  family=beta_family(),
                  data = fall_df)


fall_glmm3 <-  glmmTMB(fall_mean_NDVI ~ pulse_year + distance_from_bank + elevation_change_scaled + mean_fall_temp_scaled + total_fall_precip_scaled + metres_downstream_scaled  + (1|block/plot.id),
                  family=beta_family(),
                  data = fall_df)


fall_glmm4 <-  glmmTMB(fall_mean_NDVI ~ pulse_year * distance_from_bank + pulse_year * elevation_change_scaled + mean_fall_temp_scaled + total_fall_precip_scaled + metres_downstream_scaled + 
                    (1|block/plot.id),
                  family=beta_family(),
                  data = fall_df)


fall_glmm5 <-  glmmTMB(fall_mean_NDVI ~ distance_from_bank + elevation_change_scaled + mean_fall_temp_scaled + total_fall_precip_scaled + metres_downstream_scaled + (1|block/plot.id),
                  family=beta_family(),
                  data = fall_df)




## Save RDS
saveRDS(fall_glmm1, "FALLmod1_add-2june2021.rds")
saveRDS(fall_glmm2, "FALLmod2_int-2june2021.rds")
saveRDS(fall_glmm3, "FALLmod3_pulse-2june2021.rds")
saveRDS(fall_glmm4, "FALLmod4_pulse-2june2021.rds")
saveRDS(fall_glmm5, "FALLmod5_control-2june2021.rds")

#### Diagnostics ####

## VIF ##

check_collinearity(fall_glmm1) 
check_collinearity(fall_glmm2) 
check_collinearity(fall_glmm3) 
check_collinearity(fall_glmm4)
check_collinearity(fall_glmm5)




## AIC ##

fall_glmm1 <- readRDS( "FALLmod1_add-2june2021.rds")
fall_glmm2 <- readRDS("FALLmod2_int-2june2021.rds")
fall_glmm3 <- readRDS("FALLmod3_pulse-2june2021.rds")
fall_glmm4 <- readRDS("FALLmod4_pulse-2june2021.rds")
fall_glmm5 <- readRDS( "FALLmod5_control-2june2021.rds")

## setup a subset of models of Table 1
Cand.models <- list( )
Cand.models[[1]] <- fall_glmm1
Cand.models[[2]] <- fall_glmm2
Cand.models[[3]] <- fall_glmm3
Cand.models[[4]] <- fall_glmm4
Cand.models[[5]] <- fall_glmm5

# create a vector of names to trace back models in set
Modnames <- c("1. salmon abundance + distance_from_bank + elev_change", 
              "2. salmon abundance * distance_from_bank + salmon abundance * elev_change",
              "3. sockeye pulse + distance_from_bank + elev_change",
              "4. sockeye pulse * distance_from_bank + pulse_year * elev_change", 
              "5. distance_from_bank + elev_change")

## generate AIC table
aictab(cand.set = Cand.models, modnames = Modnames, second.ord = FALSE)
## round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames),
      digits = 4, LL = TRUE, second.ord = FALSE)


# AIC suggests model 4, which tests interactive effects of sockeye pulse, elevation change, and distance from the riverbank
summary(fall_glmm4)






## DHARMa ##


# Diagnostics for model 4, our top fall AIC model.
mod4 <- fall_glmm4


# Calculate residuals
simulationOutput <- simulateResiduals(fittedModel = mod4, plot = T)



# Plot scaled residuals
plot(simulationOutput)
# Signifiant deviation from assumed distribution (ie. KS test) was decided to be acceptable given the large number of data points (5222 NDVI observations)


# Plot residuals against predictor
# (add quantreg = T for higher diagnostics)
plotResiduals(simulationOutput, fall_df$pulse_year, quantreg = T)
plotResiduals(simulationOutput, fall_df$combined_salmon_abundance_scaled, quantreg = T)
plotResiduals(simulationOutput, fall_df$mean_fall_temp_scaled, quantreg = T)
plotResiduals(simulationOutput, fall_df$total_fall_precip_scaled, quantreg = T)
plotResiduals(simulationOutput, fall_df$distance_from_bank, quantreg = T)
plotResiduals(simulationOutput, fall_df$metres_downstream, quantreg = T)
plotResiduals(simulationOutput, fall_df$elev_change, quantreg = T)
plotResiduals(simulationOutput, fall_df$sat, quantreg = T)
plotResiduals(simulationOutput, fall_df$plot.id, quantreg = T)
# All continuous variable show quantile deviations

# Histogram of residuals
hist(simulationOutput)






#### ------------------- 4-year cohort fall model --------------------- ####


# This model tests effect of sockeye cohort year (1-dominant cohort, 2-subdominant, 3/4-off) on riparian NDVI of the same fall season
fall_df$pulse_plus <- as.factor(fall_df$pulse_plus)


fall4yr_glmm1 <-  glmmTMB(fall_mean_NDVI ~ pulse_plus + distance_from_bank + elevation_change_scaled + mean_fall_temp_scaled + total_fall_precip_scaled + metres_downstream_scaled  + (1|block/plot.id),
                       family=beta_family(),
                       data = fall_df)

## Save RDS
saveRDS(fall4yr_glmm1, "mod1_4yrfall_add-8june2021.rds")

#### Diagnostics ####

fall4yr_glmm1 <- readRDS("mod1_4yrfall_add-8june2021.rds")




## VIF ##
check_collinearity(fall4yr_glmm1) 




## DHARMa ##

# Calculate residuals
simulationOutput <- simulateResiduals(fittedModel = fall4yr_glmm1, plot = T)



# Plot scaled residuals
plot(simulationOutput)



# Plot residuals against predictor
# (add quantreg = T for higher diagnostics)
plotResiduals(simulationOutput, fall_df$pulse_plus, quantreg = T)
plotResiduals(simulationOutput, fall_df$combined_salmon_abundance_scaled, quantreg = T)
plotResiduals(simulationOutput, fall_df$mean_fall_temp_scaled, quantreg = T)
plotResiduals(simulationOutput, fall_df$total_fall_precip_scaled, quantreg = T)
plotResiduals(simulationOutput, fall_df$distance_from_bank, quantreg = T)
plotResiduals(simulationOutput, fall_df$metres_downstream, quantreg = T)
plotResiduals(simulationOutput, fall_df$elev_change, quantreg = T)
plotResiduals(simulationOutput, fall_df$sat, quantreg = T)
plotResiduals(simulationOutput, fall_df$plot.id, quantreg = T)
# All continuous variable show quantile deviations

# Histogram of residuals
hist(simulationOutput)


#### ------------------- Results / Predictions ------------------------ --------------------------------------#
# Both top seasonal models (as suggested by AIC) used the binary varible 'pulse year' as a salmon influence  measure. Here we predict riparian NDVI at the Adams River during pulse years (dominant cohort years) and non-pulse years (sub-dominant or 'off' cohort years) for both summer and fall models, as well as for our model testing each effects of each individual cohort year




#### Summer model ####

# Predictions from our top summer model (glmm3) as chosen by AIC, 
# which tests effects of a dominant cohort year (pulse year) on 
# riparian NDVI of the following summer at the Adams River


# Pulse year vs non-pulse year at plots alongside the river (CLOSE plots): 
newdat_pulse <- expand.grid(summer_mean_temp_scaled = mean(grow_df$summer_mean_temp_scaled), 
                            summer_total_precip_scaled = mean(grow_df$summer_total_precip_scaled), 
                            metres_downstream_scaled = mean(grow_df$metres_downstream_scaled),
                            elevation_change_scaled = mean(grow_df$elevation_change_scaled),
                            pulse_year = unique(grow_df$pulse_year),
                            distance_from_bank = "CLOSE",
                            block = NA,
                            plot.id = NA)

pred_pulse <- as.data.frame(cbind(newdat_pulse, 
                                  grow_mean = predict(glmm3, newdat_pulse, re.form = NULL, type = "response", se.fit = TRUE)))

pred_pulse$grow_mean <- pred_pulse$grow_mean.fit
pred_pulse$lwr <- pred_pulse$grow_mean - pred_pulse$grow_mean.se.fit * 1.96
pred_pulse$upr <- pred_pulse$grow_mean + pred_pulse$grow_mean.se.fit * 1.96

pred_pulse$grow_mean.se.fit * 1.96 # Confidence intervals (CIs)
# for both pulse and non-pulse years round to 0.008

pred_pulse 
# In a non-pulse year, predicted riparian NDVI at CLOSE plots is 0.687 +/- 0.008.
# In a pulse year,  predicted NDVI is 0.703 +/- 0.008






# Effect on NDVI of pulse years vs non pulse years for plots ~65 metres from the river bank (FAR plots)
newdat_pulse <- expand.grid(summer_mean_temp_scaled = mean(grow_df$summer_mean_temp_scaled), 
                            summer_total_precip_scaled = mean(grow_df$summer_total_precip_scaled), 
                            metres_downstream_scaled = mean(grow_df$metres_downstream_scaled),
                            elevation_change_scaled = mean(grow_df$elevation_change_scaled),
                            pulse_year = unique(grow_df$pulse_year),
                            distance_from_bank = "FAR",
                            block = NA,
                            plot.id = NA)

pred_pulse <- as.data.frame(cbind(newdat_pulse, 
                                  grow_mean = predict(glmm3, newdat_pulse, re.form = NULL, type = "response", se.fit = TRUE)))

pred_pulse$grow_mean <- pred_pulse$grow_mean.fit
pred_pulse$lwr <- pred_pulse$grow_mean - pred_pulse$grow_mean.se.fit * 1.96
pred_pulse$upr <- pred_pulse$grow_mean + pred_pulse$grow_mean.se.fit * 1.96

pred_pulse$grow_mean.se.fit * 1.96 # CIs round to 0.007

pred_pulse
# In a non-pulse year, predicted NDVI is 0.762 +/- 0.007
# In a pulse year, predicted NDVI is 0.776 +/- 0.007



#### Fall model ####

# Predictions from our top fall model (glmm3) as chosen by AIC, 
# which tested the effects of a dominant cohort year (pulse year) 
# on riparian NDVI of the same fall season as the dominant spawning event at the Adams River


# Pulse year vs non-pulse year at plots alongside the river (CLOSE plots): 
newdat_fall <- expand.grid(mean_fall_temp_scaled = mean(fall_df$mean_fall_temp_scaled), 
                           total_fall_precip_scaled = mean(fall_df$total_fall_precip_scaled), 
                           metres_downstream_scaled = mean(fall_df$metres_downstream_scaled),
                           elevation_change_scaled = mean(fall_df$elevation_change_scaled),
                           pulse_year = unique(fall_df$pulse_year),
                           distance_from_bank = "CLOSE",
                           sat = "LANDSAT7",
                           block = NA,
                           plot.id = NA)


pred_pulse <- as.data.frame(cbind(newdat_fall, 
                                  mean_NDVI = predict(fall_glmm4, newdat_fall, re.form = NULL, 
                                                      type = "response", se.fit = TRUE)))


pred_pulse$mean_NDVI <- pred_pulse$mean_NDVI.fit
pred_pulse$lwr <- pred_pulse$mean_NDVI - pred_pulse$mean_NDVI.se.fit * 1.96
pred_pulse$upr <- pred_pulse$mean_NDVI + pred_pulse$mean_NDVI.se.fit * 1.96

pred_pulse$mean_NDVI.se.fit * 1.96 

pred_pulse
# In a pulse year, NDVI is 0.619 +/- 0.009
# In a non-pulse year, NDVI is 0.632 +/- 0.007






# Pulse years vs non pulse years for plots ~65 metres from the river bank (FAR plots)
newdat_fall_far <- expand.grid(mean_fall_temp_scaled = mean(fall_df$mean_fall_temp_scaled), 
                               total_fall_precip_scaled = mean(fall_df$total_fall_precip_scaled), 
                               metres_downstream_scaled = mean(fall_df$metres_downstream_scaled),
                               elevation_change_scaled = mean(fall_df$elevation_change_scaled),
                               pulse_year = unique(fall_df$pulse_year),
                               distance_from_bank = "FAR",
                               sat = "LANDSAT7",
                               block = NA,
                               plot.id = NA)


pred_pulse <- as.data.frame(cbind(newdat_fall, 
                                  mean_NDVI = predict(fall_glmm4, newdat_fall_far, re.form = NULL, type = "response", se.fit = TRUE)))


pred_pulse$mean_NDVI <- pred_pulse$mean_NDVI.fit
pred_pulse$lwr <- pred_pulse$mean_NDVI - pred_pulse$mean_NDVI.se.fit * 1.96
pred_pulse$upr <- pred_pulse$mean_NDVI + pred_pulse$mean_NDVI.se.fit * 1.96

pred_pulse$mean_NDVI.se.fit * 1.96 

pred_pulse

# For far plots, in a pulse year, NDVI is 0.682 +/- 0.008
# In a non-pulse year, NDVI is  0.705 +/- 0.006




#### Cohort model ####

# Predictions from our fall model that tests effects of cohort year
# (1-dominant cohort, 2-subdominant, 3/4-off years) on riparian NDVI
# of same fall season at the Adams River


# Pulse year vs non-pulse year at plots alongside the river (CLOSE plots): 
fall_df$pulse_plus <- as.factor(fall_df$pulse_plus)

newdat_fall <- expand.grid(mean_fall_temp_scaled = mean(fall_df$mean_fall_temp_scaled), 
                           total_fall_precip_scaled = mean(fall_df$total_fall_precip_scaled), 
                           metres_downstream_scaled = mean(fall_df$metres_downstream_scaled),
                           elevation_change_scaled = mean(fall_df$elevation_change_scaled),
                           pulse_plus = unique(fall_df$pulse_plus),
                           distance_from_bank = "CLOSE",
                           sat = "LANDSAT7",
                           block = NA,
                           plot.id = NA)

pred_pulse <- as.data.frame(cbind(newdat_fall, 
                                  mean_NDVI = predict(fall4yr_glmm1, newdat_fall, re.form = NULL, type = "response", se.fit = TRUE)))


pred_pulse$mean_NDVI <- pred_pulse$mean_NDVI.fit
pred_pulse$lwr <- pred_pulse$mean_NDVI - pred_pulse$mean_NDVI.se.fit * 1.96
pred_pulse$upr <- pred_pulse$mean_NDVI + pred_pulse$mean_NDVI.se.fit * 1.96

pred_pulse$mean_NDVI.se.fit * 1.96 

pred_pulse
# For close plots
# In a pulse year (dominant / 1), NDVI is 0.615 +/- 0.008
# In a year following a pulse year (subdominant / 2), NDVI is 0.650 +/- 0.008
# Cohort year (off / 3), NDVI is 0.620 +/- 0.008
# Cohort year (off / 4), NDVI is 0.632 +/- 0.008





# Pulse year vs non-pulse year at plots ~65 metres from the riverbank (FAR plots): 
newdat_fall <- expand.grid(mean_fall_temp_scaled = mean(fall_df$mean_fall_temp_scaled), 
                           total_fall_precip_scaled = mean(fall_df$total_fall_precip_scaled), 
                           metres_downstream_scaled = mean(fall_df$metres_downstream_scaled),
                           elevation_change_scaled = mean(fall_df$elevation_change_scaled),
                           pulse_plus = unique(fall_df$pulse_plus),
                           distance_from_bank = "FAR",
                           sat = "LANDSAT7",
                           block = NA,
                           plot.id = NA)


pred_pulse <- as.data.frame(cbind(newdat_fall, 
                                  mean_NDVI = predict(fall4yr_glmm1, newdat_fall, re.form = NULL, type = "response", se.fit = TRUE)))


pred_pulse$mean_NDVI <- pred_pulse$mean_NDVI.fit
pred_pulse$lwr <- pred_pulse$mean_NDVI - pred_pulse$mean_NDVI.se.fit * 1.96
pred_pulse$upr <- pred_pulse$mean_NDVI + pred_pulse$mean_NDVI.se.fit * 1.96

pred_pulse$mean_NDVI.se.fit * 1.96 # rounds to 0.007 for all

pred_pulse
#  Far plots, with Landsat 7, and all other variables held at their means:
# In a pulse year (1), NDVI is 0.687 +/- 0.007
# In a year following pulse year (2), NDVI is 0.719 +/- 0.007
# Year 3, NDVI is 0.691 +/- 0.007
# Year 4, NDVI is 0.718 +/- 0.007





