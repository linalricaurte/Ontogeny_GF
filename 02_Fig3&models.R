rm(list=ls())
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Make list of packages needed to run our code
packages <- c("lubridate","tidyverse","fossil","ggplot2","sf","move","suntools","tidyterra","maptiles","rnaturalearth","glmmTMB","lme4",
              "lmerTest","performance","ggeffects", "multcomp", "gridExtra", "emmeans", "DHARMa" )
# apply ipak() function to list of packages
ipak(packages)


###function overdispertion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson_chisq <- sum(rp^2)
  dispersion <- Pearson_chisq / rdf
  return(dispersion)
}

l.stage<-read.csv('./Data/AB/long.events.dist25km.3d.v3.csv', header = TRUE, dec = ".", sep=",") #I create this file from the 01_fig_1_tracking_data_JAE
sum(duplicated(l.stage$dist.segment))

str(l.stage)
l.stage$season<-as.factor(l.stage$season)
l.stage$type2<-as.factor(l.stage$type2)
l.stage$year<-as.factor(l.stage$year)
l.stage$colony<-as.factor(l.stage$colony)
l.stage$dev<-as.factor(l.stage$dev)


### summarise by type2 and season
l.stage_summary <- l.stage %>%
  group_by(dev,type2,season) %>%
  dplyr::summarise(
    n.stops = n_distinct(dist.segment),
    sex = first(sex),
    colony = first(colony),
    year = first(year),
    .groups = "drop"
  )

l.stage.for.mod<-l.stage_summary

###############################################################################
# Note: The code below to investigate the effect of age on staying behaviour of greater flamingos
# was adapted from:
# Marta Acácio (2025). Ageing Vultures. GitHub repository:
# https://github.com/msa2015/Ageing_Vultures
# Modifications were made to meet specific project needs.
###############################################################################


########################################
### Number of resident events###
########################################
l.stage.for.mod$type2 <- as.numeric(gsub("yo", "", l.stage.for.mod$type2))
l.stage.for.mod$type2_scaled <- scale(l.stage.for.mod$type2)

summary(mfull)

### #re-scale variables

events.age.scale <- attr(l.stage.for.mod$type2_scaled,"scaled:scale")
events.age.center <- attr(l.stage.for.mod$type2_scaled,"scaled:center")
events.age.scale <- attr(l.stage.for.mod$type2_scaled,"scaled:scale")
events.age.center <- attr(l.stage.for.mod$type2_scaled,"scaled:center")

# For the models

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # 

func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

events_models <- list()

events_mod_summary <- list()

# Loop through all 4 functional relationships, 
for(i in 1:length(functional_relationships)){
  
  # Construct the model formula dynamically
  model_formula <- paste("n.stops  ~ season * ", functional_relationships[[i]], 
                         "+ (1|dev)+(1|year)+(1|colony)", sep = "")
  
  # Fit the model
  mod <- glmmTMB(formula = as.formula(model_formula),
                 REML = F,
                 family = compois(link = "log"),
                 data = l.stage.for.mod)
  
  # Store the models and their summaries
  events_models[[i]] <- mod
  events_mod_summary[[i]] <- summary(mod)
  
  # Cleanup for the next iteration
  rm(mod)
  rm(model_formula)
}

#

# Compare the performance of all models,
events_comparison <- compare_performance(events_models[[1]], events_models[[2]], 
                                         events_models[[3]],events_models[[4]])  # 

events_comparison$Model <- func_rel_names

# Order models by AIC (from lowest to highest)
events_comparison <- events_comparison %>% arrange(AIC)

events_comparison$Model <- sapply(events_comparison$Model, as.character)

write.csv(events_comparison, "./tables/forRevAB/AIC_LDstops_compois.csv")

#### best fit model or more parsimonious  
mfull <- glmmTMB(n.stops ~ type2_scaled + (1|dev)+(1|year)+(1|colony),
                 family = compois(link = "log"),
                 data = l.stage.for.mod)

summary(mfull)
# Display comparison results
print(events_comparison)

###-------------------------
### build figure 3a 
###-------------------------

# Predict values for the selected model (e.g., third-degree polynomial or age as a factor)
events_mod_df <- ggpredict(events_models[[1]], terms = c("type2_scaled[all]", "season"))
events_mod_df$type2 <- round(events_mod_df$x * events.age.scale + events.age.center)
events_mod_df <- dplyr::rename(events_mod_df,  season= group)

pd <- position_dodge(0.3)

p <- ggplot(data = events_mod_df) +
  
  geom_line(aes(x = type2, y = predicted, color = season),
            linewidth = 0.7) +
  stat_summary(data = events_mod_df,aes(x = type2, y = predicted, color = season),fun.data = "mean_sdl",fun.args = list(mult = 1),
    geom = "pointrange",size = 0.5,position = position_dodge(0.3)) +
  geom_errorbar(data = events_mod_df,
                aes(ymin = conf.low, ymax = conf.high, x = type2, colour = season),
                width = 0, position = position_dodge(width = 0.3),
                linewidth = 0.6) +
  
  scale_color_manual(values = c("Autumn-Winter"="#c7b42e","Spring-Summer"="cornflowerblue"
  )) +
  labs(title = "", x = "", y = "Number of LD\n stays") +
  geom_point(data = l.stage.for.mod,aes(y = n.stops, x = type2, colour = season),size = 2, alpha = 0.1,
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.2)) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    axis.line         = element_blank(),
    axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    text              = element_text(face = "plain"),
    axis.text         = element_text(size = 13, face = "plain"),
    axis.text.x       = element_text(size = 13, face = "plain"),
    axis.text.y       = element_text(size = 13, face = "plain"),
    axis.title        = element_text(size = 13, face = "plain"),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 13, face = "plain"),
    legend.position   = "none",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 14, face = "plain"),
    strip.text        = element_text(size = 12, face = "plain"),
    strip.background  = element_blank()
  )+
  guides(
    alpha = "none",
    size  = guide_legend(nrow = 2)
  )
p

mfull <- glmmTMB(n.stops ~ type2_scaled * season + (1|dev)+(1|year)+(1|colony),
                 family = compois(link = "log"),
                 #family = poisson,
                 data = l.stage.for.mod)

summary(mfull)
## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)

testOutliers(tmp, type = "bootstrap")
testDispersion(tmp)

# Run the outlier test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outlier test results
print(outlier_test) ### not extreme aoutliers
testDispersion(mfull) #, ok
overdisp_fun(mfull) # ok

par(mar = c(5, 4, 4, 2) + 0.1)  # default margins
check_collinearity(mfull)##ok
diagnostics.plot(mfull)# 


### select the best model using dredge in this case the interaction doest improvve model fit

d1 <- dredge(mfull, rank = "AIC") 
sw(d1)
best <- get.models(d1,subset=NA) 
Anova(mfull) 
Anova(best[[1]]) 

selection.table <- model.sel(d1) 
selection.table 
# Convert to data.frame
selection.df <- as.data.frame(selection.table)

#write.table(selection.df,"tables/forRevAB/dredge_LDresident.xlsx",quote=FALSE,row.names = TRUE)
library(openxlsx)
write.xlsx(selection.df,file = "tables/forRevAB/dredge_LDresident.xlsx",rowNames = TRUE,overwrite = TRUE)

# ### get R con and R marginal
# Pre-allocate matrix: rows = number of models, cols = R2m, R2c
r2 <- matrix(NA, ncol = 2, nrow = length(best))
colnames(r2) <- c("R2m", "R2c")

for (i in seq_along(best)) {
  result <- r.squaredGLMM(best[[i]])
  r2[i, ] <- result[1, ]  # 1st row is marginal, 2nd row conditional; or just take first row
}
r2


#############################################
###  Fig 3b stay duration in days
#############################################
l.stage<-read.csv('./Data/AB/long.events.dist25km.3d.v2.csv', header = TRUE, dec = ".", sep=",") #I create this file from the 01_fig_1_tracking_data_JAE

str(l.stage)
l.stage$season<-as.factor(l.stage$season)
l.stage$type2<-as.factor(l.stage$type2)

### compute the mean duration per device, age and season
l.stage <- l.stage %>%
  group_by(dev,sex, type2, season,colony) %>%
  dplyr::summarise(
    mean.bout.dur = mean(bout.dur, na.rm = TRUE),
    median.bout.dur = median(bout.dur, na.rm = TRUE),
    year = first(year), 
    .groups = "drop"
  )

l.stage$type2 <- as.numeric(gsub("yo", "", l.stage$type2))  # Remove 'yo' and convert to numeric
l.stage$type2_scaled <- scale(l.stage$type2)

##### check the model residuals to see normality
# ### model selection selects season and age (not the interation) 
mfull <- lmer(mean.bout.dur ~ type2_scaled + season + (1 | dev) + (1 | year) + (1 | colony),
              data = l.stage
)

## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)

### #re-scale variables

duration.age.scale <- attr(l.stage$type2_scaled,"scaled:scale")
duration.age.center <- attr(l.stage$type2_scaled,"scaled:center")


func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  

duration_models <- list()
duration_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  model_formula <- paste("mean.bout.dur ~ season * ", functional_relationships[[i]], 
                         " + (1|dev)+(1|year)+(1|colony)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,
                 family=Gamma(link = "log"),
                 data = l.stage)
  
  duration_models[[i]] <- mod
  duration_mod_summary[[i]] <- summary(mod)
  
  rm(mod)
  rm(model_formula)
  
}

# Compare the performance of all models, including the categorical age model
duration_comparison <- compare_performance(duration_models[[1]], duration_models[[2]], 
                                           duration_models[[3]],duration_models[[4]])  # Add the new model
duration_comparison$Model <- func_rel_names

# Display comparison results
print(duration_comparison)


# Order models by AIC (from lowest to highest)
duration_comparison <- duration_comparison %>% arrange(AIC)
duration_comparison$Model <- sapply(duration_comparison$Model, as.character)

write.xlsx(duration_comparison, "./tables/forRevAB/AIC_duration_AB.xls")

# List of model names
model_names <- c("Model 1: Linear", "Model 2: Exponential", "Model 3: Quadratic", "Model 4: poly")

# Compute R2 for each model
r2_values <- lapply(duration_models, r.squaredGLMM)

# Combine the results into a data frame
r2_df <- do.call(rbind, r2_values)
colnames(r2_df) <- c("Marginal R²", "Conditional R²")

# Add model names as a column
r2_df <- cbind(Model = model_names, r2_df)

# Print the result
print(r2_df)

write.csv(r2_df, "./tables/v2_JAE/r2_distance_log_JAE.csv")


###--------------------------
###figure 3b################
###--------------------------
duration_mod_df <- ggpredict(duration_models[[1]], terms = c("type2_scaled[all]", "season"))
duration_mod_df$type2 <- round(duration_mod_df$x * duration.age.scale + duration.age.center)
duration_mod_df <- dplyr::rename(duration_mod_df, season = group)

p2 <- ggplot(data = duration_mod_df) +
  geom_line(aes(x = type2, y = predicted, color = season),
            linewidth = 0.7) +
  stat_summary(data = duration_mod_df,aes(x = type2, y = predicted, color = season),fun.data = "mean_sdl",
    fun.args = list(mult = 1),geom = "pointrange",size = 0.5,position = position_dodge(0.3)) +
  
  geom_errorbar(data = duration_mod_df,
                aes(ymin = conf.low, ymax = conf.high, x = type2, colour = season),
                width = 0, position = position_dodge(width = 0.3),
                linewidth = 0.6) +
  geom_point(data = l.stage,
             aes(y = mean.bout.dur, x = type2, colour = season),
             size = 2, alpha = 0.2,
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.2)) +

  scale_color_manual(values = c("Autumn-Winter"="#c7b42e","Spring-Summer"="cornflowerblue"
  )) +
scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e","Spring-Summer" = "cornflowerblue")) +
  labs(title = "", x = "", y = "Mean LD\nstay duration") +
  theme_classic(base_size = 14) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    axis.line         = element_blank(),
    axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    text              = element_text(face = "plain"),
    axis.text         = element_text(size = 13, face = "plain"),
    axis.text.x       = element_text(size = 13, face = "plain"),
    axis.text.y       = element_text(size = 13, face = "plain"),
    axis.title        = element_text(size = 13, face = "plain"),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 13, face = "plain"),
    legend.position   = "none",
    legend.box.margin = margin(t = 0.02, r = 0, b = 0, l = 0, unit = "cm"),
    legend.text       = element_text(size = 14, face = "plain"),
    strip.text        = element_blank(),
    strip.background  = element_blank()
  ) +
  guides(
    alpha = "none",
    size  = guide_legend(nrow = 2)
  )
p2
####


# ### model selection 
mfull <- glmmTMB(
  mean.bout.dur ~ type2_scaled * season + (1 | dev) + (1 | year) + (1 | colony),
  family=Gamma(link = "log"),
  data = l.stage
)
options(scipen = 999)

summary(mfull)
### Compute some summary stadistics
tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)
summary(mfull)
anova(mfull)

## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)
# Run the outliers test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outliers test results
print(outlier_test) ### not extreme outliers

testDispersion(mfull) # ok
#overdisp_fun(mfull) #model with mean bout duration is overdispersed, 1219 with log(mean.bout.dur) I deal with overdispersion (0.35, ok)
check_overdispersion(mfull)
check_collinearity(mfull) ###ok

# Extract residuals and fitted values
residuals <- residuals(mfull, type = "pearson")  # Use Pearson residuals
fitted_values <- fitted(mfull)


### select the best model using dredge
options(na.action = "na.fail")  # required for dredge
d1 <- dredge(mfull, rank = "AIC",  # random effects locked
             trace = TRUE)
sw(d1)
summary(d1)
get.models(d1, subset = TRUE)
best <- get.models(d1,subset=NA) 
Anova(mfull) # 
Anova(best[[1]]) #
selection.table <- model.sel(d1) #
selection.table #
# Convert to data.frame
selection.df <- as.data.frame(selection.table)


#write.table(selection.df,"tables/forRevAB/dredge_LDresident.xlsx",quote=FALSE,row.names = TRUE)
library(openxlsx)
write.xlsx(selection.df,file = "tables/forRevAB/dredge_LDresidentduration.xlsx",rowNames = TRUE,overwrite = TRUE)

# ### get R con and R marginal
# Pre-allocate matrix: rows = number of models, cols = R2m, R2c
r2 <- matrix(NA, ncol = 2, nrow = length(best))
colnames(r2) <- c("R2m", "R2c")

for (i in seq_along(best)) {
  result <- r.squaredGLMM(best[[i]])
  r2[i, ] <- result[1, ]  # 1st row is marginal, 2nd row conditional; or just take first row
}
r2


###########################################
### distance to principal sites 
###############################
data2<-read.csv("./Data/AB/longest.stage_AB.csv", header = TRUE, dec = ".", sep=",") #
#remove nas
data2_clean <- data2 %>%
  filter(!is.na(distance_to_previous_km) & !is.na(prev_year_lat) & !is.na(prev_year_lon) & !is.na(prev_year))
data2_clean$type2 <- as.numeric(gsub("yo", "", data2_clean$type2))
data2_clean$type2_scaled <- scale(data2_clean$type2_num)
data2_clean$distance_scaled <- scale(data2_clean$distance_to_previous_km)

### #re-scale variables

# Scaling for distance
distance.age.scale <- attr(data2_clean$type2_scaled,"scaled:scale")
distance.age.center <- attr(data2_clean$type2_scaled,"scaled:center")

func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model


distance_models <- list()
distance_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  # Create the model formula
  model_formula <- paste("distance_to_previous_km ~ season * ", functional_relationships[[i]], 
                         " + (1|dev)+ (1|colony)+(1|year)", sep ="")
  
  # Fit the model using glmmTMB
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F, 
                 family=  Gamma(log),
                 # for model comparison
                 data = data2_clean_no_out)
  
  # Store the model and summary
  distance_models[[i]] <- mod
  distance_mod_summary[[i]] <- summary(mod)
  
  # Clean up the model and formula variables to free memory
  rm(mod)
  rm(model_formula)
  
}

# Compare the performance of all models, including the categorical age model
distance_comparison <- compare_performance(distance_models[[1]], distance_models[[2]], 
                                           distance_models[[3]],distance_models[[4]])  # Add the new model
distance_comparison$Model <- func_rel_names

# Display comparison results
print(distance_comparison)

# Order models by AIC (from lowest to highest)
distance_comparison <- distance_comparison %>% arrange(AIC)
distance_comparison$Model <- sapply(distance_comparison$Model, as.character)

write.xlsx(distance_comparison, "./tables/forRevAB/AIC_distancemainsites_AB.csv")

# List of model names
model_names <- c("Model 1: Exponential", "Model 2: Linear", "Model 3: Quadratic", "Model 4: poly")

# Compute R2 for each model
r2_values <- lapply(distance_models, r.squaredGLMM)

# Combine the results into a data frame
r2_df <- do.call(rbind, r2_values)
colnames(r2_df) <- c("Marginal R2", "Conditional R2")

# Add model names as a column
r2_df <- cbind(Model = model_names, r2_df)

# Print the result
print(r2_df)

write.csv(r2_df, "./tables/v2_JAE/r2_distance_log_JAE.csv")


###---------------------
## Figure 3d
###----------------------
distance_mod_df <- ggpredict(distance_models[[2]], terms = c("type2_scaled[all]", "season"))
distance_mod_df$type2 <- round(distance_mod_df$x * distance.age.scale + distance.age.center)

ptest <- ggplot(data = distance_mod_df) +
  geom_line(aes(x = type2, y = predicted, color = group),
            linewidth = 0.7) +
  ylim(0, 1500) +
  xlim(2, 8) +
  stat_summary(data = distance_mod_df,aes(x = type2, y = predicted, color = group),fun.data = "mean_sdl",fun.args = list(mult = 1),
    geom = "pointrange",size = 0.5,position = position_dodge(0.3)) +
  geom_errorbar(data = distance_mod_df,
                aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
                width = 0.1, position = position_dodge(width = 0.3),
                linewidth = 0.6) +  
  geom_point(data = data2,
             aes(y = distance_to_previous_km, x = type2_num, colour = season),
             size = 2, alpha = 0.1,
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.2)) +
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e","Spring-Summer" = "cornflowerblue"
  )) +scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e","Spring-Summer" = "cornflowerblue"
  )) +
  labs(title = "", x = "", y = "Dist. from previous\n year's principal site") +
  theme_classic(base_size = 14) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    axis.line         = element_blank(),
    axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    text              = element_text(face = "plain"),
    axis.text         = element_text(size = 13, face = "plain"),
    axis.text.x       = element_text(size = 13, face = "plain"),
    axis.text.y       = element_text(size = 13, face = "plain"),
    axis.title        = element_text(size = 13, face = "plain"),
    axis.title.y      = element_text(size = 13, face = "plain"),
    legend.position   = "none",
    legend.box.margin = margin(t = 0.02, r = 0, b = 0, l = 0, unit = "cm"),
    legend.text       = element_text(size = 14, face = "plain"),
    strip.text        = element_blank(),
    strip.background  = element_blank()
  ) +
  guides(
    alpha = "none",
    size  = guide_legend(nrow = 2)
  )

ptest

# Fit the model
mfull = glmmTMB(distance_to_previous_km ~  poly(type2_scaled, 2, raw = FALSE) * season + (1|dev)+(1|colony)+(1|year),
                #family=nbinom2,
                family = Gamma(log),
                data = data2_clean_no_out)

tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)
summary(mfull)
library(sjPlot)
tab_model(mfull, digits = 2)


## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)
# Run the outlier test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outlier test results
print(outlier_test)
testDispersion(mfull) #ok
check_collinearity(mfull) 
qqmath(mfull, id=0.05) #

# Extract residuals and fitted values
residuals <- residuals(mfull, type = "pearson")  # Use Pearson residuals
fitted_values <- fitted(mfull)

### select the best model using dredge

d1 <- dredge(mfull, rank = "AIC") 

best <- get.models(d1,subset=NA) 
Anova(mfull) 
Anova(best[[1]]) 
selection.table <- model.sel(d1) 
selection.table #

# Convert to data.frame
selection.df <- as.data.frame(selection.table)

library(openxlsx)
write.xlsx(selection.df,file = "tables/forRevAB/dredge_LDistance.prev.xlsx",rowNames = TRUE,overwrite = TRUE)

r2 <- matrix(NA, ncol = 2, nrow = 5)
for (i in 1:5) {
  result <- r.squaredGLMM(best[[i]])
  if (!is.null(result)) {
    r2[i, ] <- result[1, ]
  }
}
r2

write.table(r2,"tables/v2_JAE/r2_dredge_distance_log.csv",quote=FALSE,row.names = TRUE)


# Fit the model, NA values in the data will be automatically excluded ··· to redice model complexity I remove random effects with  variance of zero or close to zero.
mfull = glmmTMB(distance_to_previous_km ~  poly(type2_scaled, 2, raw = FALSE) * season + (1|dev),
                #family=nbinom2,
                family = Gamma(log),
                data = data2_clean_no_out)


sim_res <- simulateResiduals(mfull)
plot(sim_res)
res <- residuals(mfull, type = "pearson")
plot(res ~ fitted(mfull))
abline(h = c(-3, 3), col = "red", lty = 2)
summary(mfull)$varcor


tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)
summary(mfull)
library(sjPlot)
tab_model(mfull, digits = 2)

rm(list=ls())
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

###
###Dist of principal sites to the natal colony
###------------------------------------------------

data2<-read.csv("./Data/AB/longest.stage_AB.csv", header = TRUE, dec = ".", sep=",") #second attempt

library(performance)
data2 <- data2 %>%
  filter(!is.na(distance_to_previous_km) & !is.na(prev_year_lat) & !is.na(prev_year_lon) & !is.na(prev_year))

# Calculate distance of bird to center point of colony
data2$dist.to.colony <- deg.dist(lat1=data2$colony.lat,long1=data2$colony.long,long2=data2$mean.long,lat2=data2$mean.lat)
data2$type2_scaled <- scale(data2$type2_num)
data2$dist.to.colony_scaled <- scale(data2$dist.to.colony)
### #re-scale variables
dist.col.age.scale <- attr(data2$type2_scaled,"scaled:scale")
dist.col.age.center <- attr(data2$type2_scaled,"scaled:center")

func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model

dist.col_models <- list()
dist.col_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("dist.to.colony ~ season * ", functional_relationships[[i]], 
                         " + (1|dev)++ (1|colony)+ (1|year) ", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F, 
                 #family = poisson,
                 family = Gamma(log),
                 data = data2)
  
  dist.col_models[[i]] <- mod
  dist.col_mod_summary[[i]] <- summary(mod)
  
  rm(mod)
  rm(model_formula)
  
}

# Compare the performance of all models, including the categorical age model
dist.col_comparison <- compare_performance(dist.col_models[[1]], dist.col_models[[2]], 
                                           dist.col_models[[3]],dist.col_models[[4]])  # Add the new model
dist.col_comparison$Model <- func_rel_names

# Display comparison results
print(dist.col_comparison)

# Order models by AIC (from lowest to highest)
dist.col_comparison <- dist.col_comparison %>% arrange(AIC)
dist.col_comparison$Model <- sapply(dist.col_comparison$Model, as.character)

write.xlsx(dist.col_comparison, "./tables/forRevAB/AIC_dist.col.xls")
# List of model names
model_names <- c("Model 1: Linear", "Model 2: Quadratic ", "Model 3: Exponential", "Model 4: poly")

# Compute R2 for each model
r2_values <- lapply(distance_models, r.squaredGLMM)

# Combine the results into a data frame
r2_df <- do.call(rbind, r2_values)
colnames(r2_df) <- c("Marginal R²", "Conditional R²")

# Add model names as a column
r2_df <- cbind(Model = model_names, r2_df)

# Print the result
print(r2_df)

write.csv(r2_df, "./tables/v2_JAE/r2_distance_log_JAE.csv")

## Fig. 3c
dist.col_mod_df <- ggpredict(dist.col_models[[1]], terms = c("type2_scaled[all]", "season"))
dist.col_mod_df$type2 <- round(dist.col_mod_df$x * dist.col.age.scale + dist.col.age.center)

dist.col.age.scale <- attr(data2$type2_scaled,"scaled:scale")
dist.col.age.center <- attr(data2$type2_scaled,"scaled:center")

p4 <- ggplot(data = dist.col_mod_df) +
  geom_line(aes(x = type2, y = predicted, color = group),
            linewidth = 0.7) +
  stat_summary(data = dist.col_mod_df,aes(x = type2, y = predicted, color = group),
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "pointrange",
    size = 0.5,
    position = position_dodge(0.3)
  ) +
  geom_errorbar(data = dist.col_mod_df,
                aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
                width = 0.1, position = position_dodge(width = 0.3),
                linewidth = 1) +
  ylim(0, 1500) +
  geom_point(data = data2,
             aes(y = dist.to.colony, x = type2_num, colour = season),
             size = 2, alpha = 0.1,
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.2)) +
  
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e","Spring-Summer" = "cornflowerblue"
  )) +
  scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e","Spring-Summer" = "cornflowerblue"
  )) +
  labs(title = "", x = "", y = "Dist. of principal site\n to natal colony") +
  theme_classic(base_size = 14) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    axis.line         = element_blank(),
    axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    text              = element_text(face = "plain"),
    axis.text         = element_text(size = 13, face = "plain"),
    axis.text.x       = element_text(size = 13, face = "plain"),
    axis.text.y       = element_text(size = 13, face = "plain"),
    axis.title        = element_text(size = 13, face = "plain"),
    axis.title.y      = element_text(size = 13, face = "plain"),
    legend.position   = "none",
    legend.box.margin = margin(t = 0.02, r = 0, b = 0, l = 0, unit = "cm"),
    legend.text       = element_text(size = 14, face = "plain"),
    strip.text        = element_blank(),
    strip.background  = element_blank()
  ) +
  guides(
    alpha = "none",
    size  = guide_legend(nrow = 2)
  )

p4

####

# Fit the model
mfull = glmmTMB(dist.to.colony ~ type2_scaled* season + (1|dev)+(1|year)+(1|colony),
                #family=nbinom2,
                family = Gamma(log),
                data = data2)

summary(mfull)
tmp = simulateResiduals(mfull)
plot(tmp)
# Run the outlier test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outlier test results
print(outlier_test) 

testDispersion(mfull) #ok
overdisp_fun(mfull) #ok

check_collinearity(mfull2)

# Extract residuals and fitted values
residuals <- residuals(mfull, type = "pearson")  # Use Pearson residuals
fitted_values <- fitted(mfull)

# Create the residual plot
plot(fitted_values, residuals,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residual Plot",
     pch = 19, col = "blue")
abline(h = 0, col = "red", lwd = 2) 

### select the best model using dredge
options(na.action = "na.fail")  
d1 <- dredge(mfull, rank = "AIC",  
             trace = TRUE)
sw(d1)
summary(d1)
get.models(d1, subset = TRUE)
best <- get.models(d1,subset=NA) 
Anova(mfull) 
Anova(best[[1]]) 
selection.table <- model.sel(d1) 
selection.table #
# Convert to data.frame
selection.df <- as.data.frame(selection.table)

library(openxlsx)
write.xlsx(selection.df,file = "tables/forRevAB/dredge_distnatalcolony.xlsx",rowNames = TRUE,overwrite = TRUE)

# ### get R con and R marginal
# Pre-allocate matrix: rows = number of models, cols = R2m, R2c
r2 <- matrix(NA, ncol = 2, nrow = length(best))
colnames(r2) <- c("R2m", "R2c")

for (i in seq_along(best)) {
  result <- r.squaredGLMM(best[[i]])
  r2[i, ] <- result[1, ]
}

r2
 
############################################
####Distance to principal autumn-winter sites to subsequent summer principal site 
###############################
data2<-read.csv("./Data/AB/longest.stage_AB.csv", header = TRUE, dec = ".", sep=",") #second attempt
sum(duplicated(data2$dist.segment))
str(data2)

aw <- data2 %>%
  dplyr::filter(season == "Autumn-Winter") %>%
  dplyr::select(dev, year, type2, mean.lat, mean.long, dist.segment,colony, sex)

ss <- data2 %>%
  dplyr::filter(season == "Spring-Summer") %>%
  dplyr::select(dev, year,type2, mean.lat, mean.long, dist.segment,colony,sex)

ss <- ss %>%
  dplyr::mutate(year = year - 1)
pairs <- aw %>%
  inner_join(ss, by = c("dev", "year"), suffix = c("_aw", "_ss"))


pairs <- pairs %>%
  mutate(dist_km = distHaversine(
      cbind(mean.long_aw, mean.lat_aw),
      cbind(mean.long_ss, mean.lat_ss)
    ) / 1000
  )

pairs$type2_aw <- as.numeric(gsub("yo", "", pairs$type2_aw))  
pairs$type2_scaled <- scale(pairs$type2_aw)
pairs$distance_scaled <- scale(pairs$dist_km)

### #re-scale variables

# Scaling for distance
distance.age.scale <- attr(pairs$type2_scaled,"scaled:scale")
distance.age.center <- attr(pairs$type2_scaled,"scaled:center")
func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")
functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model


distance_models <- list()
distance_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  # # Create the model formula
  model_formula <- paste("dist_km  ~ ", functional_relationships[[i]],
                         " + (1|dev)+ (1|colony_ss)", sep ="")
  
  # Fit the model using glmmTMB
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F, 
                 family=  Gamma(log),
                 # for model comparison
                 data = pairs)
  
  # Store the model and summary
  distance_models[[i]] <- mod
  distance_mod_summary[[i]] <- summary(mod)
  
  # Clean up the model and formula variables to free memory
  rm(mod)
  rm(model_formula)
  
}

# Compare the performance of all models, including the categorical age model
distance_comparison <- compare_performance(distance_models[[1]], distance_models[[2]], 
                                           distance_models[[3]],distance_models[[4]])  # Add the new model
distance_comparison$Model <- func_rel_names
# Display comparison results
print(distance_comparison)

# Order models by AIC (from lowest to highest)
distance_comparison <- distance_comparison %>% arrange(AIC)
distance_comparison$Model <- sapply(distance_comparison$Model, as.character)

write.xlsx(distance_comparison, "./tables/forRevAB/AIC_distancemainsites_AB.csv")

# List of model names
model_names <- c("Model 1: Exponential", "Model 2: Linear", "Model 3: Quadratic", "Model 4: poly")

# Compute R2 for each model
r2_values <- lapply(distance_models, r.squaredGLMM)

# Combine the results into a data frame
r2_df <- do.call(rbind, r2_values)
colnames(r2_df) <- c("Marginal R2", "Conditional R2")

# Add model names as a column
r2_df <- cbind(Model = model_names, r2_df)

# Print the result
print(r2_df)


###---------------------
## Figure 3e
###----------------------

distance_mod_df <- ggpredict(distance_models[[1]], terms = c("type2_scaled[all]"))
distance_mod_df$type2 <- round(distance_mod_df$x * distance.age.scale + distance.age.center)

ptest3 <- ggplot(data = distance_mod_df) +
  
  geom_line(aes(x = type2, y = predicted),color = "black", linewidth = 0.7) +
  stat_summary(data = distance_mod_df,aes(x = type2, y = predicted),fun.data = "mean_sdl",
               fun.args = list(mult = 1),geom = "pointrange",size = 0.5,position = position_dodge(0.3),
               width = 0) +
  geom_errorbar(data = distance_mod_df,aes(ymin = conf.low, ymax = conf.high, x = type2),
                width = 0,position = position_dodge(width = 0.3),linewidth = 1) +
  geom_point(data = pairs,aes(y = dist_km, x = type2_aw),size = 2, alpha = 0.05,
             position = position_jitter(width = 0.04)
  ) +
  labs(title = "", x = "Age", y = "Dist. between winter \nand subsequent summer sites") +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    axis.line         = element_blank(),
    axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    text              = element_text(face = "plain"),
    axis.text         = element_text(size = 13, face = "plain"),
    axis.text.x       = element_text(size = 13, face = "plain"),
    axis.text.y       = element_text(size = 13, face = "plain"),
    axis.title        = element_text(size = 13, face = "plain"),
    axis.title.x      = element_text(size = 13, face = "plain"),
    axis.title.y      = element_text(size = 13, face = "plain"),
    legend.position   = "none",
    legend.text       = element_text(size = 14, face = "plain"),
    strip.text        = element_blank(),
    strip.background  = element_blank()
  ) +
  
  guides(
    alpha = "none",
    size  = guide_legend(nrow = 2)
  )

ptest3

mfull = glmmTMB(dist_km ~  type2_scaled + (1|dev)+(1|colony_ss),
                #family=nbinom2,
                family = Gamma(log),
                data = pairs) ### 

summary(mfull)
sim_res <- simulateResiduals(mfull)
plot(sim_res)
res <- residuals(mfull, type = "pearson")
simulationOutput <- simulateResiduals(fittedModel = mfull)
testDispersion(simulationOutput)
plot(simulationOutput)

tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)
summary(mfull)
library(sjPlot)
tab_model(mfull, digits = 2)
## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)###qq plot residuals is a bit skewd but not too dramatic
# Run the outlier test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outlier test results
print(outlier_test) ### 


testDispersion(mfull) #ok


check_collinearity(mfull) 
qqmath(mfull, id=0.05) #

# Extract residuals and fitted values
residuals <- residuals(mfull, type = "pearson")  # Use Pearson residuals
fitted_values <- fitted(mfull)

# Create the residual plot
plot(fitted_values, residuals,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residual Plot",
     pch = 19, col = "blue")
abline(h = 0, col = "red", lwd = 2)  # Add a horizontal line at y = 0

summary(mfull)

options(scipen=999) 


### select the best model using dredge

d1 <- dredge(mfull, rank = "AIC") 

best <- get.models(d1,subset=NA) 
Anova(mfull) 
Anova(best[[1]]) 
selection.table <- model.sel(d1) 
selection.table #

# Convert to data.frame
selection.df <- as.data.frame(selection.table)
library(openxlsx)
write.xlsx(selection.df,file = "tables/forRevAB/dredge_LDistancewinterandsummer.prev.xlsx",rowNames = TRUE,overwrite = TRUE)

r2 <- matrix(NA, ncol = 2, nrow = 5)
for (i in 1:5) {
  result <- r.squaredGLMM(best[[i]])
  if (!is.null(result)) {
    r2[i, ] <- result[1, ]
  }
}
r2

###print all panels 

top_row <- cowplot::plot_grid(p, p2, p4, ptest,ncol = 2,align = "hv",axis = "tblr")
bottom_row <- cowplot::plot_grid(ptest3, ncol = 1)

final1 <- cowplot::plot_grid(top_row, bottom_row,ncol = 1,rel_heights = c(2, 1))
final1

ggsave("Drafts/AB/third_re-submission/Figures/fig_3_fv.tiff",plot=final1,width=174,height=220,dpi=300)
