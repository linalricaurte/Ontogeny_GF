###function overdispertion
overdisp_fun <- function(model) {
  # Extract residual degrees of freedom
  rdf <- df.residual(model)
  
  # Extract Pearson residuals
  rp <- residuals(model, type = "pearson")
  
  # Calculate the sum of squared Pearson residuals
  Pearson_chisq <- sum(rp^2)
  
  # Calculate the dispersion ratio
  dispersion <- Pearson_chisq / rdf
  
  return(dispersion)
}

stop.segs<-read.csv('./Data/stop.segs.csv', header = TRUE, dec = ".", sep=",") #I create this file from the 01_fig_1_tracking_data_JAE

stop.segs$segment_type <- ifelse(stop.segs$bout.dur > 3 & stop.segs$ bout.dist.cum  > 10, "long.stage", "short.stage")

l.stage<-subset(stop.segs, stop.segs$segment_type=="long.stage")

summary.per.birdstop <- l.stage %>%
  group_by(dev, type2, season2,sex, colony) %>%  # Group by individual and age class
  dplyr::summarise(n.stops = n()) %>%         # Count the number of stop events per individual and age class
  ungroup()

unique_individuals <- summary.per.birdstop%>%
  group_by(season2) %>%
  dplyr::summarise(unique_count = sum(n.stops))### 

#write.csv(summary.per.birdstop, "./tables/v2_JAE/summary.per.birdstop_JAE.csv")


l.stage.for.mod<-summary.per.birdstop

###--------------------------
###Number of resident events###
###---------------------------

#Scale and center 'type2'(age) around its mean and standard deviation
### first, make type 2 numeric
l.stage.for.mod$type2 <- as.numeric(gsub("yo", "", l.stage.for.mod$type2))  # Remove 'yo' and convert to numeric
l.stage.for.mod$type2_scaled <- scale(l.stage.for.mod$type2)


### #re-scale variables
events.age.scale <- attr(l.stage.for.mod$type2_scaled,"scaled:scale")
events.age.center <- attr(l.stage.for.mod$type2_scaled,"scaled:center")
events.age.scale <- attr(l.stage.for.mod$type2_scaled,"scaled:scale")
events.age.center <- attr(l.stage.for.mod$type2_scaled,"scaled:center")


###############################################################################
# Note: The code below (lines 3-113) to investigate the effect of age on staging behaviour
# was adapted from:
# Marta Acácio (2025). Ageing Vultures. GitHub repository:
# https://github.com/msa2015/Ageing_Vultures
# Modifications were made to meet specific project needs.
###############################################################################

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # 

func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")


events_models <- list()

events_mod_summary <- list()

# Loop through all 4 functional relationships, 
for(i in 1:length(functional_relationships)){
  
  # Construct the model formula dynamically
  model_formula <- paste("n.stops ~ season2 + ", functional_relationships[[i]], 
                         " + (1|dev)", sep = "")
  
  # Fit the model
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                              
                 family = poisson,  
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

# Display comparison results
print(events_comparison)

###-------------------------
### build figure 2a ####
###-------------------------

# Predict values for the selected model (e.g., third-degree polynomial or age as a factor)
events_mod_df <- ggpredict(events_models[[1]], terms = c("type2_scaled[all]", "season2"))
events_mod_df$type2 <- round(events_mod_df$x * events.age.scale + events.age.center)

p <- 
  ggplot(data = events_mod_df) +
  geom_line(aes(x = type2, y = predicted, color = group), 
            linewidth = 0.7) +
  stat_summary(
    data = events_mod_df,
    aes(x = type2, y = predicted, color = group),
    fun.data = "mean_sdl",
    fun.args = list(mult = 1), 
    geom = "pointrange",
    size = 0.5,
    position = position_dodge(0.3)
  ) +
  geom_errorbar(data = events_mod_df, aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
                width = 0.1, position = position_dodge(width = 0.3), size = 1) +
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  labs(
    title = "",
    x = "", 
    y = "Number of \nstaging events"
  ) +
  geom_jitter(data = l.stage.for.mod,
              aes(y = n.stops, x = type2, colour = season2, group = type2),
              size = 1, alpha = 0.2, width = 0.1) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_blank(),
    panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title = element_text(size = 13, face = 'bold'),
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none",  # Display legend
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12, face = 'bold'),
    strip.background = element_blank()
  ) +
  guides(alpha = FALSE, size = guide_legend(nrow = 2))
p

# get stadistics

#extract means for each season
means <- aggregate(l.stage.for.mod$n.stops,by=list(l.stage.for.mod$season2),mean)
# #Now same for sd's
sds <- aggregate(l.stage.for.mod$n.stops,by=list(l.stage.for.mod$season2),sd)



# Fit the model with   n .stops using a poisson familly. 
## poisson family fit better the data
###first full model with interaction between age and season


### the most parsimonious/final model according to model selection (dredge) is without the interaction 
mfull <- glmmTMB(n.stops ~ type2_scaled * season2 + (1|dev),
               family = poisson,
               data = l.stage.for.mod)

summary(mfull)

mfull <- glmmTMB(n.stops ~ type2_scaled + season2 + (1|dev),
                 family = poisson,
                 data = l.stage.for.mod)

summary(mfull)

tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)


###############
##########

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


testDispersion(mfull) #1.37, ok
overdisp_fun(mfull) #1.10, ok


check_collinearity(mfull)##ok
diagnostics.plot(mfull)


### select the best model using dredge

d1 <- dredge(mfull, rank = "AIC") # 

best <- get.models(d1,subset=NA) #get full model set
Anova(mfull) # anova for full model
Anova(best[[1]]) # run anova on best model from set
selection.table <- model.sel(d1) # pull out model selection table
selection.table #
write.table(selection.table,"tables/v2_JAE/dredge_stops_poisson4d.csv",quote=FALSE,row.names = TRUE)

### get R con and R marginal

r2 <- matrix(NA, ncol = 2, nrow = 5)  # Pre-allocate matrix
for (i in 1:5) {
  result <- r.squaredGLMM(best[[i]])
  if ("delta" %in% rownames(result)) {
    r2[i, ] <- result["delta", ]
  } else {
    warning(paste("Model", i, "does not include 'delta'."))
  }
}

write.table(r2,"tables/v2_JAE/r2_dredge_stops_poisson4d.csv",quote=FALSE,row.names = TRUE)

#############################################
###---------Average resident duration in days
#############################################

l.stage<-subset(stop.segs, stop.segs$segment_type=="long.stage")



### compute the mean duration per device, age and season
 l.stage <- l.stage %>%
   group_by(dev, type2, season2) %>%
   dplyr::summarise(
     mean.bout.dur = mean(bout.dur, na.rm = TRUE),
     year = first(year),   # Assuming year is consistent for each dev, type2, season2
     colony = first(colony),  # Assuming colony is consistent for each dev, type2, season2
     .groups = "drop"
   )


l.stage$type2 <- as.numeric(gsub("yo", "", l.stage$type2))  # Remove 'yo' and convert to numeric
l.stage$type2_scaled <- scale(l.stage$type2)


### Compute some summary stadistics
# #extract means for each season
means <- aggregate(l.stage$mean.bout.dur,by=list(l.stage$season2),mean)

# 
# #Now same for sd's
sds <- aggregate(l.stage$mean.bout.dur,by=list(l.stage$season2),sd)



## to get the same scale of analysis in figure 2 I will compute the mean the average staging duration per individual per year per season
write.csv(l.stage,"Data/Av_res_dur.csv",quote=FALSE,row.names = TRUE)


### #re-scale variables

duration.age.scale <- attr(l.stage$type2_scaled,"scaled:scale")
duration.age.center <- attr(l.stage$type2_scaled,"scaled:center")


func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model


duration_models <- list()
duration_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("mean.bout.dur ~ season2 * ", functional_relationships[[i]], 
                         " + (1|dev)+(1|year)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,
                 family=nbinom1,
                 # for model comparison
                 #family = poisson,
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

write.table(duration_comparison, "./tables/v2_JAE/AIC_duration_log_3dJAE.csv")

## Figure
duration_mod_df <- ggpredict(duration_models[[1]], terms = c("type2_scaled[all]", "season2"))
duration_mod_df$type2 <- round(duration_mod_df$x * duration.age.scale + duration.age.center)
#duration_mod_df <- rename(duration_mod_df, season2 = group)



###--------------------------
###figure 2b################
###--------------------------
p2 <- ggplot(data = duration_mod_df) +
  #Remove facet_grid and plot both groups in the same plot
  # geom_ribbon(aes(x = type2, y = predicted,
  #                 ymin = conf.low, ymax = conf.high, fill = group),
  #             alpha = 0.2, linewidth = 0.6) +
  geom_line(aes(x = type2, y = predicted, color = group),
            linewidth = 0.7) +

  stat_summary(
    data = duration_mod_df,
    aes(x = type2, y = predicted, color = group),
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "pointrange",
    size = 0.5,
    position = position_dodge(0.3)
  ) +
  geom_errorbar(data = duration_mod_df, aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
                width =0.1, position = position_dodge(width = 0.3), size=1) +
  #facet_grid(.~group)+
  #ylim(0,150)+
  geom_jitter(data=l.stage, aes(y = mean.bout.dur, x = type2, colour=season2, group=type2), size=1,alpha=.1,width = 0.1) +
  #Custom colors for the groups
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  labs(
    title = "",  # Add your title here
    x = "", 
    y = "Average staging\n duration (days)"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        axis.line = element_line(colour = "black"),
        legend.position   = "none",
        #legend.direction  = 'vertical',
        legend.box.margin = margin(t=0.02,r=0,b=0,l=0,unit="cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  
  # Removing legend if you don't want it
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

p2

####


# Fit the full model as before
mfull = lmer(mean.bout.dur ~ type2_scaled  *season2 + (1|dev) + (1|year) + (1|colony),
             data = l.stage)


# ### model selection selects season and age (not the interation) 
# mfull = glmmTMB(mean.bout.dur ~ expm1(type2_scaled)  + season2 + (1|dev) + (1|year) + (1|colony),
#                 data = l.stage)


tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)
summary(mfull)

## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)
# Run the outliers test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outliers test results
print(outlier_test) ### not extreme outliers


testDispersion(mfull) #1.07, ok
#overdisp_fun(mfull) #model with mean bout duration is overdispersed, 1219 with log(mean.bout.dur) I deal with overdispersion (0.35, ok)
check_overdispersion(mfull)


check_collinearity(mfull) ###ok

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

### select the best model using dredge

d1 <- dredge(mfull, rank = "AIC") # Do not include REML here

best <- get.models(d1,subset=NA) #get full model set
Anova(mfull) # anova for full model
Anova(best[[1]]) # run anova on best model from set
selection.table <- model.sel(d1) # pull out model selection table
selection.table #
write.table(selection.table,"tables/v2_JAE/dredge_duration_logv2.csv",quote=FALSE,row.names = TRUE)


r2 <- matrix(NA, ncol = 2, nrow = 5)  # Pre-allocate matrix
colnames(r2) <- c("R2m", "R2c")      # Marginal and Conditional R²
for (i in 1:5) {
  result <- r.squaredGLMM(best[[i]])
  if (!is.null(result)) {
    r2[i, ] <- result[1, ]  # Assuming the first row contains the R² values
  } else {
    warning(paste("Model", i, "did not return R² values."))
  }
}

write.table(r2,"tables/v2_JAE/r2_dredge_duration_log.csv",quote=FALSE,row.names = TRUE)

###########################################
### -------distance to principal sites 
###############################
#data2<-read.csv("./Data/long.stage_REV_JAE.csv", header = TRUE, dec = ".", sep=",") #I create this txt using script 02_Fig_stage_events_dur_num_JAE

data2<-read.csv("./Data/long.stage_REV6_JAE.csv", header = TRUE, dec = ".", sep=",") #second attempt


#data2$distance_to_previous_km <- ifelse(is.na(data2$distance_to_previous_km) == TRUE,0,data2$distance_to_previous_km)

#remove nas
data2_clean <- data2 %>%
  filter(!is.na(distance_to_previous_km) & !is.na(prev_year_lat) & !is.na(prev_year_lon) & !is.na(prev_year))


data2_clean$type2 <- as.numeric(gsub("yo", "", data2_clean$type2))  # Remove 'yo' and convert to numeric
data2_clean$type2_scaled <- scale(data2_clean$type2_num)
data2_clean$distance_scaled <- scale(data2_clean$distance_to_previous_km)
write.csv(data2_clean, "./tables/v2_JAE/dist_prin_sites.csv")

### #re-scale variables

# Scaling for distance
distance.age.scale <- attr(data2_clean$type2_scaled,"scaled:scale")
distance.age.center <- attr(data2_clean$type2_scaled,"scaled:center")


func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model


library(MuMIn)    # for calculating R²
distance_models <- list()
distance_mod_summary <- list()


for(i in 1:length(functional_relationships)){
  
  
  # Create the model formula
  model_formula <- paste("distance_to_previous_km ~ season2 * ", functional_relationships[[i]], 
                         " + (1|dev)", sep ="")
  
  
  # Fit the model using glmmTMB
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F, 
                 family=  Gamma(log),
                 # for model comparison
                 data = data2_clean)
  
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

write.table(distance_comparison, "./tables/v2_JAE/AIC_distance_4dJAE.csv")

# List of model names
model_names <- c("Model 1: Exponential", "Model 2: Linear", "Model 3: Quadratic", "Model 4: poly")

# Compute R² for each model
r2_values <- lapply(distance_models, r.squaredGLMM)

# Combine the results into a data frame
r2_df <- do.call(rbind, r2_values)
colnames(r2_df) <- c("Marginal R²", "Conditional R²")

# Add model names as a column
r2_df <- cbind(Model = model_names, r2_df)

# Print the result
print(r2_df)

write.csv(r2_df, "./tables/v2_JAE/r2_distance_log_JAE.csv")


###---------------------
## Figure 2D
###----------------------

distance_mod_df <- ggpredict(distance_models[[1]], terms = c("type2_scaled[all]", "season2"))
distance_mod_df$type2 <- round(distance_mod_df$x * distance.age.scale + distance.age.center)
#duration_mod_df <- rename(duration_mod_df, season2 = group)

ptest <- ggplot(data = distance_mod_df) +
  #Remove facet_grid and plot both groups in the same plot
  # geom_ribbon(aes(x = type2, y = predicted,
  #                 ymin = conf.low, ymax = conf.high, fill = group),
  #             alpha = 0.2, linewidth = 0.6) +
  geom_line(aes(x = type2, y = predicted, color = group),
            linewidth = 0.7) +
  ylim(0,1500)+
  xlim(2,8)+
  stat_summary(
    data = distance_mod_df,
    aes(x = type2, y = predicted, color = group),
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "pointrange",
    size = 0.5,
    position = position_dodge(0.3)
  ) +
  geom_errorbar(data = distance_mod_df, aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
                width =0.1, position = position_dodge(width = 0.3), size=1) +
  #facet_grid(.~group)+

  geom_jitter(data=data2, aes(y =distance_to_previous_km, x = type2_num, colour=season2, group=type2_num), size=1,alpha=.1,width = 0.1) +
  #Custom colors for the groups
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  labs(
    title = "",  # Add your title here
    x = "Age", 
    y = "Dist. from previous\n main sites (km)"
  ) +
  
  # Y-axis label
  #scale_y_continuous(name = "Staging duration\n (days)") +
  #xlab("Age (years)") +  # Updated x-axis label for clarity
  # Customizing the theme
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 13),
        axis.text = element_text(size = 13),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        axis.line = element_line(colour = "black"),
        legend.position   = "none",
        #legend.direction  = 'vertical',
        legend.box.margin = margin(t=0.02,r=0,b=0,l=0,unit="cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  
  # Removing legend if you don't want it
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

ptest


####estimates of teh linear model
# Fit the model as you've done
# Remove rows with missing values in the relevant columns


# Fit the full model with cleaned data

# mfull = lmer(distance_scaled ~  poly(type2_scaled, 3, raw = TRUE) * season2 + (1|dev) + (1|year) + (1|colony),
#              data = data2_clean)

# mfull = lmer(distance_to_previous_km ~  type2_scaled * season2 + (1|dev) + (1|year),
#               data = data2_clean)

### after model selection, the best model includes the interaction 
# mfull = lmer(log(distance_to_previous_km) ~  poly(type2_scaled, 3, raw = TRUE) * season2 + (1|dev) + (1|year) + (1|colony),
#              data = data2_clean)

unique_individuals <- data2_clean %>%
  group_by(dev) %>%
  dplyr::summarise(unique_count = n_distinct(dev))


# Fit the model, NA values in the data will be automatically excluded ··· to redice model complexity I remove random effects with  variance of zero or close to zero.
mfull = glmmTMB(distance_to_previous_km ~ type2_scaled * season2 + (1|dev),
                #family=nbinom2,
                family = Gamma(log),
             data = data2_clean)
outliers <- check_outliers(mfull, method = "cook")
outliers #ok

r2_values <- performance::r2(mfull)
print(r2_values)

check_overdispersion(mfull) #ok


# Print the model summary
summary(mfull)
summ(mfull,ddf = "Satterthwaite", digits=2)
#mfull = lmer(log(distance_km) ~type2_scaled *season2 + (1|dev) + (1|year) + (1|colony),
 #            data = data2)


# poly(type2_scaled, 2, raw = TRUE)",
#                                  "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")
# 
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
print(outlier_test) ### not extreme aoutliers


testDispersion(mfull) #0.93, ok


check_collinearity(mfull) ## MODERATE
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

d1 <- dredge(mfull, rank = "AIC") # Do not include REML here

best <- get.models(d1,subset=NA) #get full model set
Anova(mfull) # anova for full model
Anova(best[[1]]) # run anova on best model from set
selection.table <- model.sel(d1) # pull out model selection table
selection.table #
write.table(selection.table,"tables/v2_JAE/dredge_distance_cons.csv",quote=FALSE,row.names = TRUE)


r2 <- matrix(NA, ncol = 2, nrow = 5)  # Pre-allocate matrix
#colnames(r2) <- c("R2m", "R2c")      # Marginal and Conditional R²
for (i in 1:5) {
  result <- r.squaredGLMM(best[[i]])
  if (!is.null(result)) {
    r2[i, ] <- result[1, ]  # Assuming the first row contains the R² values
  } else {
    warning(paste("Model", i, "did not return R² values."))
  }
}

write.table(r2,"tables/v2_JAE/r2_dredge_distance_log.csv",quote=FALSE,row.names = TRUE)

###----------------------------
####### figure 2E #############
###----------------------------

unique_individuals <- l.stage.for.mod %>%
  group_by(dev, type2, season2) %>%
  dplyr::summarise(unique_count = n_distinct(dev))



ordered_segments <- sort(unique(l.stage.for.mod$n.stops))
print(ordered_segments)
l.stage.for.mod$n.segments.cat <- recode_factor(l.stage.for.mod$n.stops, 
                                         "1" = "1-3", 
                                         "2" = "1-3",
                                         "3" = "1-3",
                                         "4" = "4-6",
                                         "5" = "4-6",
                                         "6" = "4-6",
                                         "7" = "7-9",
                                         "8" = "7-9", 
                                         "9" = "7-9",
                                         "10" = ">10",
                                         "11" = ">10",
                                         "12" = ">10",
                                         "13" = ">10",
                                         "14" = ">10",
                                         "15" = ">10")


#dev.off()

pc <- l.stage.for.mod %>%
  ggplot(aes(x = factor(type2), fill = factor(n.segments.cat))) +
  geom_bar(position = "fill", alpha = 0.8) +
  scale_fill_viridis(discrete = TRUE) +
  facet_grid(. ~ season2) +
  labs(
    title = "",                # Add your title here
    x = "Age",                 # Add your x-axis name here
    y = "Proportion\nof individuals",
    fill = "Num. of staging events"  # Set legend title here
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13, face = 'bold'),
    axis.line = element_line(colour = "black"),
    legend.position = "none",  # If you want no legend, leave this as "none"
    legend.text = element_text(size = 14),  # Style legend text
    legend.title = element_text(size = 14, face = 'bold'),  # Style legend title
    strip.text = element_text(size = 13, face = 'bold'),
    strip.background = element_rect(fill = "white", colour = NA)
  ) +
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

pc

###------------------------------------------------
###-----max of principal sites  to the natal colony
###------------------------------------------------
#data2<-read.csv("./Data/long.stage_REV_JAE.csv", header = TRUE, dec = ".", sep=",") #I create this txt using script 02_Fig_stage_events_dur_num_JAE

data2<-read.csv("./Data/long.stage_REV6_JAE.csv", header = TRUE, dec = ".", sep=",") #second attempt

library(performance)
data2 <- data2 %>%
  filter(!is.na(distance_to_previous_km) & !is.na(prev_year_lat) & !is.na(prev_year_lon) & !is.na(prev_year))

# Calculate distance of bird to center point of colony
data2$dist.to.colony <- deg.dist(lat1=data2$colony.lat,long1=data2$colony.long,long2=data2$mean.long,lat2=data2$mean.lat)

#data2$type2 <- as.numeric(gsub("yo", "", data2$type2))  # Remove 'yo' and convert to numeric
data2$type2_scaled <- scale(data2$type2_num)
data2$dist.to.colony_scaled <- scale(data2$dist.to.colony)


# str(dist.breed)
### filter some missclassifications
#data2<-subset(data2,data2$bout.dur<200)

### #re-scale variables

dist.col.age.scale <- attr(data2$type2_scaled,"scaled:scale")
dist.col.age.center <- attr(data2$type2_scaled,"scaled:center")


func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model

# # For the models
# functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)", 
#                                  "expm1(type2_scaled)")  # Add the categorical age model
# 
# 
# func_rel_names <- list("Linear", "Quadratic", "Exponential")
# 


###-------- duration ----
dist.col_models <- list()
dist.col_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("dist.to.colony ~ season2 * ", functional_relationships[[i]], 
                         " + (1|dev) ", sep ="")
  
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

write.table(dist.col_comparison, "./tables/v2_JAE/AIC_dist.col_log_3dJAE.csv")


## extract rcon manually
# mfull = lmer(log(dist.to.colony) ~ poly(type2_scaled, 3, raw = TRUE)*season2 + (1|dev) + (1|year) + (1|colony),
#              data = data2)
# 



## Figure
dist.col_mod_df <- ggpredict(dist.col_models[[1]], terms = c("type2_scaled[all]", "season2"))
dist.col_mod_df$type2 <- round(dist.col_mod_df$x * dist.col.age.scale + dist.col.age.center)
#duration_mod_df <- rename(duration_mod_df, season2 = group)
write.csv(dist.col_mod_df, "./tables/v2_JAE/dist.col_mod_JAE.csv")

dist.col.age.scale <- attr(data2$type2_scaled,"scaled:scale")
dist.col.age.center <- attr(data2$type2_scaled,"scaled:center")


####plot
p4 <- ggplot(data = dist.col_mod_df) +
  #Remove facet_grid and plot both groups in the same plot
  # geom_ribbon(aes(x = type2, y = predicted,
  #                 ymin = conf.low, ymax = conf.high, fill = group),
  #             alpha = 0.2, linewidth = 0.6) +
  geom_line(aes(x = type2, y = predicted, color = group),
            linewidth = 0.7) +
  
  stat_summary(
    data = dist.col_mod_df,
    aes(x = type2, y = predicted, color = group),
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "pointrange",
    size = 0.5,
    position = position_dodge(0.3)
  ) +
  geom_errorbar(data = dist.col_mod_df, aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
                width =0.1, position = position_dodge(width = 0.3), size=1) +
  #facet_grid(.~group)+
  ylim(0,1500)+
  geom_jitter(data=data2, aes(y = dist.to.colony, x = type2_num, colour=season2, group=type2_num), size=1,alpha=.1,width = 0.1) +
  #Custom colors for the groups
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  labs(
    title = "",  # Add your title here
    x = "Age", 
    y = "Dist. of main sites\n to natal colony (km)"
  ) +
  
  # Y-axis label
  #scale_y_continuous(name = "Staging duration\n (days)") +
  #xlab("Age (years)") +  # Updated x-axis label for clarity
  # Customizing the theme
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 13),
        axis.text = element_text(size = 13),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        axis.line = element_line(colour = "black"),
        legend.position   = "none",
        #legend.direction  = 'vertical',
        legend.box.margin = margin(t=0.02,r=0,b=0,l=0,unit="cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  
  # Removing legend if you don't want it
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

p4
ptest

####

####estimates of teh linear model
# Fit the model as you've done
mfull = lmer(dist.to.colony ~  type2_scaled*season2 + (1|dev),
             data = data2)

#mfull = lmer(log(distance_to_previous_km) ~ poly(type2_scaled, 3, raw = TRUE)*season2 + (1|dev) + (1|year),
 #            data = data2)

# Fit the model, NA values in the data will be automatically excluded
mfull = glmmTMB(dist.to.colony ~ type2_scaled * season2 + (1|dev),
                #family=nbinom2,
                family = Gamma(log),
                data = data2)

mfull = glmmTMB(dist.to.colony ~ type2_scaled + season2 + (1|dev),
                #family=nbinom2,
                family = Gamma(log),
                data = data2)


#summ(mfull,ddf = "Satterthwaite", digits=3)
tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)
tab_model(mfull, transform = NULL, auto.label = FALSE, collapse.ci=TRUE, file = "./tables/results.xls")
summary(mfull)

## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)
# Run the outlier test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outlier test results
print(outlier_test) ### not extreme aoutliers


testDispersion(mfull) #1.00, ok
overdisp_fun(mfull) #0.87, ok


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




### select the best model using dredge

d1 <- dredge(mfull, rank = "AIC") # Do not include REML here

best <- get.models(d1,subset=NA) #get full model set
Anova(mfull) # anova for full model
Anova(best[[1]]) # run anova on best model from set
selection.table <- model.sel(d1) # pull out model selection table
selection.table #
write.table(selection.table,"tables/v2_JAE/dredge_dist.col_logv2.csv",quote=FALSE,row.names = TRUE)


r2 <- matrix(NA, ncol = 2, nrow = 5)  # Pre-allocate matrix
colnames(r2) <- c("R2m", "R2c")      # Marginal and Conditional R²
for (i in 1:5) {
  result <- r.squaredGLMM(best[[i]])
  if (!is.null(result)) {
    r2[i, ] <- result[1, ]  # Assuming the first row contains the R² values
  } else {
    warning(paste("Model", i, "did not return R² values."))
  }
}

write.table(r2,"tables/v2_JAE/r2_dredge_dist.col_log.csv",quote=FALSE,row.names = TRUE)


# Create the top row with four panels
top_row <- cowplot::plot_grid(p, p2, p4, ptest, labels = c("A", "B", "C", "D"),
                              label_size = 16, ncol = 2, align = "vh", axis = "lr")

top_row
# Combine the top row and the bottom panel
pfin <- cowplot::plot_grid(top_row, pc, labels = c("", "E"), label_size = 16, 
                           ncol = 1, rel_heights = c(2, 1)) # Adjust heights as needed

# Display the final plot
pfin


legend_b1 <- get_legend(
  p + 
    guides(fill = guide_legend(nrow = 1, title = "Season")) +  # Set legend title for fill
    theme(
      legend.position = "bottom",  # Place legend at the bottom
      legend.title.position = "top"  # Ensure the title is at the top of the legend
    )
)

# Create the second legend (legend_b2)
legend_b2 <- get_legend(
  pc + 
    guides(color = guide_legend(nrow = 1, title = "Num. of staging events")) +  # Set legend title for color
    theme(
      legend.position = "bottom",  # Place legend at the bottom
      legend.title.position = "top"  # Ensure the title is at the top of the legend
    )
)


# Combine the legends horizontally
combined_legend <- cowplot::plot_grid(legend_b1, legend_b2, ncol = 2, rel_widths = c(1, 1))

# Combine the main plot grid with the combined legends
pfin2 <- cowplot::plot_grid(pfin, combined_legend, ncol = 1, rel_heights = c(5, 1))
pfin2

ggsave(plot=pfin2,filename='./Figures_2025/Fig2_v12_JAE.tiff',dpi=300,width=8.5,height=9.5)


# # Calculate distance of bird to center point of colony
# source('Rfunctions/sidescript_pt2pt_fxns_v20230911.R')
# 
# # The pt2pt.range function returns distance in km 
# 
# # Calculate distance of bird to center point of colony
# data2$dist.to.colony <- deg.dist(lat1=data2$colony.lat,long1=data2$colony.long,long2=data2$mean.long,lat2=data2$mean.lat)

#-----------------------------------------------------
### check sex differences with the subset of data
#---------------------------------------------------
## subset birds with known age
table(l.stage.for.mod$sex)

l.stage.for.mod.sex<-subset(l.stage.for.mod, sex %in% c ("Male", "Female"))

unique_individuals <- l.stage.for.mod.sex%>%
  group_by(type2, sex) %>%
  dplyr::summarise(unique_count = n_distinct(dev))



#write.csv(unique_individuals, "./tables/v2_JAE/sex_dev_age_JAE.csv")


### I remove age 8 and 7 and 6  because I dint have enogh male and females to make comparisons
## filter adult birds during may and june
l.stage.for.mod.sex <- subset(l.stage.for.mod.sex, type2 %in% c("1", "2", "3", "4"))

#Scale and center 'type2' around its mean and standard deviation
### first, make type 2 numeric
l.stage.for.mod.sex$type2 <- as.numeric(gsub("yo", "", l.stage.for.mod.sex$type2))  # Remove 'yo' and convert to numeric
l.stage.for.mod.sex$type2_scaled <- scale(l.stage.for.mod.sex$type2)

### #re-scale variables

events.age.scale <- attr(l.stage.for.mod.sex$type2_scaled,"scaled:scale")
events.age.center <- attr(l.stage.for.mod.sex$type2_scaled,"scaled:center")
events.age.scale <- attr(l.stage.for.mod.sex$type2_scaled,"scaled:scale")
events.age.center <- attr(l.stage.for.mod.sex$type2_scaled,"scaled:center")


# For the models
# Add a model where 'age' is a factor
functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model


func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

# functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)", 
#                                  "expm1(type2_scaled)")  # Add the categorical age model
# 
# 
# func_rel_names <- list("Linear", "Quadratic", "Exponential")

events_models <- list()
events_mod_summary <- list()

# Loop through all functional relationships, including age as a factor
for(i in 1:length(functional_relationships)){
  
  # Construct the model formula dynamically
  model_formula <- paste("n.stops ~ season2+sex * ", functional_relationships[[i]], 
                         " + (1|dev)", sep = "")
  
  # Fit the model
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                              
                 family = poisson,  
                 data = l.stage.for.mod.sex)
  
  # Store the models and their summaries
  events_models[[i]] <- mod
  events_mod_summary[[i]] <- summary(mod)
  
  # Cleanup for the next iteration
  rm(mod)
  rm(model_formula)
}

# Compare the performance of all models, including the categorical age model

# Compare the performance of all models, including the categorical age model
events_comparison <- compare_performance(events_models[[1]], events_models[[2]], 
                                         events_models[[3]],events_models[[4]])  # Add the new model

events_comparison$Model <- func_rel_names

# Order models by AIC (from lowest to highest)
events_comparison <- events_comparison %>% arrange(AIC)

events_comparison$Model <- sapply(events_comparison$Model, as.character)


write.csv(events_comparison, "./tables/v2_JAE/AIC_stops_poisson_sex_JAE.csv")



# Display comparison results
print(events_comparison)

#components <- get_plot_component(p, "guide-box", return_all = TRUE)
#components


# Predict values for the selected model (e.g., third-degree polynomial or age as a factor)
events_mod_df <- ggpredict(events_models[[2]], terms = c("type2_scaled[all]", "season2", "sex"))
events_mod_df$type2 <- round(events_mod_df$x * events.age.scale + events.age.center)



p <- 
  ggplot(data = events_mod_df) +
  # geom_line(aes(x = type2, y = predicted, color = group),
  #           linewidth = 0.7) +
  facet_wrap(~sex) +
  stat_summary(
    data = events_mod_df,
    aes(x = type2, y = predicted, color = group),
    fun.data = "mean_sdl",
    fun.args = list(mult = 1), 
    geom = "pointrange",
    size = 0.5,
    position = position_dodge(0.3)
  ) +
  geom_errorbar(data = events_mod_df, aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
                width = 0.1, position = position_dodge(width = 0.3), size = 1) +
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  #scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  labs(
    title = "",
    x = "", 
    y = "Number of \nstaging events"
  ) +
  geom_jitter(data = l.stage.for.mod.sex,
              aes(y = n.stops, x = type2, colour = season2, group = type2),
              size = 1, alpha = 0.2, width = 0.1) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_blank(),
    panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title = element_text(size = 13, face = 'bold'),
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none",  # Display legend
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12, face = 'bold'),
    strip.background = element_blank()
  ) +
  guides(alpha = FALSE, size = guide_legend(nrow = 2))
p

# get stadistics

#extract means for each season
means <- aggregate(l.stage.for.mod$n.stops,by=list(l.stage.for.mod$type2,l.stage.for.mod$season2),mean)
means <- aggregate(l.stage.for.mod$n.segments,by=list(l.stage.for.mod$type2,l.stage2$season2),mean)

winter <- subset(means,Group.2=="Autumn-Winter")
summer <- subset(means,Group.2=="Spring-Summer")

#Now same for sd's
sds <- aggregate(l.stage2$n.segments,by=list(l.stage2$age_gp,l.stage2$season2),sd)
sdwinter <- subset(sds,Group.2=="Autumn-Winter")
sdsummer <- subset(sds,Group.2=="Spring-Summer")


# Fit the model with log transformed n.stops and a normal distribution and n .stops using a poisson familly. 
## poisson familly fit bettwer the data
# mfull <- lmer(log(n.stops) ~ type2_scaled * season2 + (1|dev) + (1|year) + (1|colony),
#               data = l.stage.for.mod, REML = FALSE) # Set REML to FALSE here

# mfull1 <- glmer(n.stops ~ sex*type2_scaled + (1|dev) + (1|year) + (1|colony),
#                family = poisson,
#                data = l.stage.for.mod.sex)
mfull1 <- glmmTMB(n.stops ~ sex*poly(type2_scaled, 2, raw = TRUE) + season2 + (1|dev),
                 family = poisson,
                 data = l.stage.for.mod.sex)




tab_model(mfull1, show.se = TRUE, show.stat = TRUE, transform = NULL) ### no sex differences

## check model assumtions
tmp = simulateResiduals(mfull1)
plot(tmp)
summary(mfull1)

testOutliers(tmp, type = "bootstrap")
testDispersion(tmp)

# Run the outlier test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outlier test results
print(outlier_test) ### not extreme aoutliers


testDispersion(mfull1) #0.83, ok
overdisp_fun(mfull1) #0.77, ok


check_collinearity(mfull1)
diagnostics.plot(mfull1)
qqmath(mfull, id=0.05) #

# Extract residuals and fitted values
residuals <- residuals(mfull, type = "pearson")  # Use Pearson residuals
fitted_values <- fitted(mfull1)

# Create the residual plot
plot(fitted_values, residuals,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residual Plot",
     pch = 19, col = "blue")
abline(h = 0, col = "red", lwd = 2)  # Add a horizontal line at y = 0



### select the best model using dredge

d1 <- dredge(mfull1, rank = "AIC") # Do not include REML here

best <- get.models(d1,subset=NA) #get full model set
Anova(mfull1) # anova for full model
Anova(best[[1]]) # run anova on best model from set
selection.table <- model.sel(d1) # pull out model selection table
selection.table #
write.table(selection.table,"tables/v2_JAE/dredge_stops_poisson_sex.csv",quote=FALSE,row.names = TRUE)



r2 <- matrix(NA, ncol = 2, nrow = 10)  # Pre-allocate matrix
for (i in 1:10) {
  result <- r.squaredGLMM(best[[i]])
  if ("delta" %in% rownames(result)) {
    r2[i, ] <- result["delta", ]
  } else {
    warning(paste("Model", i, "does not include 'delta'."))
  }
}

r2
write.table(r2,"tables/v2_JAE/r2_dredge_stops_poisson_sex.csv",quote=FALSE,row.names = TRUE)


#######################
###---------duartion
#########################
### filter some missclassifications


l.stage<-subset(stop.segs, stop.segs$segment_type=="long.stage")

l.stage<-subset(l.stage,l.stage$bout.dur<200)

l.stage <- l.stage %>%
  group_by(dev, type2, season2, sex) %>%
  dplyr::summarise(
    mean.bout.dur = mean(bout.dur, na.rm = TRUE),
    year = first(year),   # Assuming year is consistent for each dev, type2, season2
    colony = first(colony),  # Assuming colony is consistent for each dev, type2, season2
    .groups = "drop"
  )

### filter sex data

table(l.stage.sex$sex)

l.stage.sex<-subset(l.stage, sex %in% c ("Male", "Female"))
l.stage.sex$type2 <- as.numeric(gsub("yo", "", l.stage.sex$type2))  # Remove 'yo' and convert to numeric
l.stage.sex <- subset(l.stage.sex, type2 %in% c("1", "2", "3", "4"))

unique_individuals <- l.stage.sex%>%
  group_by(type2, sex) %>%
  dplyr::summarise(unique_count = n_distinct(dev))

# #extract means for each season
# means <- aggregate(l.stage$bout.dur,by=list(l.stage$type2,l.stage$season2),mean)
# means <- aggregate(l.stage.for.mod$n.segments,by=list(l.stage.for.mod$type2,l.stage2$season2),mean)
# 
# winter <- subset(means,Group.2=="Autumn-Winter")
# summer <- subset(means,Group.2=="Spring-Summer")
# 
# #Now same for sd's
# sds <- aggregate(l.stage2$n.segments,by=list(l.stage2$age_gp,l.stage2$season2),sd)
# sdwinter <- subset(sds,Group.2=="Autumn-Winter")
# sdsummer <- subset(sds,Group.2=="Spring-Summer")

# Assuming your dataset is named 'data' and the age column is named 'age'

# Calculate the mean of age
mean_age <- mean(l.stage$type2, na.rm = TRUE)

# Calculate the standard deviation of age
sd_age <- sd(l.stage$type2, na.rm = TRUE)

# Print the results
cat("Mean of age:", mean_age, "\n")
cat("Standard deviation of age:", sd_age, "\n")



l.stage.sex$type2_scaled <- scale(l.stage.sex$type2)


# Calculate the standard deviation of the original variable
original_sd <- sd(l.stage$type2, na.rm = TRUE)

# The standard deviation of the scaled variable
scaled_sd <- sd(l.stage$type2_scaled, na.rm = TRUE)


### #re-scale variables

duration.age.scale <- attr(l.stage.sex$type2_scaled,"scaled:scale")
duration.age.center <- attr(l.stage.sex$type2_scaled,"scaled:center")


func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model

# # For the models
# functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)", 
#                                  "expm1(type2_scaled)")  # Add the categorical age model
# 
# 
# func_rel_names <- list("Linear", "Quadratic", "Exponential")
# 


###-------- duration ----
duration_models <- list()
duration_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("mean.bout.dur ~ season2 +sex* ", functional_relationships[[i]], 
                         " + (1|dev) + (1|year)++ (1|colony)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                              # for model comparison
                 #family = poisson,
                 data = l.stage.sex)
  
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

write.table(duration_comparison, "./tables/v2_JAE/AIC_duration_3dJAE_sex.csv")

## Figure
duration_mod_df <- ggpredict(duration_models[[2]], terms = c("type2_scaled[all]", "season2", "sex"))
duration_mod_df$type2 <- round(duration_mod_df$x * duration.age.scale + duration.age.center)
#duration_mod_df <- rename(duration_mod_df, season2 = group)



####plot
p2 <- ggplot(data = duration_mod_df) +
  facet_wrap(~sex) +

  stat_summary(
    data = duration_mod_df,
    aes(x = type2, y = predicted, color = group),
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "pointrange",
    size = 0.5,
    position = position_dodge(0.3)
  ) +
  geom_errorbar(data = duration_mod_df, aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
                width =0.1, position = position_dodge(width = 0.3), size=1) +
  #facet_grid(.~group)+
  ylim(0,150)+
  geom_jitter(data=l.stage.sex, aes(y = mean.bout.dur, x = type2, colour=season2, group=type2), size=1,alpha=.1,width = 0.1) +
  #Custom colors for the groups
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  labs(
    title = "",  # Add your title here
    x = "", 
    y = "Staging duration\n (days)"
  ) +
  
  # Y-axis label
  #scale_y_continuous(name = "Staging duration\n (days)") +
  #xlab("Age (years)") +  # Updated x-axis label for clarity
  # Customizing the theme
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        axis.line = element_line(colour = "black"),
        legend.position   = "none",
        #legend.direction  = 'vertical',
        legend.box.margin = margin(t=0.02,r=0,b=0,l=0,unit="cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  
  # Removing legend if you don't want it
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

p2




####
####estimates of teh linear model
# Fit the model as you've done
mfull = lmer(mean.bout.dur~ season2+sex*poly(type2_scaled, 2, raw = TRUE) + (1|dev)+(1|year),
             data = l.stage.sex)





tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)
summary(mfull)

## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)
# Run the outlier test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outlier test results
print(outlier_test) ### not extreme aoutliers


testDispersion(mfull) #1.00, ok
overdisp_fun(mfull) #0.87, ok


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




### select the best model using dredge

d1 <- dredge(mfull, rank = "AIC") # Do not include REML here

best <- get.models(d1,subset=NA) #get full model set
Anova(mfull) # anova for full model
Anova(best[[1]]) # run anova on best model from set
selection.table <- model.sel(d1) # pull out model selection table
selection.table #
write.table(selection.table,"tables/v2_JAE/dredge_duration_logv2_sex.csv",quote=FALSE,row.names = TRUE)


r2 <- matrix(NA, ncol = 2, nrow = 10)  # Pre-allocate matrix
colnames(r2) <- c("R2m", "R2c")      # Marginal and Conditional R²
for (i in 1:10) {
  result <- r.squaredGLMM(best[[i]])
  if (!is.null(result)) {
    r2[i, ] <- result[1, ]  # Assuming the first row contains the R² values
  } else {
    warning(paste("Model", i, "did not return R² values."))
  }
}

write.table(r2,"tables/v2_JAE/r2_dredge_duration_log_sex.csv",quote=FALSE,row.names = TRUE)


###########################################
### -------distance to consecutive sites
###############################

data2<-read.csv("./Data/long.stage_REV6_JAE.csv", header = TRUE, dec = ".", sep=",") #second attempt

#data2$distance_to_previous_km <- ifelse(is.na(data2$distance_to_previous_km) == TRUE,0,data2$distance_to_previous_km)

#remove nas
data2_clean <- data2 %>%
  filter(!is.na(distance_to_previous_km) & !is.na(prev_year_lat) & !is.na(prev_year_lon) & !is.na(prev_year))


#data2$type2 <- as.numeric(gsub("yo", "", data2$type2))  # Remove 'yo' and convert to numeric
data2_clean$type2_scaled <- scale(data2_clean$type2_num)
data2_clean$distance_scaled <- scale(data2_clean$distance_to_previous_km)

data2_clean<-subset(data2_clean, sex %in% c ("Male", "Female"))
data2_clean$type2 <- as.numeric(gsub("yo", "", data2_clean$type2))  # Remove 'yo' and convert to numeric
data2_clean <- subset(data2_clean, type2 %in% c("1", "2", "3", "4"))



data2<-data2_clean

unique_individuals <- data2%>%
  group_by(type2, sex) %>%
  dplyr::summarise(unique_count = n_distinct(dev))


data2$type2_scaled <- scale(data2$type2)



### #re-scale variables

duration.age.scale <- attr(l.stage.sex$type2_scaled,"scaled:scale")
duration.age.center <- attr(l.stage.sex$type2_scaled,"scaled:center")


### #re-scale variables

# Scaling for distance
distance.age.scale <- attr(data2$type2_scaled,"scaled:scale")
distance.age.center <- attr(data2$type2_scaled,"scaled:center")


func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model


###-------- distance ----
library(MuMIn)    # for calculating R²
distance_models <- list()
distance_mod_summary <- list()


##----------all data

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("distance_to_previous_km ~ season2+sex * ", functional_relationships[[i]], 
                         " + (1|dev) + (1|year)", sep ="")
  
  mod <- glmmTMB(formula = as.formula(model_formula), 
                 REML = F,                              # for model comparison
                 #family = poisson,
                 family = Gamma(log),
                 data = data2)
  
  distance_models[[i]] <- mod
  distance_mod_summary[[i]] <- summary(mod)
  
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

write.table(distance_comparison, "./tables/v2_JAE/AIC_distance_sex_dJAE.csv")

# List of model names
model_names <- c("Model 1: Exponential", "Model 2: Linear", "Model 3: Quadratic", "Model 4: poly")

# Compute R² for each model
r2_values <- lapply(distance_models, r.squaredGLMM)

# Combine the results into a data frame
r2_df <- do.call(rbind, r2_values)
colnames(r2_df) <- c("Marginal R²", "Conditional R²")

# Add model names as a column
r2_df <- cbind(Model = model_names, r2_df)

# Print the result
print(r2_df)

write.csv(distance_comparison, "./tables/v2_JAE/AIC_distance_log_JAE.csv")

## Figure
distance_mod_df <- ggpredict(distance_models[[1]], terms = c("type2_scaled[all]","season2","sex"))
distance_mod_df$type2 <- round(distance_mod_df$x * distance.age.scale + distance.age.center)


####plot
ptest <- ggplot(data = distance_mod_df) +
  stat_summary(
    data = distance_mod_df,
    aes(x = type2, y = predicted, color = group),
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "pointrange",
    size = 0.5,
    position = position_dodge(0.3)
  ) +
  # geom_errorbar(data = distance_mod_df, aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
  #               width =0.1, position = position_dodge(width = 0.3), size=1) +
  facet_grid(.~sex)+
  #ylim(0,1000)+
  geom_jitter(data=data2, aes(y = distance_to_previous_km, x = type2, colour=season2, group=type2), size=1,alpha=.1,width = 0.1) +
  #Custom colors for the groups
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  labs(
    title = "",  # Add your title here
    x = "", 
    y = "Distance to\n consecutive sites (km)"
  ) +
  
  # Y-axis label
  #scale_y_continuous(name = "Staging duration\n (days)") +
  #xlab("Age (years)") +  # Updated x-axis label for clarity
  # Customizing the theme
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        axis.line = element_line(colour = "black"),
        legend.position   = "none",
        #legend.direction  = 'vertical',
        legend.box.margin = margin(t=0.02,r=0,b=0,l=0,unit="cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  
  # Removing legend if you don't want it
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

ptest



options(scipen = 999)  # Suppresses scientific notation

####estimates of teh linear model
# Fit the model as you've done
# mfull = lmer(log(bee_line_distance_km) ~ sex*type2_scaled+season2 + (1|dev) + (1|year) + (1|colony),
#              data = data2)

mfull = glmmTMB(distance_to_previous_km ~ sex*type2_scaled + season2 + (1|dev) + (1|year),
                #family=nbinom2,
                family = Gamma(log),
                data = data2)



# poly(type2_scaled, 2, raw = TRUE)",
#                                  "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")
# 
tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)
summary(mfull)

## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)
# Run the outlier test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outlier test results
print(outlier_test) ### not extreme aoutliers


testDispersion(mfull) #0.7487, ok
overdisp_fun(mfull) #0.65, ok


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




### select the best model using dredge

d1 <- dredge(mfull, rank = "AIC") # Do not include REML here

best <- get.models(d1,subset=NA) #get full model set
Anova(mfull) # anova for full model
Anova(best[[1]]) # run anova on best model from set
selection.table <- model.sel(d1) # pull out model selection table
selection.table #
write.table(selection.table,"tables/v2_JAE/dredge_distance_log.csv",quote=FALSE,row.names = TRUE)


r2 <- matrix(NA, ncol = 2, nrow = 5)  # Pre-allocate matrix
#colnames(r2) <- c("R2m", "R2c")      # Marginal and Conditional R²
for (i in 1:5) {
  result <- r.squaredGLMM(best[[i]])
  if (!is.null(result)) {
    r2[i, ] <- result[1, ]  # Assuming the first row contains the R² values
  } else {
    warning(paste("Model", i, "did not return R² values."))
  }
}

write.table(r2,"tables/v2_JAE/r2_dredge_distance_log.csv",quote=FALSE,row.names = TRUE)

#####distance to the natal colony
data2<-read.csv("./Data/long.stage_REV6_JAE.csv", header = TRUE, dec = ".", sep=",") #second attempt

library(performance)
data2 <- data2 %>%
  filter(!is.na(distance_to_previous_km) & !is.na(prev_year_lat) & !is.na(prev_year_lon) & !is.na(prev_year))


data2_clean<-subset(data2_clean, sex %in% c ("Male", "Female"))
data2_clean$type2 <- as.numeric(gsub("yo", "", data2_clean$type2))  # Remove 'yo' and convert to numeric
data2_clean <- subset(data2_clean, type2 %in% c("1", "2", "3", "4"))
data2<-data2_clean

#data2$type2 <- as.numeric(gsub("yo", "", data2$type2))  # Remove 'yo' and convert to numeric
# Calculate distance of bird to center point of colony
data2$dist.to.colony <- deg.dist(lat1=data2$colony.lat,long1=data2$colony.long,long2=data2$mean.long,lat2=data2$mean.lat)


table(data2$sex)
# str(dist.breed)
### filter some missclassifications
#data2<-subset(data2,data2$bout.dur<200)

### #re-scale variables

dist.col.age.scale <- attr(data2$type2_scaled,"scaled:scale")
dist.col.age.center <- attr(data2$type2_scaled,"scaled:center")


func_rel_names <- list("Linear", "Quadratic", "Exponential", "Third-degree polynomial")

functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)",
                                 "expm1(type2_scaled)","poly(type2_scaled, 3, raw = TRUE)")  # Add the categorical age model

# # For the models
# functional_relationships <- list("type2_scaled", "poly(type2_scaled, 2, raw = TRUE)", 
#                                  "expm1(type2_scaled)")  # Add the categorical age model
# 
# 
# func_rel_names <- list("Linear", "Quadratic", "Exponential")
# 


###-------- duration ----
dist.col_models <- list()
dist.col_mod_summary <- list()

for(i in 1:length(functional_relationships)){
  
  model_formula <- paste("dist.to.colony ~ season2 + sex* ", functional_relationships[[i]], 
                         " + (1|dev)", sep ="")
  
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

write.table(dist.col_comparison, "./tables/v2_JAE/AIC_dist.col_log_3dJAE.csv")


## extract rcon manually
# mfull = lmer(log(dist.to.colony) ~ poly(type2_scaled, 3, raw = TRUE)*season2 + (1|dev) + (1|year) + (1|colony),
#              data = data2)
# 



## Figure
dist.col_mod_df <- ggpredict(dist.col_models[[1]], terms = c("type2_scaled[all]", "season2", "sex"))
dist.col_mod_df$type2 <- round(dist.col_mod_df$x * dist.col.age.scale + dist.col.age.center)
#duration_mod_df <- rename(duration_mod_df, season2 = group)
write.csv(dist.col_mod_df, "./tables/v2_JAE/dist.col_mod_JAE.csv")

dist.col.age.scale <- attr(data2$type2_scaled,"scaled:scale")
dist.col.age.center <- attr(data2$type2_scaled,"scaled:center")


####plot
p4 <- ggplot(data = dist.col_mod_df) +
  #Remove facet_grid and plot both groups in the same plot
  # geom_ribbon(aes(x = type2, y = predicted,
  #                 ymin = conf.low, ymax = conf.high, fill = group),
  #             alpha = 0.2, linewidth = 0.6) +
  geom_line(aes(x = type2, y = predicted, color = group),
            linewidth = 0.7) +
  
  stat_summary(
    data = dist.col_mod_df,
    aes(x = type2, y = predicted, color = group),
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "pointrange",
    size = 0.5,
    position = position_dodge(0.3)
  ) +
  geom_errorbar(data = dist.col_mod_df, aes(ymin = conf.low, ymax = conf.high, x = type2, colour = group),
                width =0.1, position = position_dodge(width = 0.3), size=1) +
  #facet_grid(.~group)+
  ylim(0,1500)+
  geom_jitter(data=data2, aes(y = dist.to.colony, x = type2_num, colour=season2, group=type2_num), size=1,alpha=.1,width = 0.1) +
  #Custom colors for the groups
  scale_color_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  scale_fill_manual(values = c("Autumn-Winter" = "#c7b42e", "Spring-Summer" = "#059748")) +
  labs(
    title = "",  # Add your title here
    x = "Age", 
    y = "Dist. of principal sites\n to natal colony (km)"
  ) +
  
  # Y-axis label
  #scale_y_continuous(name = "Staging duration\n (days)") +
  #xlab("Age (years)") +  # Updated x-axis label for clarity
  # Customizing the theme
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 13),
        axis.text = element_text(size = 13),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        axis.line = element_line(colour = "black"),
        legend.position   = "none",
        #legend.direction  = 'vertical',
        legend.box.margin = margin(t=0.02,r=0,b=0,l=0,unit="cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  
  # Removing legend if you don't want it
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

p4
ptest

####

####estimates of teh linear model
# Fit the model as you've done
# mfull = lmer(dist.to.colony ~  type2_scaled*season2 + (1|dev),
#              data = data2)

#mfull = lmer(log(distance_to_previous_km) ~ poly(type2_scaled, 3, raw = TRUE)*season2 + (1|dev) + (1|year),
#            data = data2)

# Fit the model, NA values in the data will be automatically excluded
mfull = glmmTMB(dist.to.colony ~ type2_scaled *sex+season2 + (1|dev),
                #family=nbinom2,
                family = Gamma(log),
                data = data2)

# mfull = glmmTMB(dist.to.colony ~ type2_scaled + season2 + (1|dev),
#                 #family=nbinom2,
#                 family = Gamma(log),
#                 data = data2)


#summ(mfull,ddf = "Satterthwaite", digits=3)
tab_model(mfull, show.se = TRUE, show.stat = TRUE, transform = NULL)
tab_model(mfull, transform = NULL, auto.label = FALSE, collapse.ci=TRUE, file = "./tables/results.xls")
summary(mfull)

## check model assumtions
tmp = simulateResiduals(mfull)
plot(tmp)
# Run the outlier test using the bootstrap method
outlier_test <- testOutliers(tmp, type = "bootstrap")
# Print the outlier test results
print(outlier_test) ### not extreme aoutliers


testDispersion(mfull) #1.00, ok
overdisp_fun(mfull) #0.87, ok


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




### select the best model using dredge

d1 <- dredge(mfull, rank = "AIC") # Do not include REML here

best <- get.models(d1,subset=NA) #get full model set
Anova(mfull) # anova for full model
Anova(best[[1]]) # run anova on best model from set
selection.table <- model.sel(d1) # pull out model selection table
selection.table #
write.table(selection.table,"tables/v2_JAE/dredge_dist.col_logv2.csv",quote=FALSE,row.names = TRUE)


r2 <- matrix(NA, ncol = 2, nrow = 5)  # Pre-allocate matrix
colnames(r2) <- c("R2m", "R2c")      # Marginal and Conditional R²
for (i in 1:5) {
  result <- r.squaredGLMM(best[[i]])
  if (!is.null(result)) {
    r2[i, ] <- result[1, ]  # Assuming the first row contains the R² values
  } else {
    warning(paste("Model", i, "did not return R² values."))
  }
}

write.table(r2,"tables/v2_JAE/r2_dredge_dist.col_log.csv",quote=FALSE,row.names = TRUE)

##################################
####individual differences
##################################

l.stage.for.mod<-read.csv('./Data/Num_res_events_for.mod.csv', header = TRUE, dec = ".", sep=",") #I create this file from the 01_fig_1_tracking_data_JAE

### 25 ids selected randomly
l.stage.for.mod_filt<- l.stage.for.mod %>% plotly::filter(dev %in% sample(unique(dev),25))

p3 <- ggplot(data = l.stage.for.mod_filt, aes(y = n.stops, x = type2, group = dev, color = dev)) +
  geom_blank() +
  ylim(0, 15) +
  scale_fill_viridis_c(option = "viridis", na.value = "white") +
  # Apply a linear smoothing method
  geom_smooth(method = 'lm', se = FALSE, fullrange = FALSE) +
  labs(
    title = "",
    x = "", 
    y = "Number of \nresident events"
  ) +
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_blank(),
    panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 13),
    axis.text = element_text(size = 13),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 13, face = 'bold'),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    legend.box.margin = margin(t = 0.02, r = 0, b = 0, l = 0, unit = "cm"),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  # Removing legend if you don't want it
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

p3



l.stage<-read.csv('./Data/Av_res_dur.csv', header = TRUE, dec = ".", sep=".") #I create this file from the 01_fig_1_tracking_data_JAE

l.stage_filt<- l.stage %>% plotly::filter(dev %in% sample(unique(dev),25))


p4<- ggplot(data=l.stage_filt,aes(y=mean.bout.dur,x=type2, group=dev, color = dev))+
  geom_blank()+
  ylim(0,190)+
  
  scale_fill_viridis_c(option="viridis",na.value="white")+
  #geom_smooth(method='loess', se = FALSE)+
  #geom_smooth(method=lm,se=FALSE,fullrange=TRUE,
  #aes(color=ID))+
  geom_smooth(method = 'lm', se = FALSE, fullrange = FALSE) +
  
  ylab("Average resident\n duration (days)") + xlab("")+ 
  theme_bw()+
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        axis.line = element_line(colour = "black"),
        legend.position   = "none",
        #legend.direction  = 'vertical',
        legend.box.margin = margin(t=0.02,r=0,b=0,l=0,unit="cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  
  # Removing legend if you don't want it
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

p4



data2<-read.csv('./Data/dist_prin_sites.csv', header = TRUE, dec = ".", sep=".") #I create this file from the 01_fig_1_tracking_data_JAE
data2_filt<- data2_clean %>% plotly::filter(dev %in% sample(unique(dev),25))



p5 <- ggplot(data = data2_filt, aes(y = distance_to_previous_km, x = type2, group = dev, color = dev)) +
  geom_blank() +
  ylim(0, 1500) +
  xlim(2,8)+
  scale_fill_viridis_c(option = "viridis", na.value = "white") +
  # Apply a linear smoothing method
  geom_smooth(method = 'lm', se = FALSE, fullrange = FALSE) +
  labs(
    title = "",  # Add your title here
    x = "Age", 
    y = "Dist. from previous\n principal sites (km)"
  ) +
  
  theme_bw() +
 
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 13),
        axis.text = element_text(size = 13),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        axis.line = element_line(colour = "black"),
        legend.position   = "none",
        #legend.direction  = 'vertical',
        legend.box.margin = margin(t=0.02,r=0,b=0,l=0,unit="cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) 

p5


data2<-read.csv("./Data/long.stage_REV6_JAE.csv", header = TRUE, dec = ".", sep=",") #second attempt

library(performance)
data2 <- data2 %>%
  filter(!is.na(distance_to_previous_km) & !is.na(prev_year_lat) & !is.na(prev_year_lon) & !is.na(prev_year))

# Calculate distance of bird to center point of colony
data2$dist.to.colony <- deg.dist(lat1=data2$colony.lat,long1=data2$colony.long,long2=data2$mean.long,lat2=data2$mean.lat)

data2$type2 <- as.numeric(gsub("yo", "", data2$type2))  # Remove 'yo' and convert to numeric
data2$type2_scaled <- scale(data2$type2_num)
data2$dist.to.colony_scaled <- scale(data2$dist.to.colony)

data2_filt<- data2%>% plotly::filter(dev %in% sample(unique(dev),25))


p6 <- ggplot(data = data2_filt, aes(y = distance_to_previous_km, x = type2, group = dev, color = dev)) +
  geom_blank() +
  ylim(0, 1500) +
  xlim(2,8)+
  scale_fill_viridis_c(option = "viridis", na.value = "white") +
  # Apply a linear smoothing method
  geom_smooth(method = 'lm', se = FALSE, fullrange = FALSE) +
  labs(
    title = "",  # Add your title here
    x = "Age", 
    y = "Dist. of principal sites\n to natal colony (km)"
  ) +
  
  theme_bw() +
  theme_bw(base_size = 14) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = alpha("#FFFFFF", 0.2)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(size = 13),
        axis.text = element_text(size = 13),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = 13, face = 'bold'),
        axis.line = element_line(colour = "black"),
        legend.position   = "none",
        #legend.direction  = 'vertical',
        legend.box.margin = margin(t=0.02,r=0,b=0,l=0,unit="cm"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  # Removing legend if you don't want it
  guides(alpha = FALSE, size = guide_legend(nrow = 2))

p6



# Create the top row with four panels
top_row <- cowplot::plot_grid(p3, p4,p6,p5, labels = c("A", "B", "C","D"),
                              label_size = 16, ncol = 2, align = "vh", axis = "lr")

top_row


ggsave(plot=top_row,filename='./Figures_2025/Fig_sup_mat_JAE.tiff',dpi=300,width=8.5,height=6)
