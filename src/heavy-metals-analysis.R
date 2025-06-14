##############################################################################
# R Script for Heavy metals Concentration Analysis
# Purpose: Perform Exploratory Data Analysis (EDA) and apply Generalized Linear Models by robust method (GLMs) or Hurdle Models (sdmTMB) to investigate the relationship between heavy metals concentrations in larval food of solitary bees between samples from yellow passion fruit crop (Crop) and Brazilian savanna fragments (Savanna).

######################## Elaborate by Thayane N. Araújo (Araujo TN et.al. 2025)
#########################################   - thayane.n.a@gmail.com -

##############################################################################
######################## Load Necessary Libraries ############################

library(vegan)
library(robustbase)
library(car)
library(sdmTMB)
library(glmmTMB)
library(dplyr)
library(DHARMa)
library(ggplot2)
library(ggpubr)
library(gridExtra)

#############################################################################
########################## Data Preparation #################################

data <-read.csv2("data_metals.csv", 
                 dec = ",", 
                 header =  T)

data$group <- as.factor(data$group)
data$sample <- as.factor(data$sample)
data$year <- as.factor(data$year)

summary(data)


###############################################################################
############### Permutation Multivariate Analysis of Variance #################


## Extracting the multivariate matrix
multi <- data[,6:ncol(data)]

## PERMANOVA analysis
set.seed(123) # For reproducible permutations
resu.adonis= adonis2(multi ~ data$group,
                     permutations = 999,
                     na.rm = TRUE, 
                     by= NULL)

## Displays the PERMANOVA output table 
print(resu.adonis)


################### Individually heavy metals concentration ##################
################################# Aluminum (Al) ##############################

## Exploratory Data Analysis (EDA) for Al

# Summary statistics
summary(data$Al)

# Histogram of Al
ggplot(data, aes(x = Al)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Al") +
  theme_minimal()

# Density Plot of Al
ggplot(data, aes(x = Al)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  ggtitle("Density Plot of Al") +
  theme_minimal()

# Histogram of log-transformed Al
ggplot(data, aes(x = log(Al))) +
  geom_histogram(binwidth = 0.2, fill = "orange", color = "black") +
  ggtitle("Histogram of log(Al)") +
  theme_minimal()

##### 

# -Generalized Linear Model by robust method (glmrob) for Al using  Gamma family-

# Robust Gamma Model
model_Al_robGamma <- glmrob(Al~group, 
                            family = Gamma(link = "log"), 
                            data = data)

# Check for model convergence. 'TRUE' ensures reliable parameter estimates (coefficients, standard errors and p-values). 'FALSE' means the optimization failed, making results unreliable.
model_Al_robGamma$converged

# Function to plot glmrob residuals
# This function is designed to be general for any glmrob model
plot_glmrob_diagnostics <- function(model_Al_robGamma) {
  # Extract deviance residuals
  dev_residuals <- residuals(model_Al_robGamma, type = "deviance")
  
  # Extract fitted (predicted) values
  fitted_values <- fitted(model_Al_robGamma)
  
  # Extract response variable name, model family, and link function for titles
  response_var_name <- as.character(formula(model_Al_robGamma))[2]
  model_family_name <- model_Al_robGamma$family$family
  model_link_name <- model_Al_robGamma$family$link
  
  # --- 1. Deviance Residuals vs. Fitted Values Plot (ggplot) ---
  df_residuals_fitted <- data.frame(Fitted = fitted_values, Deviance = dev_residuals)
  
  p1 <- ggplot(df_residuals_fitted, aes(x = Fitted, y = Deviance)) +
    geom_point(alpha = 0.7, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(title = paste0("Deviance Residuals vs. Fitted Values for ",response_var_name, model_family_name, "(link=", model_link_name, ")"),
         x = paste0("Fitted Values (link=", model_link_name, ")"),
         y = "Deviance Residuals") +
    theme_minimal()
  
  # --- 2. Q-Q Plot of Deviance Residuals (ggplot) ---
  df_qq_deviance <- data.frame(Deviance = dev_residuals)
  
  p2 <- ggplot(df_qq_deviance, aes(sample = Deviance)) +
    stat_qq(color = "darkblue", alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(title = "Q-Q Plot of Deviance Residuals",
         x = "Theoretical Normal Quantiles",
         y = "Deviance Residuals Quantiles") +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
}


# Analyze the residuals of the glmrob model
plot_glmrob_diagnostics(model_Al_robGamma)

# Print dispersion parameter for Gamma family
summary(model_Al_robGamma)

# Print ANOVA (Wald test) result
Anova(model_Al_robGamma, 
      type="II", 
      test.statistic = "Wald")



############################## Cadmium (Cd) ##################################

## Exploratory Data Analysis (EDA) for Cd

# Summary statistics
summary(data$Cd)

# Calculating the percentage of zero values
print(paste0("The percentage of zero values is: ", 
             round(sum(data$Cd == 0) / length(data$Cd) * 100, 2), 
             "%"))

# Histogram of Cd (only values > 0)
ggplot(subset(data, Cd > 0), aes(x = Cd)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Cd (only values > 0)") +
  theme_minimal()

# Density Plot of Cd (only values > 0)
ggplot(subset(data, Cd > 0), aes(x = Cd)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  ggtitle("Density Plot of Cd (only values > 0)") +
  theme_minimal()

# Histogram of log-transformed Cd (only values > 0)
ggplot(subset(data, Cd > 0), aes(x = log(Cd))) +
  geom_histogram(binwidth = 0.2, fill = "orange", color = "black") +
  ggtitle("Histogram of log(Cd) (only values > 0)") +
  theme_minimal()

##### 

# -Generalized Linear Model by robust method (glmrob) for Cd using Gamma family-

# Adding a constant value of 0.01 to enable the use of the Gamma distribution, which requires strictly positive values
data$Cd_constant <- data$Cd + 0.01 

# Robust Gamma Model with Cd + 0.01 constant (Cd_constant)
model_Cd_robGamma <- glmrob(Cd_constant~group, 
                            family = Gamma(link = "log"), 
                            data = data)

# Check for model convergence. 'TRUE' ensures reliable parameter estimates. 'FALSE' means the optimization failed, making results unreliable.
model_Cd_robGamma$converged


# -Fitting a Hurdle GLMM (sdmTMB) using the delta_gamma family (non-spatial, non-spatiotemporal)-

# Fits the sdmTMB model, which handles zero-inflated continuous data using a delta_gamma distribution
model_hurdleCd_deltagamma <- sdmTMB(Cd ~ group,
                             family = delta_gamma(),
                             spatial = "off", 
                             spatiotemporal = "off",
                             data = data)

# Plots a histogram of the raw residuals from the sdmTMB model
hist(residuals(model_hurdleCd_deltagamma))

# Generates 250 sets of simulated responses from the fitted sdmTMB model for residual diagnostics
simulated_response <- simulate(model_hurdleCd_deltagamma, nsim = 250)

# Creates a DHARMa object for comprehensive residual diagnostics
res_hurdle_gamma <- createDHARMa(
   simulatedResponse = simulated_response,
   observedResponse = data$Cd,
)

# Analyze the residuals of delta_gamma model using  DHARMa
plot(res_hurdle_gamma)

# Displays a detailed summary of the sdmTMB model results, including coefficients and convergence information
summary(model_hurdleCd_deltagamma)

# Getting the tidy() result
# For the Binomial component (presence or absence):
tidy_binomial <- tidy(model_hurdleCd_deltagamma, model = 1)

print(tidy_binomial)


results_binomial <- data.frame(
  # Data from tidy_binomial output
  Term = c("(Intercept)", "groupSavanna"),
  Estimate = c(0.405, -0.405),
  Std.Error = c(0.913, 1.08)
) %>%
  dplyr::mutate(
    # Calculate the Z-value
    Z_Value = Estimate / Std.Error,
    # Calculate the two-tailed p-value
    P_Value = 2 * pnorm(-abs(Z_Value))
  )

# Print results for the Binomial Component (Presence/Absence)
print(results_binomial)



# For the Gamma component (positive values):
tidy_gamma <- tidy(model_hurdleCd_deltagamma, model = 2)
print(tidy_gamma)



results_gamma <- data.frame(
  Term = c("(Intercept)", "groupSavanna"),
  Estimate = c(1.93, 0.559),
  Std.Error = c(0.568, 0.695)
) %>%
  dplyr::mutate(
    Z_Value = Estimate / Std.Error,
    P_Value = 2 * pnorm(-abs(Z_Value))
  )

# Print results for the Gamma Component (Positive Values)
print(results_gamma)


########################### Chromium (Cr) ###################################

## Exploratory Data Analysis (EDA) for Cr

# Summary statistics
summary(data$Cr)

# Calculating the percentage of zero
print(paste0("The percentage of zero values is: ", 
             round(sum(data$Cr == 0) / length(data$Cr) * 100, 2), 
             "%"))

# Histogram of Cr 
ggplot(data, aes(x = Cr)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Cr") +
  theme_minimal()

# Density Plot of Cr
ggplot(data, aes(x = Cr)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  ggtitle("Density Plot of Cr") +
  theme_minimal()

# Histogram of log-transformed Cr
ggplot(data, aes(x = log(Cr))) +
  geom_histogram(binwidth = 0.2, fill = "orange", color = "black") +
  ggtitle("Histogram of log(Cr)") +
  theme_minimal()



# -Generalized Linear Model by robust method (glmrob) for Cr with gaussian family-

# Robust gaussian Model 
model_Cr_robGAU <- glmrob(Cr~group, 
                            family = gaussian(), 
                            data = data)


# Check for model convergence. 'TRUE' ensures reliable parameter estimates. 'FALSE' means the optimization failed, making results unreliable.
model_Cr_robGAU$converged


# Function to plot GLM residuals
plot_glmrob_diagnostics <- function(model_Cr_robGAU) {
  dev_residuals <- residuals(model_Cr_robGAU, type = "deviance")
  fitted_values <- fitted(model_Cr_robGAU)
  
  response_var_name <- as.character(formula(model_Cr_robGAU))[2]
  model_family_name <- model_Cr_robGAU$family$family
  model_link_name <- model_Cr_robGAU$family$link
  
  # --- 1. Deviance Residuals vs. Fitted Values Plot (ggplot) ---
  df_residuals_fitted <- data.frame(Fitted = fitted_values, Deviance = dev_residuals)
  
  p1 <- ggplot(df_residuals_fitted, aes(x = Fitted, y = Deviance)) +
    geom_point(alpha = 0.7, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(title = paste0("Deviance Residuals vs. Fitted Values for ",response_var_name, model_family_name, "(link=", model_link_name, ")"),
         x = paste0("Fitted Values (link=", model_link_name, ")"),
         y = "Deviance Residuals") +
    theme_minimal()
  
  # --- 2. Q-Q Plot of Deviance Residuals (ggplot) ---
  df_qq_deviance <- data.frame(Deviance = dev_residuals)
  
  p2 <- ggplot(df_qq_deviance, aes(sample = Deviance)) +
    stat_qq(color = "darkblue", alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(title = "Q-Q Plot of Deviance Residuals",
         x = "Theoretical Normal Quantiles",
         y = "Deviance Residuals Quantiles") +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
}


# Analyze the residuals of the glmrob model
plot_glmrob_diagnostics(model_Cr_robGAU)

# Print dispersion parameter for gaussian family
summary(model_Cr_robGAU)

# Print ANOVA (Wald test) result
Anova(model_Cr_robGAU, 
      type="II", 
      test.statistic = "Wald")




################################ Copper (Cu) ##################################

## Exploratory Data Analysis (EDA) for Cu

# Summary statistics
summary(data$Cu)

# Histogram of Cu
ggplot(data, aes(x = Cu)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Cu") +
  theme_minimal()

# Density Plot of Cu
ggplot(data, aes(x = Cu)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  ggtitle("Density Plot of Cu") +
  theme_minimal()

# Histogram of log-transformed Cu
ggplot(data, aes(x = log(Cu))) +
  geom_histogram(binwidth = 0.2, fill = "orange", color = "black") +
  ggtitle("Histogram of log(Cu)") +
  theme_minimal()



# -Generalized Linear Model by robust method (glmrob) for Cu with Gamma family-

# Robust Gamma Model
model_Cu_robGamma <- glmrob(Cu~group, 
                            family = Gamma(link = "log"), 
                            data = data)

# Check for model convergence. 'TRUE' ensures reliable parameter estimates. 'FALSE' means the optimization failed, making results unreliable.
model_Cu_robGamma$converged

# Function to plot glmrob residuals
plot_glmrob_diagnostics <- function(model_Cu_robGamma) {
  dev_residuals <- residuals(model_Cu_robGamma, type = "deviance")
  fitted_values <- fitted(model_Cu_robGamma)
  
  response_var_name <- as.character(formula(model_Cu_robGamma))[2]
  model_family_name <- model_Cu_robGamma$family$family
  model_link_name <- model_Cu_robGamma$family$link
  
  # --- 1. Deviance Residuals vs. Fitted Values Plot (ggplot) ---
  df_residuals_fitted <- data.frame(Fitted = fitted_values, Deviance = dev_residuals)
  
  p1 <- ggplot(df_residuals_fitted, aes(x = Fitted, y = Deviance)) +
    geom_point(alpha = 0.7, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(title = paste0("Deviance Residuals vs. Fitted Values for ",response_var_name, model_family_name, "(link=", model_link_name, ")"),
         x = paste0("Fitted Values (link=", model_link_name, ")"),
         y = "Deviance Residuals") +
    theme_minimal()
  
  # --- 2. Q-Q Plot of Deviance Residuals (ggplot) ---
  df_qq_deviance <- data.frame(Deviance = dev_residuals)
  
  p2 <- ggplot(df_qq_deviance, aes(sample = Deviance)) +
    stat_qq(color = "darkblue", alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(title = "Q-Q Plot of Deviance Residuals",
         x = "Theoretical Normal Quantiles",
         y = "Deviance Residuals Quantiles") +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
}


# Analyze the residuals of the glmrob model
plot_glmrob_diagnostics(model_Cu_robGamma)


# Print dispersion parameter for Gamma family
summary(model_Cu_robGamma)

# Print ANOVA (Wald test) result
Anova(model_Cu_robGamma, 
      type="II", 
      test.statistic = "Wald")



############################### Iron (Fe) #####################################

## Exploratory Data Analysis (EDA) for Fe

# Summary statistics
summary(data$Fe)

# Histogram of Fe
ggplot(data, aes(x = Fe)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Fe") +
  theme_minimal()

# Density Plot of Fe
ggplot(data, aes(x = Fe)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  ggtitle("Density Plot of Fe") +
  theme_minimal()

# Histogram of log-transformed Fe
ggplot(data, aes(x = log(Fe))) +
  geom_histogram(binwidth = 0.2, fill = "orange", color = "black") +
  ggtitle("Histogram of log(Fe)") +
  theme_minimal()



# -Generalized Linear Model by robust method (glmrob) for Fe with Gamma family-

# Robust Gamma Model
model_Fe_robGamma <- glmrob(Fe~group, 
                            family = Gamma(link = "log"), 
                            data = data)

# Check for model convergence. 'TRUE' ensures reliable parameter estimates. 'FALSE' means the optimization failed, making results unreliable.
model_Fe_robGamma$converged

# Function to plot glmrob residuals
plot_glmrob_diagnostics <- function(model_Fe_robGamma) {
  dev_residuals <- residuals(model_Fe_robGamma, type = "deviance")
  fitted_values <- fitted(model_Fe_robGamma)
  
  response_var_name <- as.character(formula(model_Fe_robGamma))[2]
  model_family_name <- model_Fe_robGamma$family$family
  model_link_name <- model_Fe_robGamma$family$link
  
  # --- 1. Deviance Residuals vs. Fitted Values Plot (ggplot) ---
  df_residuals_fitted <- data.frame(Fitted = fitted_values, Deviance = dev_residuals)
  
  p1 <- ggplot(df_residuals_fitted, aes(x = Fitted, y = Deviance)) +
    geom_point(alpha = 0.7, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(title = paste0("Deviance Residuals vs. Fitted Values for ",response_var_name, model_family_name, "(link=", model_link_name, ")"),
         x = paste0("Fitted Values (link=", model_link_name, ")"),
         y = "Deviance Residuals") +
    theme_minimal()
  
  # --- 2. Q-Q Plot of Deviance Residuals (ggplot) ---
  df_qq_deviance <- data.frame(Deviance = dev_residuals)
  
  p2 <- ggplot(df_qq_deviance, aes(sample = Deviance)) +
    stat_qq(color = "darkblue", alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(title = "Q-Q Plot of Deviance Residuals",
         x = "Theoretical Normal Quantiles",
         y = "Deviance Residuals Quantiles") +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
}


# Analyze the residuals of the glmrob model
plot_glmrob_diagnostics(model_Fe_robGamma)

# Print dispersion parameter for Gamma family
summary(model_Fe_robGamma)

# Print ANOVA (Wald test) result
Anova(model_Fe_robGamma, 
      type="II", 
      test.statistic = "Wald")



############################### Nickel  (Ni) #################################

## Exploratory Data Analysis (EDA) for Ni

# Summary statistics
summary(data$Ni)

# Calculating the percentage of zero
print(paste0("The percentage of zero values is: ", 
             round(sum(data$Ni == 0) / length(data$Ni) * 100, 2), 
             "%"))

# Histogram of Ni (only values > 0)
ggplot(subset(data, Ni > 0), aes(x = Ni)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Ni (only values > 0)") +
  theme_minimal()

# Density Plot of Ni (only values > 0)
ggplot(subset(data, Ni > 0), aes(x = Ni)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  ggtitle("Density Plot of Ni (only values > 0)") +
  theme_minimal()

# Histogram of log-transformed Ni (only values > 0)
ggplot(subset(data, Ni > 0), aes(x = log(Ni))) +
  geom_histogram(binwidth = 0.2, fill = "orange", color = "black") +
  ggtitle("Histogram of log(Ni) (only values > 0)") +
  theme_minimal()



# -Generalized Linear Model by robust method (glmrob) for Ni with Binomial family-

# Transform in presence and absence to allow use the Binomial distribution

data <- data %>%
  mutate(Ni_bi = ifelse(Ni > 0, 1, 0))


summary(data)


# GLM with binomial family
model_Ni_Binomial <- ?glm(Ni_bi~group, 
                            family = binomial(), 
                            data = data)

# Check for model convergence. 'TRUE' ensures reliable parameter estimates. 'FALSE' means the optimization failed, making results unreliable.
model_Ni_Binomial$converged

# Function to plot glmrob residuals
plot_glmrob_diagnostics <- function(model_Ni_Binomial) {
  dev_residuals <- residuals(model_Ni_Binomial, type = "deviance")
  fitted_values <- fitted(model_Ni_Binomial)
  
  response_var_name <- as.character(formula(model_Ni_Binomial))[2]
  model_family_name <- model_Ni_Binomial$family$family
  model_link_name <- model_Ni_Binomial$family$link
  
  # --- 1. Deviance Residuals vs. Fitted Values Plot (ggplot) ---
  df_residuals_fitted <- data.frame(Fitted = fitted_values, Deviance = dev_residuals)
  
  p1 <- ggplot(df_residuals_fitted, aes(x = Fitted, y = Deviance)) +
    geom_point(alpha = 0.7, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(title = paste0("Deviance Residuals vs. Fitted Values for ",response_var_name, model_family_name, "(link=", model_link_name, ")"),
         x = paste0("Fitted Values (link=", model_link_name, ")"),
         y = "Deviance Residuals") +
    theme_minimal()
  
  # --- 2. Q-Q Plot of Deviance Residuals (ggplot) ---
  df_qq_deviance <- data.frame(Deviance = dev_residuals)
  
  p2 <- ggplot(df_qq_deviance, aes(sample = Deviance)) +
    stat_qq(color = "darkblue", alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(title = "Q-Q Plot of Deviance Residuals",
         x = "Theoretical Normal Quantiles",
         y = "Deviance Residuals Quantiles") +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
}


# Analyze the residuals of the glmrob model
plot_glmrob_diagnostics(model_Ni_Binomial)

summary(model_Ni_Binomial)

# Print ANOVA (Wald test) result
Anova(model_Ni_Binomial, 
      type="II", 
      test.statistic = "Wald")

############################### Lead  (Pb) ###################################

## Exploratory Data Analysis (EDA) for Pb

# Summary statistics
summary(data$Pb)

# Calculating the percentage of zero
print(paste0("The percentage of zero values is: ", 
             round(sum(data$Pb == 0) / length(data$Pb) * 100, 2), 
             "%"))

# Histogram of Pb (only values > 0)
ggplot(subset(data, Pb > 0), aes(x = Pb)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Pb (only values > 0)") +
  theme_minimal()

# Density Plot of Pb (only values > 0)
ggplot(subset(data, Pb > 0), aes(x = Pb)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  ggtitle("Density Plot of Pb (only values > 0)") +
  theme_minimal()

# Histogram of log-transformed Pb (only values > 0)
ggplot(subset(data, Pb > 0), aes(x = log(Pb))) +
  geom_histogram(binwidth = 0.2, fill = "orange", color = "black") +
  ggtitle("Histogram of log(Pb) (only values > 0)") +
  theme_minimal()



# -Generalized Linear Model by robust method (glmrob) for Pb with Gamma family-

# Adding a constant value to allow use the Gamma distribution
data$Pb_constant <- data$Pb + 0.01 

# Fits the robust GLM with Pb_constant using the Gamma family with a log link function
model_Pb_robGamma <- glmrob(Pb_constant~group, 
                            family = Gamma(link = "log"), 
                            data = data)

# Check for model convergence. 'TRUE' ensures reliable parameter estimates. 'FALSE' means the optimization failed, making results unreliable.
model_Pb_robGamma$converged

# Function to plot glmrob residuals
plot_glmrob_diagnostics <- function(model_Pb_robGamma) {
  dev_residuals <- residuals(model_Pb_robGamma, type = "deviance")
  fitted_values <- fitted(model_Pb_robGamma)
  
  response_var_name <- as.character(formula(model_Pb_robGamma))[2]
  model_family_name <- model_Pb_robGamma$family$family
  model_link_name <- model_Pb_robGamma$family$link
  
  # --- 1. Deviance Residuals vs. Fitted Values Plot (ggplot) ---
  df_residuals_fitted <- data.frame(Fitted = fitted_values, Deviance = dev_residuals)
  
  p1 <- ggplot(df_residuals_fitted, aes(x = Fitted, y = Deviance)) +
    geom_point(alpha = 0.7, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(title = paste0("Deviance Residuals vs. Fitted Values for ",response_var_name, model_family_name, "(link=", model_link_name, ")"),
         x = paste0("Fitted Values (link=", model_link_name, ")"),
         y = "Deviance Residuals") +
    theme_minimal()
  
  # --- 2. Q-Q Plot of Deviance Residuals (ggplot) ---
  df_qq_deviance <- data.frame(Deviance = dev_residuals)
  
  p2 <- ggplot(df_qq_deviance, aes(sample = Deviance)) +
    stat_qq(color = "darkblue", alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(title = "Q-Q Plot of Deviance Residuals",
         x = "Theoretical Normal Quantiles",
         y = "Deviance Residuals Quantiles") +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
}


# Analyze the residuals of the glmrob model
plot_glmrob_diagnostics(model_Pb_robGamma)


# -Fitting a Hurdle GLMM (sdmTMB) using the delta_gamma family (non-spatial, non-spatiotemporal)-

# Fits the sdmTMB model, which handles zero-inflated continuous data using a delta_gamma distribution
model_hurdlePb_deltagamma <- sdmTMB(Pb ~ group,
                                    family = delta_gamma(),
                                    spatial = "off", 
                                    spatiotemporal = "off",
                                    data = data)

# Plots a histogram of the raw residuals from the sdmTMB model
hist(residuals(model_hurdlePb_deltagamma))

# Generates 250 sets of simulated responses from the fitted sdmTMB model for residual diagnostics
simulated_response <- simulate(model_hurdlePb_deltagamma, nsim = 250)

# Creates a DHARMa object for comprehensive residual diagnostics
res_hurdlePb_gamma <- createDHARMa(
  simulatedResponse = simulated_response,
  observedResponse = data$Pb,
)

# Analyze the residuals of delta_gamma model using  DHARMa
plot(res_hurdlePb_gamma)

# Displays a detailed summary of the sdmTMB model results, including coefficients and convergence information
summary(model_hurdlePb_deltagamma)

# Getting the tidy() result
# For the Binomial component (presence or absence):
tidy_binomial <- tidy(model_hurdlePb_deltagamma, model = 1)
print(tidy_binomial)


results_binomial <- data.frame(
  Term = c("(Intercept)", "groupSavanna"),
  Estimate = c(21.7, -21.4),
  Std.Error = c(23148., 23148.)
) %>%
  dplyr::mutate(
    Z_Value = Estimate / Std.Error,
    P_Value = 2 * pnorm(-abs(Z_Value))
  )

# Print results for the Binomial Component (Presence/Absence)
print(results_binomial)



# For the Gamma component (positive values):
tidy_gamma <- tidy(model_hurdlePb_deltagamma, model = 2)
print(tidy_gamma)


results_gamma <- data.frame(
  Term = c("(Intercept)", "groupSavanna"),
  Estimate = c(1.32, 0.812),
  Std.Error = c(0.455, 0.595)
) %>%
  dplyr::mutate(
    Z_Value = Estimate / Std.Error,
    P_Value = 2 * pnorm(-abs(Z_Value))
  )

# Print results for the Gamma Component (Positive Values)
print(results_gamma)




################################ Tin (Sn) #####################################

## Exploratory Data Analysis (EDA) for Sn

# Summary statistics
summary(data$Sn)

# Histogram of Sn
ggplot(data, aes(x = Sn)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Sn") +
  theme_minimal()

# Density Plot of Sn
ggplot(data, aes(x = Sn)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  ggtitle("Density Plot of Sn") +
  theme_minimal()

# Histogram of log-transformed Sn
ggplot(data, aes(x = log(Sn))) +
  geom_histogram(binwidth = 0.2, fill = "orange", color = "black") +
  ggtitle("Histogram of log(Sn)") +
  theme_minimal()



# -Generalized Linear Model by robust method (glmrob) for Sn with Gamma family-

# Fits the robust GLM using the Gamma family with a log link function
model_Sn_robGamma <- glmrob(Sn ~ group, 
                            family = Gamma(link = "log"), 
                            data = data)

# Check for model convergence. 'TRUE' ensures reliable parameter estimates. 'FALSE' means the optimization failed, making results unreliable.
model_Sn_robGamma$converged

# Function to plot glmrob residuals
plot_glmrob_diagnostics <- function(model_Sn_robGamma) {
  dev_residuals <- residuals(model_Sn_robGamma, type = "deviance")
  fitted_values <- fitted(model_Sn_robGamma)
  
  response_var_name <- as.character(formula(model_Sn_robGamma))[2]
  model_family_name <- model_Sn_robGamma$family$family
  model_link_name <- model_Sn_robGamma$family$link
  
  # --- 1. Deviance Residuals vs. Fitted Values Plot (ggplot) ---
  df_residuals_fitted <- data.frame(Fitted = fitted_values, Deviance = dev_residuals)
  
  p1 <- ggplot(df_residuals_fitted, aes(x = Fitted, y = Deviance)) +
    geom_point(alpha = 0.7, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(title = paste0("Deviance Residuals vs. Fitted Values for ",response_var_name, model_family_name, "(link=", model_link_name, ")"),
         x = paste0("Fitted Values (link=", model_link_name, ")"),
         y = "Deviance Residuals") +
    theme_minimal()
  
  # --- 2. Q-Q Plot of Deviance Residuals (ggplot) ---
  df_qq_deviance <- data.frame(Deviance = dev_residuals)
  
  p2 <- ggplot(df_qq_deviance, aes(sample = Deviance)) +
    stat_qq(color = "darkblue", alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(title = "Q-Q Plot of Deviance Residuals",
         x = "Theoretical Normal Quantiles",
         y = "Deviance Residuals Quantiles") +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
}


# Analyze the residuals of the glmrob model
plot_glmrob_diagnostics(model_Sn_robGamma)

# Print dispersion parameter for Gamma family
summary(model_Sn_robGamma)

# Print ANOVA (Wald test) result
Anova(model_Sn_robGamma, 
      type="II", 
      test.statistic = "Wald")




############################## Zinc (Zn) #####################################

## Exploratory Data Analysis (EDA) for Zn

# Summary statistics
summary(data$Zn)

# Histogram of Zn
ggplot(data,  aes(x = Zn)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Zn") +
  theme_minimal()

# Density Plot of Z
ggplot(data,  aes(x = Zn)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  ggtitle("Density Plot of Zn") +
  theme_minimal()

# Histogram of log-transformed Zn
ggplot(data, aes(x = log(Zn))) +
  geom_histogram(binwidth = 0.2, fill = "orange", color = "black") +
  ggtitle("Histogram of log(Zn)") +
  theme_minimal()



# -Generalized Linear Model by robust method (glmrob) for Zn with Gamma family-

# Fits the robust GLM using the Gamma family with a log link function
model_Zn_robGamma <- glmrob(Zn ~ group, 
                            family = Gamma(link = "log"), 
                            data = data)

# Check for model convergence. 'TRUE' ensures reliable parameter estimates. 'FALSE' means the optimization failed, making results unreliable.
model_Zn_robGamma$converged

# Function to plot glmrob residuals
plot_glmrob_diagnostics <- function(model_Zn_robGamma) {
  dev_residuals <- residuals(model_Zn_robGamma, type = "deviance")
  fitted_values <- fitted(model_Zn_robGamma)
  
  response_var_name <- as.character(formula(model_Zn_robGamma))[2]
  model_family_name <- model_Zn_robGamma$family$family
  model_link_name <- model_Zn_robGamma$family$link
  
  # --- 1. Deviance Residuals vs. Fitted Values Plot (ggplot) ---
  df_residuals_fitted <- data.frame(Fitted = fitted_values, Deviance = dev_residuals)
  
  p1 <- ggplot(df_residuals_fitted, aes(x = Fitted, y = Deviance)) +
    geom_point(alpha = 0.7, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(title = paste0("Deviance Residuals vs. Fitted Values for ",response_var_name, model_family_name, "(link=", model_link_name, ")"),
         x = paste0("Fitted Values (link=", model_link_name, ")"),
         y = "Deviance Residuals") +
    theme_minimal()
  
  # --- 2. Q-Q Plot of Deviance Residuals (ggplot) ---
  df_qq_deviance <- data.frame(Deviance = dev_residuals)
  
  p2 <- ggplot(df_qq_deviance, aes(sample = Deviance)) +
    stat_qq(color = "darkblue", alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(title = "Q-Q Plot of Deviance Residuals",
         x = "Theoretical Normal Quantiles",
         y = "Deviance Residuals Quantiles") +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
}


# Analyze the residuals of the glmrob model
plot_glmrob_diagnostics(model_Zn_robGamma)

# Print dispersion parameter for Gamma family
summary(model_Zn_robGamma)

# Print ANOVA (Wald test) result
Anova(model_Zn_robGamma, 
      type="II", 
      test.statistic = "Wald")


